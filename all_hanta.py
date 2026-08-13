#!/usr/bin/env python3
"""Orthohantavirus NCBI Datasets production pipeline.

Ingests an NCBI Datasets virus package for Orthohantavirus (taxon 1980442) and produces:
- accession-level records
- canonical species records (ICTV binomial, resolved from virus.lineage)
- isolate-level records reconstructed from S/M/L segments
- host-level records
- pathogenicity / surveillance-risk labels
- reassortment candidate calls
- ML-ready feature tables
- QC, suppression and provenance tables

Reproducibility model
---------------------
Integrity is established against NCBI's own published checksums, not against digests
this script invents:

- file level   : ``md5sum.txt`` shipped inside the package is parsed and verified.
- content level: ``nucleotide.sequenceHash`` (CRC32, emitted per record by NCBI) is
                 used for sequence identity and de-duplication. It is independent of
                 line endings and platform, so it behaves identically on Linux and
                 Windows even when a checkout has mangled the files.

``hashlib`` is still used, but only to hash *our own* outputs and configuration --
never to re-derive something NCBI already published.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import traceback
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from functools import lru_cache
from pathlib import Path
from typing import Any, DefaultDict, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------

VERSION = "2.0.0"

ORTHOHANTAVIRUS_TAXID = 1980442

DATA_REPORT = "data_report.jsonl"
ANNOTATION_REPORT = "annotation_report.jsonl"
BIOSAMPLE_REPORT = "biosample_report.jsonl"
GENOMIC_FASTA = "genomic.fna"
PROTEIN_FASTA = "protein.faa"
CDS_FASTA = "cds.fna"
MD5SUM_FILE = "md5sum.txt"

YEAR_RE = re.compile(r"(19|20)\d{2}")
# Hot-path patterns. clean_text alone runs ~930k times per invocation, so these are
# compiled once here rather than re-resolved from re's cache on every call.
WHITESPACE_RE = re.compile(r"\s+")
NON_ALNUM_RE = re.compile(r"[^a-z0-9]+")
HEX_HASH_RE = re.compile(r"[0-9A-F]{1,8}")
MD5_RE = re.compile(r"[0-9a-fA-F]{32}")
SEGMENT_STRIP_RE = re.compile(r"[^a-z0-9;\s()/-]+")
SEGMENT_PUNCT_RE = re.compile(r"[;()/_-]+")
ATTR_NAME_RE = re.compile(r"[^a-z0-9]+")

# ICTV Hantaviridae segment sizes: S 1.2-2.1 kb, M 3.4-3.8 kb, L 6.5-6.6 kb.
# https://ictv.global/report/chapter/hantaviridae/hantaviridae/mammantavirinae/orthohantavirus
#
# S and M are used as ICTV states them. L is widened to 6300-6800 rather than 6500-6600:
# submitted L segments cluster just outside the ICTV figure at both ends, and tightening to
# the literal bounds would drop a 6,462 nt record -- unambiguously an L at that length --
# back to 'unresolved'. It changes no length_anomaly flag either way.
SEGMENT_LENGTH_RANGE = {"S": (1200, 2100), "M": (3400, 3800), "L": (6300, 6800)}
# Minimum length at which a segment can plausibly be coding-complete. Slightly below the
# ICTV lower bound because several validated coding-complete records are trimmed to the ORF.
COMPLETE_LENGTH_MIN = {"S": 1100, "M": 3200, "L": 6000}
# Beyond these a record is almost certainly a multi-segment concatenation or a mis-label.
SEGMENT_LENGTH_ABSURD_MAX = {"S": 3000, "M": 5000, "L": 8000}
# How far below the ICTV lower bound a COMPLETE record may fall before it is flagged.
COMPLETE_LENGTH_TOLERANCE = 0.85

KMER_K = 3

# Ranks pulled out of virus.lineage, ordered from broad to narrow. This is the
# "until there is no bigger group" hierarchy -- it is already present in every record.
LINEAGE_RANKS = [
    ("realm", "Riboviria"),
    ("kingdom", "Orthornavirae"),
    ("phylum", "Negarnaviricota"),
    ("subphylum", "Polyploviricotina"),
    ("class", "Bunyaviricetes"),
    ("order", "Elliovirales"),
    ("family", "Hantaviridae"),
    ("subfamily", "Mammantavirinae"),
    ("genus", "Orthohantavirus"),
]

# --------------------------------------------------------------------------------------
# Segment vocabulary
# --------------------------------------------------------------------------------------

# Values actually observed in the ``segment`` field of this package, beyond plain S/M/L:
# 'Partial L', 'S segment', 'M-segment', 'M; medium', 'small (S)', 'medium (M)',
# 'G1', 'G2', 'ANONYMOUS', 's', 'N', 'Small', 'Medium', 'Large', 'middle'.
SEGMENT_FIELD_RULES: List[Tuple[re.Pattern, str]] = [
    (re.compile(r"\b(l|large)\b"), "L"),
    (re.compile(r"\b(m|medium|middle)\b"), "M"),
    (re.compile(r"\b(s|small)\b"), "S"),
    # G1/Gn and G2/Gc are both encoded on the M segment.
    (re.compile(r"\bg[12]\b|\bg[nc]\b"), "M"),
    # 'N' means nucleocapsid, which is the S segment product.
    (re.compile(r"^n$"), "S"),
]

# Derived from the 122 distinct CDS product names in annotation_report.jsonl.
#
# Two related but deliberately different questions are asked of a product name:
#   - which SEGMENT does it imply?  (loose: any nonstructural mention implies S)
#   - which PROTEIN ROLE is it?     (strict: only a true NSs product counts as NSs)
# The shared fragments are defined once here so the two rule tables cannot drift, while
# the intentional differences stay visible as explicit additions.
_RDRP = r"polymerase|\brdrp\b|\bl protein\b|^l$"
_GPC = r"glyco|\bgpc\b|\bg1\b|\bg2\b|\bgn\b|\bgc\b|envelope|polyprotein"
_N = r"nucleocapsid|nucleoprotein|\bn protein\b"
_NSS_STRICT = r"\bnss\b|small non-?structural|non-?structural protein s"

# Order matters: polymerase is checked before glycoprotein because a few deflines
# mention both.
PRODUCT_SEGMENT_RULES: List[Tuple[re.Pattern, str]] = [
    (re.compile(_RDRP + r"|\blarge protein\b"), "L"),
    (re.compile(_GPC + r"|m segment"), "M"),
    # Looser than the NSs role rule on purpose: a bare "nonstructural protein" pins the
    # record to S without asserting it is the NSs product.
    (re.compile(_N + r"|\bnss?\b|nonstructural|non-structural"), "S"),
]

# Which protein a product name represents, used for per-segment protein selection.
PRODUCT_ROLE_RULES: List[Tuple[re.Pattern, str]] = [
    (re.compile(_NSS_STRICT), "NSs"),
    (re.compile(_RDRP), "RdRp"),
    (re.compile(_GPC), "GPC"),
    (re.compile(_N), "N"),
]

SEGMENT_PRIMARY_PRODUCT = {"S": "N", "M": "GPC", "L": "RdRp"}

# --------------------------------------------------------------------------------------
# Pathogenicity reference table
# --------------------------------------------------------------------------------------
# Every entry carries its source. Species absent from this table are labelled "unknown"
# and are NEVER silently treated as non-pathogenic.
#
# CFR figures genuinely conflict between sources (HTNV is reported both as ~1% and as
# 5-15%), so the reported range is stored verbatim alongside an ordinal band used for
# modelling. Do not collapse these to a point estimate.

CIT_ZOONOTIC = "https://pmc.ncbi.nlm.nih.gov/articles/PMC10459939/"
CIT_ECDC = "https://www.ecdc.europa.eu/en/infectious-disease-topics/hantavirus-infection/factsheet-orthohantavirus-infections"
CIT_REASSORT = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6153690/"
CIT_ICTV = "https://ictv.global/report/chapter/hantaviridae/hantaviridae/mammantavirinae/orthohantavirus"

# cfr_band ordinal: 0 = no human disease known, 1 = <1%, 2 = 1-15%, 3 = >30%
CFR_BAND_LABEL = {0: "none_known", 1: "low_lt_1pct", 2: "moderate_1_15pct", 3: "high_gt_30pct", -1: "unknown"}

SPECIES_RISK: Dict[str, Dict[str, Any]] = {
    "Orthohantavirus hantanense": {
        "abbrev": "HTNV", "syndrome": "HFRS", "cfr_band": 2,
        "cfr_reported": "~1% (PMC10459939); 5-15% (ECDC/other reviews)",
        "reservoir": "Apodemus agrarius", "region": "East Asia", "citation": CIT_ZOONOTIC},
    "Orthohantavirus dobravaense": {
        "abbrev": "DOBV", "syndrome": "HFRS", "cfr_band": 2,
        "cfr_reported": "0-15% depending on genotype (Dobrava severe, Kurkino/Saaremaa mild)",
        "reservoir": "Apodemus spp.", "region": "Europe", "citation": CIT_ZOONOTIC},
    "Orthohantavirus seoulense": {
        "abbrev": "SEOV", "syndrome": "HFRS", "cfr_band": 1,
        "cfr_reported": "<1%", "reservoir": "Rattus norvegicus", "region": "Worldwide",
        "citation": CIT_ZOONOTIC},
    "Orthohantavirus puumalaense": {
        "abbrev": "PUUV", "syndrome": "HFRS/NE", "cfr_band": 1,
        "cfr_reported": "0.1-0.4%", "reservoir": "Myodes glareolus", "region": "Europe",
        "citation": CIT_ZOONOTIC},
    "Orthohantavirus andesense": {
        "abbrev": "ANDV", "syndrome": "HCPS", "cfr_band": 3,
        "cfr_reported": "~35-40%; only orthohantavirus with documented person-to-person spread",
        "reservoir": "Oligoryzomys longicaudatus", "region": "Argentina, Chile",
        "citation": CIT_ZOONOTIC},
    "Orthohantavirus sinnombreense": {
        "abbrev": "SNV", "syndrome": "HCPS", "cfr_band": 3,
        "cfr_reported": "30-50%", "reservoir": "Peromyscus maniculatus",
        "region": "USA, Canada, Mexico", "citation": CIT_ZOONOTIC},
    "Orthohantavirus chocloense": {
        "abbrev": "CHOV", "syndrome": "HCPS", "cfr_band": 2,
        "cfr_reported": "12-15%", "reservoir": "Oligoryzomys fulvescens", "region": "Panama",
        "citation": CIT_ZOONOTIC},
    "Orthohantavirus negraense": {
        "abbrev": "LANV", "syndrome": "HCPS", "cfr_band": 2,
        "cfr_reported": "12-15%", "reservoir": "Calomys spp.", "region": "South America",
        "citation": CIT_ZOONOTIC},
    "Orthohantavirus mamorense": {
        "abbrev": "RIOMV", "syndrome": "HCPS", "cfr_band": -1,
        "cfr_reported": "human cases reported; CFR not well established",
        "reservoir": "Oligoryzomys microtis", "region": "South America", "citation": CIT_ZOONOTIC},
    "Orthohantavirus bayoui": {
        "abbrev": "BAYV", "syndrome": "HCPS", "cfr_band": -1,
        "cfr_reported": "human cases reported; CFR not well established",
        "reservoir": "Oryzomys palustris", "region": "USA", "citation": CIT_ZOONOTIC},
    "Orthohantavirus montanoense": {
        "abbrev": "MTNV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Peromyscus spp.",
        "region": "North America", "citation": CIT_ZOONOTIC},
    "Orthohantavirus tulaense": {
        "abbrev": "TULV", "syndrome": "rare_opportunistic", "cfr_band": 1,
        "cfr_reported": "a handful of mild human cases, mostly immunocompromised; encodes NSs",
        "reservoir": "Microtus spp.", "region": "Europe", "citation": CIT_ZOONOTIC},
    "Orthohantavirus prospectense": {
        "abbrev": "PHV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Microtus pennsylvanicus",
        "region": "North America", "citation": CIT_ZOONOTIC},
    "Orthohantavirus khabarovskense": {
        "abbrev": "KHAV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Microtus fortis",
        "region": "Russia", "citation": CIT_ZOONOTIC},
    "Orthohantavirus thailandense": {
        "abbrev": "THAIV", "syndrome": "suspected_HFRS", "cfr_band": -1,
        "cfr_reported": "serological evidence of human infection; CFR unknown",
        "reservoir": "Bandicota indica", "region": "Southeast Asia", "citation": CIT_ZOONOTIC},
    "Orthohantavirus sangassouense": {
        "abbrev": "SANGV", "syndrome": "suspected", "cfr_band": -1,
        "cfr_reported": "serological evidence only", "reservoir": "Hylomyscus simus",
        "region": "West Africa", "citation": CIT_ZOONOTIC},
    "Orthohantavirus asikkalaense": {
        "abbrev": "ASIV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Sorex minutus",
        "region": "Europe", "citation": CIT_ZOONOTIC},
    "Orthohantavirus dabieshanense": {
        "abbrev": "DBSV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Niviventer confucianus",
        "region": "China", "citation": CIT_ZOONOTIC},
    "Orthohantavirus jejuense": {
        "abbrev": "JJUV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "shrew-borne; regarded as non-pathogenic", "reservoir": "Crocidura shantungensis",
        "region": "South Korea", "citation": CIT_ZOONOTIC},
    "Orthohantavirus caobangense": {
        "abbrev": "CBNV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "shrew-borne; regarded as non-pathogenic", "reservoir": "Anourosorex squamipes",
        "region": "Vietnam, China", "citation": CIT_ZOONOTIC},
    "Orthohantavirus brugesense": {
        "abbrev": "BRGV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "mole-borne; regarded as non-pathogenic", "reservoir": "Talpa europaea",
        "region": "Europe", "citation": CIT_ZOONOTIC},
    "Orthohantavirus asamaense": {
        "abbrev": "ASAV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "mole-borne; regarded as non-pathogenic", "reservoir": "Urotrichus talpoides",
        "region": "Japan", "citation": CIT_ZOONOTIC},
    "Orthohantavirus rockportense": {
        "abbrev": "RKPV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "mole-borne; regarded as non-pathogenic", "reservoir": "Scalopus aquaticus",
        "region": "USA", "citation": CIT_ZOONOTIC},
    "Orthohantavirus artybashense": {
        "abbrev": "ARTV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "shrew-borne; regarded as non-pathogenic", "reservoir": "Sorex spp.",
        "region": "Russia", "citation": CIT_ZOONOTIC},
    "Orthohantavirus kenkemeense": {
        "abbrev": "KKMV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "shrew-borne; regarded as non-pathogenic", "reservoir": "Sorex roboratus",
        "region": "Russia", "citation": CIT_ZOONOTIC},
    "Orthohantavirus luxiense": {
        "abbrev": "LUXV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Eothenomys miletus",
        "region": "China", "citation": CIT_ZOONOTIC},
    "Orthohantavirus fugongense": {
        "abbrev": "FUGV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Eothenomys cachinus",
        "region": "China", "citation": CIT_ZOONOTIC},
    "Orthohantavirus lankaense": {
        "abbrev": "LNKV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Mus booduga",
        "region": "Sri Lanka", "citation": CIT_ZOONOTIC},
    "Orthohantavirus tatenalense": {
        "abbrev": "TATV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Microtus agrestis",
        "region": "Europe", "citation": CIT_ZOONOTIC},
    "Orthohantavirus delgaditoense": {
        "abbrev": "DELV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Sigmodon alstoni",
        "region": "Venezuela", "citation": CIT_ZOONOTIC},
    "Orthohantavirus maporalense": {
        "abbrev": "MAPV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Oligoryzomys delicatus",
        "region": "Venezuela", "citation": CIT_ZOONOTIC},
    "Orthohantavirus carrizalense": {
        "abbrev": "CARV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Reithrodontomys spp.",
        "region": "Mexico", "citation": CIT_ZOONOTIC},
    "Orthohantavirus ozarkense": {
        "abbrev": "OZKV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Peromyscus leucopus",
        "region": "USA", "citation": CIT_ZOONOTIC},
    "Orthohantavirus sagercreekense": {
        "abbrev": "SGCV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Microtus ochrogaster",
        "region": "USA", "citation": CIT_ZOONOTIC},
    "Orthohantavirus nigrorivense": {
        "abbrev": "BCCV", "syndrome": "HCPS", "cfr_band": -1,
        "cfr_reported": "Black Creek Canal virus; human cases reported, CFR not established",
        "reservoir": "Sigmodon hispidus", "region": "USA", "citation": CIT_ZOONOTIC},
    "Orthohantavirus boweense": {
        "abbrev": "BOWV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "shrew-borne; regarded as non-pathogenic", "reservoir": "Crocidura spp.",
        "region": "Africa", "citation": CIT_ZOONOTIC},
    "Orthohantavirus tigrayense": {
        "abbrev": "TIGV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "no confirmed human disease", "reservoir": "Stenocephalemys albipes",
        "region": "Ethiopia", "citation": CIT_ZOONOTIC},
    "Orthohantavirus wufangense": {
        "abbrev": "WFSV", "syndrome": "none_known", "cfr_band": 0,
        "cfr_reported": "shrew-borne; regarded as non-pathogenic", "reservoir": "Soriculus spp.",
        "region": "China", "citation": CIT_ZOONOTIC},
}

# Taxonomic order of each reservoir genus named in SPECIES_RISK. Used to detect records
# where a virus was found outside its expected reservoir order -- e.g. SEOV in Chinese
# shrews and bats, or ANDV-clade sequences in Brazilian bats. Both are real, documented
# detections in this dataset, and a pathogenic rodent-borne virus recovered from a
# non-rodent host is a host-jump signal worth surfacing rather than an error.
RESERVOIR_GENUS_ORDER = {
    # Rodentia
    "Apodemus": "Rodentia", "Rattus": "Rodentia", "Myodes": "Rodentia",
    "Oligoryzomys": "Rodentia", "Peromyscus": "Rodentia", "Calomys": "Rodentia",
    "Oryzomys": "Rodentia", "Microtus": "Rodentia", "Bandicota": "Rodentia",
    "Hylomyscus": "Rodentia", "Niviventer": "Rodentia", "Sigmodon": "Rodentia",
    "Reithrodontomys": "Rodentia", "Eothenomys": "Rodentia", "Mus": "Rodentia",
    "Stenocephalemys": "Rodentia", "Alexandromys": "Rodentia",
    # Eulipotyphla (shrews and moles)
    "Crocidura": "Eulipotyphla", "Anourosorex": "Eulipotyphla", "Talpa": "Eulipotyphla",
    "Urotrichus": "Eulipotyphla", "Scalopus": "Eulipotyphla", "Sorex": "Eulipotyphla",
    "Episoriculus": "Eulipotyphla", "Suncus": "Eulipotyphla", "Soriculus": "Eulipotyphla",
    "Chodsigoa": "Eulipotyphla",
}


def audit_reservoir_curation(accessions: Sequence["AccessionRecord"]) -> List[Dict[str, str]]:
    """Check every curated reservoir against the modal observed host order in the data.

    SPECIES_RISK reservoirs were curated by hand, and six of the original 38 entries
    disagreed with the package's own host annotations at the taxonomic-order level
    (asikkalaense, fugongense, lankaense, luxiense, sagercreekense, tigrayense). Those
    were corrected against the data, which comes from the submitting studies and is the
    better source. This audit runs on every invocation so the same class of curation
    error cannot silently return when SPECIES_RISK is next edited.
    """
    observed: DefaultDict[str, Counter] = defaultdict(Counter)
    hosts: DefaultDict[str, Counter] = defaultdict(Counter)
    for record in accessions:
        if record.host_order and record.host_order != "Primates":
            observed[record.species_name][record.host_order] += 1
            hosts[record.species_name][record.host_name] += 1

    findings: List[Dict[str, str]] = []
    for species, entry in SPECIES_RISK.items():
        counts = observed.get(species)
        if not counts:
            continue
        curated = RESERVOIR_GENUS_ORDER.get(entry["reservoir"].split()[0], "")
        modal = mode(counts)
        if not curated:
            findings.append({"species": species, "issue": "reservoir genus not in RESERVOIR_GENUS_ORDER",
                             "curated_reservoir": entry["reservoir"], "modal_observed_host": ""})
        elif curated != modal:
            findings.append({
                "species": species, "issue": "curated reservoir order disagrees with data",
                "curated_reservoir": f"{entry['reservoir']} ({curated})",
                "modal_observed_host": f"{mode(hosts[species])} ({modal})"})
    return findings


# Reservoir-host phylogroup, one level finer than RESERVOIR_GENUS_ORDER's taxonomic
# order. Used to group candidate isolates for tree-building at scope="phylogroup" (see
# group_candidates): a species-only tree cannot show reassortment acquired from a
# DIFFERENT species, because the donor segment has no possible match in a single-species
# tree (Laenen et al. 2019, 10.3390/v11090788, on interspecies exchange in hantavirid
# evolution). Grouping by phylogroup instead of raw host order keeps groups small enough
# to be biologically meaningful -- "all Rodentia" would still merge Old World murine and
# New World sigmodontine reservoirs that share no recent common ancestor among the
# reservoir hosts, even though their viruses are all nominally "rodent-borne".
#
# Subfamily split is standard rodent taxonomy (Wilson & Reeder, Mammal Species of the
# World, 3rd ed., 2005), not something specific to this dataset or citable to a single
# hantavirus paper -- textbook classification, applied here rather than re-derived.
GENUS_PHYLOGROUP = {
    # Murinae (Old World rats and mice)
    "Apodemus": "murinae_borne", "Rattus": "murinae_borne", "Bandicota": "murinae_borne",
    "Hylomyscus": "murinae_borne", "Niviventer": "murinae_borne", "Mus": "murinae_borne",
    "Stenocephalemys": "murinae_borne",
    # Arvicolinae (voles, lemmings and relatives)
    "Myodes": "arvicolinae_borne", "Microtus": "arvicolinae_borne",
    "Alexandromys": "arvicolinae_borne", "Eothenomys": "arvicolinae_borne",
    # Sigmodontinae (New World rats and mice -- the HCPS reservoirs)
    "Oligoryzomys": "sigmodontinae_borne", "Peromyscus": "sigmodontinae_borne",
    "Calomys": "sigmodontinae_borne", "Oryzomys": "sigmodontinae_borne",
    "Sigmodon": "sigmodontinae_borne", "Reithrodontomys": "sigmodontinae_borne",
    # Eulipotyphla (shrews and moles). One phylogroup for the whole order, unlike
    # Rodentia's three-way subfamily split: every Eulipotyphla genus in this dataset is
    # shrew- or mole-borne with no pathogenic representative, so a finer split would add
    # groups without adding a distinction the risk/reassortment analysis uses.
    "Crocidura": "soricomorpha_borne", "Anourosorex": "soricomorpha_borne",
    "Talpa": "soricomorpha_borne", "Urotrichus": "soricomorpha_borne",
    "Scalopus": "soricomorpha_borne", "Sorex": "soricomorpha_borne",
    "Episoriculus": "soricomorpha_borne", "Suncus": "soricomorpha_borne",
    "Soriculus": "soricomorpha_borne", "Chodsigoa": "soricomorpha_borne",
}

# species_name -> phylogroup, derived from SPECIES_RISK's curated reservoir genus rather
# than hand-duplicated, so the two tables cannot drift the way SPECIES_RISK and
# RESERVOIR_GENUS_ORDER already once did (see audit_reservoir_curation). A species whose
# reservoir genus isn't in GENUS_PHYLOGROUP falls to 'unassigned' rather than guessing --
# none currently do (every SPECIES_RISK reservoir genus is covered above), but a future
# SPECIES_RISK entry with an uncovered genus should surface as 'unassigned', not silently
# join the wrong group. No species maps to 'chiroptera_borne' via this route: no curated
# reservoir in SPECIES_RISK is bat-borne. It stays a valid phylogroup for Task 3's
# host-order-based fallback, where unclassified isolates observed in bat hosts need
# somewhere to go.
PHYLOGROUP: Dict[str, str] = {
    species: GENUS_PHYLOGROUP.get(entry["reservoir"].split()[0], "unassigned")
    for species, entry in SPECIES_RISK.items()
}

# Single source of truth for the mammalian orders this dataset contains: the host role
# each implies, and the baseline human-risk prior. Shrew (Eulipotyphla), mole and bat
# borne hantaviruses are regarded as non-pathogenic to humans (PMC10459939).
# 'Soricomorpha' and 'Insectivora' are superseded names for Eulipotyphla that still
# appear in older NCBI lineages.
HOST_ORDER_ROLE = {
    "Rodentia": "reservoir_rodent",
    "Eulipotyphla": "reservoir_shrew_mole",
    "Soricomorpha": "reservoir_shrew_mole",
    "Insectivora": "reservoir_shrew_mole",
    "Chiroptera": "reservoir_bat",
    "Primates": "spillover_host_primate",
    "Carnivora": "",
    "Lagomorpha": "",
    "Artiodactyla": "",
}

HOST_ORDER_RISK_PRIOR = {
    "Rodentia": "rodent_borne_risk_possible",
    "Eulipotyphla": "shrew_mole_borne_low_risk",
    # Same superseded names as in HOST_ORDER_ROLE. They carry no records in this package
    # but must resolve identically to Eulipotyphla if an older lineage ever appears,
    # otherwise the same shrew is 'low risk' in one table and 'unknown' in another.
    "Soricomorpha": "shrew_mole_borne_low_risk",
    "Insectivora": "shrew_mole_borne_low_risk",
    "Chiroptera": "bat_borne_low_risk",
    "Primates": "human_or_primate_derived",
}

# Disease syndrome -> pathogenicity class. TULV's 'rare_opportunistic' is kept out of
# pathogenic_hfrs deliberately: lumping its handful of mild case reports with HTNV/DOBV
# would put 2,123 records into the same class as the viruses that drive HFRS mortality.
SYNDROME_TO_PATHOGENICITY = {
    "HCPS": "pathogenic_hcps",
    "HFRS": "pathogenic_hfrs",
    "HFRS/NE": "pathogenic_hfrs",
    "rare_opportunistic": "rare_or_opportunistic",
    "suspected": "suspected",
    "suspected_HFRS": "suspected",
    "none_known": "no_known_human_disease",
}

# Weights for surveillance_priority. Named here because the score is explicitly meant to
# be audited and re-weighted downstream.
SURVEILLANCE_WEIGHTS = {
    "cfr_by_band": {3: 4.0, 2: 2.5, 1: 1.0, 0: 0.0, -1: 0.5},
    "human_spillover": 2.0,
    "taxonomic_incongruence": 1.5,
    "host_jump": 1.0,
    "complete_triplet": 1.0,
    "nss_orf": 0.5,
    # Provisional -- not fitted to anything, chosen to sit below taxonomic_incongruence
    # (1.5, a metadata FACT) and above complete_triplet (1.0), since bootstrap-supported
    # tree-topology evidence is a real inference but a noisier one than a fact NCBI itself
    # recorded. Only contributes for reassortment_call == 'supported_reassortment_candidate'
    # -- see surveillance_priority. Re-tune once enough --phylogenetics runs exist to know
    # how often this term actually fires.
    "supported_reassortment": 1.5,
}

# Reassortment frequency reported in the literature, used only as a prior for reporting
# and never as a substitute for the per-isolate signal. Source: PMC6153690.
REPORTED_REASSORTMENT = {
    "Orthohantavirus puumalaense": "frequent (19.1-32% in Finnish studies)",
    "Orthohantavirus sinnombreense": "frequent (M and S segment)",
    "Orthohantavirus dobravaense": "frequent",
    "Orthohantavirus hantanense": "documented (SAMALB pattern in Korean genomes)",
    "Orthohantavirus seoulense": "documented (between lineages, central-south China)",
    "Orthohantavirus tulaense": "frequent",
    "Orthohantavirus andesense": "absent or sparse",
}


# --------------------------------------------------------------------------------------
# Dataclasses
# --------------------------------------------------------------------------------------


@dataclass
class InputPaths:
    root: Path
    package_dir: Path
    data_dir: Path
    data_report: Path
    annotation_report: Path
    biosample_report: Path
    dataset_catalog: Path
    md5sum_file: Optional[Path]
    genomic_fasta: Optional[Path]
    protein_fasta: Optional[Path]
    cds_fasta: Optional[Path]
    output_dir: Path
    log_file: Path


@dataclass
class FastaRecord:
    accession: str
    accession_base: str
    header: str
    sequence: str
    length: int


@dataclass
class ChecksumResult:
    file_name: str
    expected_md5: str
    observed_md5: str
    verified: bool
    matched_form: str          # "raw" | "crlf_normalized" | "none" | "unlisted" | "missing"
    diagnosis: str


@dataclass
class SpeciesAssignment:
    canonical_name: str
    canonical_taxid: str
    observed_taxid: str
    observed_name: str
    status: str                # "classified" | "unclassified" | "genus_only" | "outside_genus"
    ranks: Dict[str, str]


@dataclass
class ProteinHit:
    role: str
    product: str
    protein_accession: str
    sequence: str


@dataclass
class BioSampleInfo:
    accession: str = ""
    host: Optional[str] = None
    host_taxid: Optional[str] = None
    isolation_source: Optional[str] = None
    tissue: Optional[str] = None
    collection_date: Optional[str] = None
    geographic_location: Optional[str] = None
    strain: Optional[str] = None
    isolate: Optional[str] = None
    raw_attributes: Dict[str, str] = field(default_factory=dict)


@dataclass
class AccessionRecord:
    accession: str
    sequence_hash: str
    species_name: str
    species_taxid: str
    species_status: str
    observed_organism: str
    observed_taxid: str
    ranks: Dict[str, str]
    segment: str
    segment_source: str
    segment_confidence: float
    completeness: str
    completeness_source: str
    is_coding_complete: bool
    nt_length: int
    biosample: str
    bioproject: str
    sra_accessions: List[str]
    isolate: str
    strain: str
    collection_date: str
    year: str
    location: str
    country: str
    region: str
    host_name: str
    host_taxid: str
    host_order: str
    host_family: str
    host_source: str
    host_role: str
    isolation_source: str
    tissue_source: str
    is_lab_host: bool
    lab_host: str
    is_vaccine_strain: bool
    protein_accessions: List[str]
    annotation_products: List[str]
    protein_aa: str
    protein_product: str
    has_nss_orf: bool
    group_key: str
    grouping_tier: str
    group_confidence: float
    is_refseq: bool
    source_database: str
    warning_flags: List[str]


@dataclass
class IsolateRecord:
    isolate_key: str
    species_name: str
    species_taxid: str
    species_status: str
    biosample: str
    isolate: str
    strain: str
    host_name: str
    host_taxid: str
    host_order: str
    host_family: str
    host_role: str
    collection_date: str
    year: str
    location: str
    country: str
    region: str
    grouping_tier: str
    grouping_confidence: float
    s_accession: str
    m_accession: str
    l_accession: str
    s_completeness: str
    m_completeness: str
    l_completeness: str
    s_nt_length: int
    m_nt_length: int
    l_nt_length: int
    s_protein_aa: str
    m_protein_aa: str
    l_protein_aa: str
    concatenated_protein_aa: str
    segments_present: List[str]
    coding_complete_segments: List[str]
    triplet_status: str
    has_complete_triplet: bool
    is_lab_derived: bool
    is_vaccine_strain: bool
    has_nss_orf: bool
    metadata_conflict: bool
    taxonomic_incongruence: bool
    cross_species_partners: List[str]
    kmer_cluster_id: str
    s_protein_accessions: List[str]
    m_protein_accessions: List[str]
    l_protein_accessions: List[str]
    all_accessions: List[str]
    warnings: List[str]


@dataclass
class SegmentDivergenceProfile:
    """Per-segment nearest-neighbour diagnostics for one isolate.

    Deliberately NOT a reassortment call. See compute_segment_divergence for why the
    reassortment classifier that used to live here was withdrawn.
    """
    isolate_key: str
    species_name: str
    assessable: bool
    reason: str
    s_neighbour: str
    m_neighbour: str
    l_neighbour: str
    s_distance: float
    m_distance: float
    l_distance: float
    sl_rank_separation: float
    m_rank_offset: float
    taxonomic_incongruence: bool
    reassortment_status: str
    reported_species_prior: str




# --------------------------------------------------------------------------------------
# Generic utilities
# --------------------------------------------------------------------------------------


def configure_logging(log_file: Path, verbose: bool = False) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        # mode="w": one run, one log. Appending made this a 428-line accretion of every run
        # ever done in the checkout, and since it is a tracked file it carried absolute paths
        # (and so usernames) from more than one developer's machine into git.
        handlers=[logging.FileHandler(log_file, mode="w", encoding="utf-8"),
                  logging.StreamHandler(sys.stdout)],
    )


def clean_text(value: Any) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, (dict, list)):
        text = json.dumps(value, ensure_ascii=False)
    else:
        text = str(value).strip()
    if not text:
        return None
    if text.lower() in {"none", "null", "nan", "na", "n/a", "unknown", "missing",
                        "not collected", "not provided", "not applicable", "restricted access"}:
        return None
    return WHITESPACE_RE.sub(" ", text).strip()


def first_nonempty(*values: Any) -> Optional[str]:
    for value in values:
        if isinstance(value, list):
            for item in value:
                found = clean_text(item)
                if found:
                    return found
        else:
            found = clean_text(value)
            if found:
                return found
    return None


def uniq(values: Iterable[Any]) -> List[Any]:
    seen = set()
    out = []
    for value in values:
        key = json.dumps(value, ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list, tuple, set)) else str(value)
        if key not in seen:
            seen.add(key)
            out.append(value)
    return out


def mode(counter: Counter, default: str = "") -> str:
    """Most frequent key, or `default` when the counter is empty."""
    return counter.most_common(1)[0][0] if counter else default


def accession_base(accession: str) -> str:
    return accession.split(".")[0] if accession else ""


def indexed_record_count(index: Dict[str, Any]) -> int:
    """Number of distinct records in an index that files each one under several keys.

    build_annotation_index and parse_fasta both key a record on its versioned accession
    and on its base, so len(index) reports roughly twice the number of records read.
    """
    return len({id(value) for value in index.values()})


def iter_jsonl(path: Path) -> Iterator[Dict[str, Any]]:
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"JSON decode error in {path.name} line {line_no}: {exc}") from exc
            if isinstance(obj, dict):
                yield obj


def extract_year(text: Optional[str]) -> str:
    if not text:
        return ""
    match = YEAR_RE.search(text)
    return match.group(0) if match else ""


def split_geo(geo: str) -> Tuple[str, str]:
    if not geo:
        return "", ""
    if ":" in geo:
        country, region = geo.split(":", 1)
        return country.strip(), region.strip()
    return geo.strip(), ""


def normalize_key_text(value: str) -> str:
    """Collapse case and punctuation so 'PUUV/Foo-01' and 'puuv foo 01' agree.

    Separator runs become a single space rather than being deleted. Deleting them merges
    52 genuine spelling variants in this package but also collapses 68 pairs that differ
    only at a digit boundary -- '11-4' vs '114', 'B-1' vs 'B1', 'Aa-08-1111' vs
    'Aa08/1111' -- which are not safe to treat as the same sample.
    """
    return WHITESPACE_RE.sub(" ", NON_ALNUM_RE.sub(" ", (value or "").lower())).strip()


CROSS_STUDY_LABEL_MIN_CHARS = 4


def is_cross_study_label(label: str) -> bool:
    """Is an isolate label distinctive enough to identify one sample across submissions?

    Submitters number their own samples independently, so labels like '2', '43' or 'B1'
    collide between unrelated studies. Treating such a collision as evidence that one
    isolate was assigned to two species manufactures incongruence: 31 of the 70 isolates
    previously flagged carried a purely numeric label, including five distinct isolates
    named simply '2' spread across PUUV, SEOV, SNV, CBNV and unclassified records.

    Grouping is not subject to this: a group key also contains the species and any
    over-merge it causes surfaces as metadata_conflict (4 of 254 numeric-label isolates).
    A cross-species check has no such backstop, and its output feeds the risk score.
    """
    key = normalize_key_text(label).replace(" ", "")
    return len(key) >= CROSS_STUDY_LABEL_MIN_CHARS and not key.isdigit()


def _atomic_write_bytes(path: Path, data: bytes) -> None:
    """Write `data` to `path` without ever leaving a truncated file at that name.

    Writes to a sibling temp file first (same directory, so the final `os.replace` is a
    single filesystem-level rename -- atomic on POSIX), fsyncing before the rename so the
    bytes are actually on disk rather than just buffered. A crash, kill, or power loss
    mid-write can only ever destroy the temp file; `path` itself is either the previous
    good version or absent, never a partial one.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, path)
    except BaseException:
        Path(tmp_name).unlink(missing_ok=True)
        raise


def safe_write_csv(path: Path, rows: Sequence[Dict[str, Any]]) -> str:
    """Write a CSV table atomically; returns its sha256 for the output checksum manifest."""
    frame = pd.DataFrame(list(rows)) if rows else pd.DataFrame()
    data = frame.to_csv(index=False).encode("utf-8")
    _atomic_write_bytes(path, data)
    logging.info("Wrote %s rows x %s cols to %s", len(frame), len(frame.columns), path.name)
    return hashlib.sha256(data).hexdigest()


def safe_write_json(path: Path, payload: Any) -> str:
    """Write JSON atomically; returns its sha256 for the output checksum manifest."""
    data = json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8")
    _atomic_write_bytes(path, data)
    return hashlib.sha256(data).hexdigest()


def write_output_checksums(path: Path, digests: Dict[str, str]) -> None:
    """Write a sha256sum(1)-compatible manifest so output integrity can be checked later.

    Verify with: sha256sum -c ortho_output_checksums.txt (run from the output directory).
    """
    lines = [f"{digest}  {name}\n" for name, digest in sorted(digests.items())]
    _atomic_write_bytes(path, "".join(lines).encode("utf-8"))


def _iter_fasta_records(path: Optional[Path]) -> Iterator[Tuple[str, str]]:
    """Yield (header-without-'>', concatenated sequence) pairs from a FASTA file."""
    if path is None or not path.exists():
        return
    header: Optional[str] = None
    parts: List[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
        if header is not None:
            yield header, "".join(parts)


def parse_fasta(path: Optional[Path]) -> Dict[str, FastaRecord]:
    index: Dict[str, FastaRecord] = {}
    for desc, seq in _iter_fasta_records(path):
        if not desc:
            continue
        acc = desc.split()[0]
        # Protein/CDS deflines carry a coordinate suffix, e.g. 'UCR75372.1:1-338'.
        acc = acc.split(":")[0]
        base = accession_base(acc)
        record = FastaRecord(accession=acc, accession_base=base, header=desc, sequence=seq, length=len(seq))
        index.setdefault(acc, record)
        if base:
            index.setdefault(base, record)
    return index


# --------------------------------------------------------------------------------------
# Integrity gate -- verify against NCBI's published checksums
# --------------------------------------------------------------------------------------


def file_md5(path: Path, normalize_crlf: bool = False) -> str:
    """MD5 of a file. With normalize_crlf, CRLF is collapsed to LF before hashing.

    The normalized form exists only to *diagnose* a checkout whose line endings were
    rewritten (e.g. by git core.autocrlf on Windows). It is never treated as success.
    """
    digest = hashlib.md5()
    with path.open("rb") as handle:
        if not normalize_crlf:
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                digest.update(chunk)
        else:
            tail = b""
            for chunk in iter(lambda: handle.read(1 << 20), b""):
                data = tail + chunk
                # Hold back a trailing CR in case the LF starts the next chunk.
                if data.endswith(b"\r"):
                    tail, data = b"\r", data[:-1]
                else:
                    tail = b""
                digest.update(data.replace(b"\r\n", b"\n"))
            if tail:
                digest.update(tail)
    return digest.hexdigest()


def parse_md5sum_file(path: Path) -> Dict[str, str]:
    """Parse NCBI's md5sum.txt ('<md5>  <relative/path>' per line) into {basename: md5}."""
    expected: Dict[str, str] = {}
    if not path.exists():
        return expected
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        digest, rel_path = parts[0].strip(), parts[1].strip().lstrip("*")
        if MD5_RE.fullmatch(digest):
            expected[Path(rel_path).name] = digest.lower()
    return expected


def verify_inputs(paths: InputPaths) -> List[ChecksumResult]:
    """Verify every shipped file against md5sum.txt, diagnosing line-ending corruption."""
    if paths.md5sum_file is None or not paths.md5sum_file.exists():
        logging.warning("No md5sum.txt found in package; integrity cannot be verified against NCBI.")
        return []

    expected = parse_md5sum_file(paths.md5sum_file)
    if not expected:
        logging.warning("md5sum.txt present but no usable entries were parsed.")
        return []

    results: List[ChecksumResult] = []
    for candidate in sorted(paths.data_dir.iterdir()):
        if not candidate.is_file():
            continue
        want = expected.get(candidate.name)
        if want is None:
            results.append(ChecksumResult(candidate.name, "", "", True, "unlisted",
                                          "not listed in md5sum.txt; skipped"))
            continue
        observed = file_md5(candidate)
        if observed == want:
            results.append(ChecksumResult(candidate.name, want, observed, True, "raw", "ok"))
            continue
        normalized = file_md5(candidate, normalize_crlf=True)
        if normalized == want:
            results.append(ChecksumResult(
                candidate.name, want, observed, False, "crlf_normalized",
                "on-disk bytes are CRLF but NCBI published the LF form; the working copy has been "
                "rewritten (typically git core.autocrlf=true without a .gitattributes rule). "
                "Fix with the .gitattributes shipped alongside this script, then re-checkout the package."))
        else:
            results.append(ChecksumResult(
                candidate.name, want, observed, False, "none",
                "content does not match NCBI's published checksum in either raw or LF-normalized "
                "form; the file is truncated, modified or from a different download."))

    # Walking the directory only ever sees files that exist, so a file listed in md5sum.txt
    # but deleted from the checkout would produce no result at all and pass the gate in
    # silence. resolve_inputs already hard-fails on the four required reports, but the
    # FASTAs are optional there and parse_fasta treats a missing file as an empty index --
    # the ML tables would quietly shrink instead of aborting.
    seen = {result.file_name for result in results}
    for name, want in sorted(expected.items()):
        if name not in seen:
            results.append(ChecksumResult(
                name, want, "", False, "missing",
                "listed in NCBI's md5sum.txt but not present on disk; the package is "
                "incomplete. Re-download or re-checkout it."))
    return results


def report_checksums(results: Sequence[ChecksumResult], allow_mismatch: bool) -> None:
    failures = [r for r in results if not r.verified]
    for result in results:
        if result.verified and result.matched_form == "raw":
            logging.debug("checksum ok: %s", result.file_name)
    listed = [r for r in results if r.matched_form != "unlisted"]
    unlisted = len(results) - len(listed)
    if not failures:
        # Report against the number of files md5sum.txt actually covers. Counting only the
        # files that passed on both sides of the ratio made this a tautological 'N/N' that
        # stayed silent about files shipped without a published checksum.
        logging.info("Integrity: %s/%s files verified against NCBI md5sum.txt%s",
                     sum(1 for r in listed if r.verified), len(listed),
                     f"; {unlisted} shipped file(s) not listed and not checked" if unlisted else "")
        return

    for result in failures:
        logging.error("CHECKSUM MISMATCH %s", result.file_name)
        logging.error("    expected (NCBI): %s", result.expected_md5)
        logging.error("    observed (disk): %s", result.observed_md5)
        logging.error("    diagnosis      : %s", result.diagnosis)

    if not allow_mismatch:
        raise SystemExit(
            f"\nAborting: {len(failures)} input file(s) do not match NCBI's published checksums.\n"
            "Results computed from these bytes would not be reproducible on another platform.\n"
            "Re-run with --allow-checksum-mismatch to proceed anyway (the mismatch is recorded "
            "in ortho_summary.json)."
        )
    logging.warning("Continuing despite %s checksum failure(s) because --allow-checksum-mismatch was set.",
                    len(failures))


def normalize_sequence_hash(value: Any) -> str:
    """NCBI emits nucleotide.sequenceHash as CRC32 with leading zeros stripped.

    Observed widths in this package are 5-8 characters, so zero-pad to 8 before any
    comparison or grouping -- otherwise '5F2C0FB' and '05F2C0FB' look like different
    sequences.
    """
    text = clean_text(value)
    if not text:
        return ""
    text = text.strip().upper()
    if not HEX_HASH_RE.fullmatch(text):
        return text
    return text.zfill(8)


# --------------------------------------------------------------------------------------
# Taxonomy -- resolve the ICTV species from virus.lineage
# --------------------------------------------------------------------------------------


def resolve_species(record: Dict[str, Any]) -> SpeciesAssignment:
    """Resolve a record to its canonical ICTV species.

    The organism-level ``virus.taxId`` is unusable as a species key: in this package it
    takes 553 distinct values, mixing strain-level taxa ('Hantaan virus Lee', 'Four
    Corners hantavirus') with legacy/current taxid pairs for the same species
    ('Puumala orthohantavirus' 1980486 vs 'Orthohantavirus puumalaense' 3052493).

    The node immediately below genus Orthohantavirus in ``virus.lineage`` *is* the ICTV
    species, and canonicalising on its binomial name merges the legacy/current pairs.
    That collapses 553 organism values to 38 species plus 'unclassified Orthohantavirus'.
    """
    virus = record.get("virus") or {}
    observed_name = clean_text(virus.get("organismName")) or ""
    observed_taxid = str(virus.get("taxId") or "")
    lineage = virus.get("lineage") or []

    ranks: Dict[str, str] = {}
    known_ranks = {name: rank for rank, name in LINEAGE_RANKS}
    genus_index = -1
    for index, node in enumerate(lineage):
        if not isinstance(node, dict):
            continue
        name = clean_text(node.get("name")) or ""
        if name in known_ranks:
            ranks[known_ranks[name]] = name
        if node.get("taxId") == ORTHOHANTAVIRUS_TAXID:
            genus_index = index

    if genus_index < 0:
        return SpeciesAssignment(observed_name or "unassigned", observed_taxid, observed_taxid,
                                 observed_name, "outside_genus", ranks)

    if genus_index + 1 >= len(lineage):
        return SpeciesAssignment("Orthohantavirus (genus-level only)", str(ORTHOHANTAVIRUS_TAXID),
                                 observed_taxid, observed_name, "genus_only", ranks)

    species_node = lineage[genus_index + 1]
    species_name = clean_text(species_node.get("name")) or ""
    species_taxid = str(species_node.get("taxId") or "")

    status = "classified"
    if species_name.lower().startswith("unclassified"):
        status = "unclassified"

    ranks["species"] = species_name
    return SpeciesAssignment(species_name, species_taxid, observed_taxid, observed_name, status, ranks)


# Cache for resolve_host_taxonomy, keyed on host taxid. Verified safe against this
# package: all 253 distinct host taxids carry a consistent organismName and lineage, so
# one walk per taxid replaces one walk per record (15,991 -> 253).
#
# Deliberately NOT done for resolve_species: 3 of 551 virus taxids carry more than one
# organismName or lineage (e.g. taxid 1980456 appears as both 'Orthohantavirus andesense'
# and 'Andes orthohantavirus'), and even (taxid, organismName) still collides for
# 'Lanka hantavirus'. Caching there would silently change the observed_organism column.
_HOST_TAXONOMY_CACHE: Dict[str, Tuple[str, str, str, str]] = {}


def resolve_host_taxonomy(record: Dict[str, Any]) -> Tuple[str, str, str, str]:
    """Return (host_name, host_taxid, order, family) from host.lineage."""
    host = record.get("host") or {}
    taxid = str(host.get("taxId") or "")
    cached = _HOST_TAXONOMY_CACHE.get(taxid) if taxid else None
    if cached is not None:
        return cached

    name = clean_text(host.get("organismName")) or ""
    order = ""
    family = ""
    order_names = HOST_ORDER_ROLE.keys()
    for node in host.get("lineage") or []:
        if not isinstance(node, dict):
            continue
        node_name = clean_text(node.get("name")) or ""
        if node_name in order_names:
            order = node_name
        elif node_name.endswith("idae"):
            # lineage runs broad -> narrow, so the last -idae node is the family rank.
            family = node_name

    result = (name, taxid, order, family)
    if taxid:
        _HOST_TAXONOMY_CACHE[taxid] = result
    return result


# --------------------------------------------------------------------------------------
# Segment resolution
# --------------------------------------------------------------------------------------


def normalize_segment_field(value: Any) -> str:
    """Normalize NCBI's free-text ``segment`` field to S / M / L.

    Values observed in this package beyond plain S/M/L include 'Partial L', 'S segment',
    'M-segment', 'M; medium', 'small (S)', 'medium (M)', 'G1', 'G2', 'ANONYMOUS', 's'
    and 'N'. The previous exact-match check ({'S','M','L','SMALL','MEDIUM','LARGE'})
    silently dropped all of these into the length heuristic.
    """
    text = clean_text(value)
    if not text:
        return ""
    low = SEGMENT_STRIP_RE.sub(" ", text.lower())
    low = SEGMENT_PUNCT_RE.sub(" ", low)
    low = WHITESPACE_RE.sub(" ", low).strip()
    if not low or low in {"anonymous", "unknown", "unnamed"}:
        return ""
    for pattern, segment in SEGMENT_FIELD_RULES:
        if pattern.search(low):
            return segment
    return ""


def segment_from_products(products: Sequence[str]) -> str:
    """Infer segment from annotated CDS product names, by majority vote."""
    votes: Counter = Counter()
    for product in products:
        low = (product or "").lower()
        for pattern, segment in PRODUCT_SEGMENT_RULES:
            if pattern.search(low):
                votes[segment] += 1
                break
    if not votes:
        return ""
    top = votes.most_common()
    if len(top) > 1 and top[0][1] == top[1][1]:
        return ""
    return top[0][0]


def segment_from_length(nt_len: int) -> str:
    """Last-resort segment call, only when the length falls in exactly one ICTV band."""
    if not nt_len:
        return ""
    matches = [seg for seg, (lo, hi) in SEGMENT_LENGTH_RANGE.items() if lo <= nt_len <= hi]
    return matches[0] if len(matches) == 1 else ""


def classify_protein_role(product: str) -> str:
    low = (product or "").lower()
    for pattern, role in PRODUCT_ROLE_RULES:
        if pattern.search(low):
            return role
    return ""


def choose_segment(
    raw_segment: Any,
    products: Sequence[str],
    defline: str,
    nt_len: int,
) -> Tuple[str, str, float, List[str]]:
    """Resolve segment with an explicit evidence hierarchy and a confidence score."""
    warnings: List[str] = []

    explicit = normalize_segment_field(raw_segment)
    if explicit:
        return explicit, "explicit", 1.0, warnings

    from_products = segment_from_products(products)
    if from_products:
        return from_products, "annotation", 0.90, warnings

    from_defline = segment_from_products([defline]) if defline else ""
    if from_defline:
        return from_defline, "defline", 0.75, warnings

    from_length = segment_from_length(nt_len)
    if from_length:
        warnings.append("segment_from_length_only")
        return from_length, "length_heuristic", 0.40, warnings

    warnings.append("segment_unresolved")
    return "", "unresolved", 0.0, warnings


# --------------------------------------------------------------------------------------
# Completeness
# --------------------------------------------------------------------------------------


def determine_completeness(record: Dict[str, Any], segment: str, nt_len: int) -> Tuple[str, str, bool]:
    """Return (completeness, source, is_coding_complete).

    NCBI's structured ``completeness`` is present on 100% of records in this package
    (PARTIAL 14,652 / COMPLETE 1,339) and is authoritative. The previous implementation
    carried an elaborate cascade of definition-line text matching, ambiguous-base counting
    and terminal-panhandle motif checks that never executed. Length is used only to
    qualify a COMPLETE call for coding-completeness, never to override NCBI.
    """
    raw = clean_text(record.get("completeness")) or ""
    low = raw.lower()
    if low == "complete":
        if not segment:
            return "complete", "ncbi", False
        coding_complete = nt_len >= COMPLETE_LENGTH_MIN.get(segment, 10 ** 9)
        return "complete", "ncbi", coding_complete
    if low == "partial":
        return "partial", "ncbi", False
    if low:
        return low, "ncbi", False
    return "unknown", "absent", False


def length_anomaly(segment: str, nt_len: int, completeness: str) -> str:
    """Flag lengths that contradict the ICTV segment size ranges.

    Only complete records are checked against the lower bound: a partial sequence is
    *expected* to be short, and the previous rule ('M' and not 1500 <= n <= 5000) fired
    on 5,856 records, almost all of them legitimately partial.
    """
    if not segment or not nt_len:
        return ""
    absurd_max = SEGMENT_LENGTH_ABSURD_MAX.get(segment)
    if absurd_max and nt_len > absurd_max:
        return "length_exceeds_segment_maximum"
    if completeness == "complete":
        lo, _hi = SEGMENT_LENGTH_RANGE[segment]
        if nt_len < lo * COMPLETE_LENGTH_TOLERANCE:
            return "complete_but_shorter_than_ictv_range"
    return ""


# --------------------------------------------------------------------------------------
# Annotation and BioSample indexes
# --------------------------------------------------------------------------------------


def build_annotation_index(records: Iterable[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Index annotation records by their own nucleotide accession.

    The previous implementation serialised each record to JSON, regex-scanned it for
    anything accession-shaped, and filed the record under *every* hit -- including the
    protein accessions -- which both mis-attributed annotations and made the index far
    larger than the input.
    """
    index: Dict[str, Dict[str, Any]] = {}
    for record in records:
        accession = clean_text(record.get("accession"))
        if not accession:
            continue
        index[accession] = record
        base = accession_base(accession)
        if base:
            index.setdefault(base, record)
    return index


def extract_annotation_details(
    annotation: Optional[Dict[str, Any]],
    protein_index: Dict[str, FastaRecord],
) -> Tuple[List[str], List[str], List[ProteinHit]]:
    """Return (product_names, protein_accessions, protein_hits) for one nucleotide record."""
    if not annotation:
        return [], [], []
    products: List[str] = []
    protein_accessions: List[str] = []
    hits: List[ProteinHit] = []
    for gene in annotation.get("genes") or []:
        if not isinstance(gene, dict):
            continue
        for cds in gene.get("cds") or []:
            if not isinstance(cds, dict):
                continue
            product = clean_text(cds.get("name")) or ""
            if product:
                products.append(product)
            protein = cds.get("protein") or {}
            protein_acc = clean_text(protein.get("accessionVersion")) or ""
            if protein_acc:
                protein_accessions.append(protein_acc)
            role = classify_protein_role(product)
            if not role:
                continue
            fasta = protein_index.get(protein_acc) or protein_index.get(accession_base(protein_acc))
            hits.append(ProteinHit(role=role, product=product, protein_accession=protein_acc,
                                   sequence=fasta.sequence if fasta else ""))
    return uniq(products), uniq(protein_accessions), hits


def select_primary_protein(hits: Sequence[ProteinHit], segment: str) -> Tuple[str, str, str]:
    """Pick the longest protein matching the segment's primary product."""
    wanted = SEGMENT_PRIMARY_PRODUCT.get(segment)
    if not wanted:
        return "", "", ""
    candidates = [hit for hit in hits if hit.role == wanted and hit.sequence]
    if not candidates:
        return "", "", ""
    best = max(candidates, key=lambda hit: len(hit.sequence))
    return best.sequence, best.product, best.protein_accession


def normalize_attr_name(name: str) -> str:
    return ATTR_NAME_RE.sub("_", name.strip().lower()).strip("_")


def parse_biosample_record(record: Dict[str, Any]) -> BioSampleInfo:
    accession = clean_text(record.get("accession")) or ""
    info = BioSampleInfo(accession=accession)
    attr_map: Dict[str, str] = {}
    for item in record.get("attributes") or []:
        if not isinstance(item, dict):
            continue
        name = normalize_attr_name(str(item.get("name", "")))
        value = clean_text(item.get("value"))
        if name and value:
            attr_map[name] = value
    info.raw_attributes = attr_map
    info.host = first_nonempty(attr_map.get("host"), attr_map.get("host_scientific_name"))
    info.host_taxid = first_nonempty(attr_map.get("host_taxid"))
    info.isolation_source = first_nonempty(attr_map.get("isolation_source"), attr_map.get("source"))
    info.tissue = first_nonempty(attr_map.get("tissue"), attr_map.get("host_tissue_sampled"))
    info.collection_date = first_nonempty(attr_map.get("collection_date"))
    info.geographic_location = first_nonempty(attr_map.get("geo_loc_name"), attr_map.get("geographic_location"))
    info.strain = first_nonempty(attr_map.get("strain"))
    info.isolate = first_nonempty(attr_map.get("isolate"), attr_map.get("submitter_id"))
    return info


def build_biosample_index(records: Iterable[Dict[str, Any]]) -> Dict[str, BioSampleInfo]:
    index: Dict[str, BioSampleInfo] = {}
    for record in records:
        info = parse_biosample_record(record)
        if info.accession:
            index[info.accession.upper()] = info
    return index


# --------------------------------------------------------------------------------------
# Host role
# --------------------------------------------------------------------------------------


def determine_host_role(host_name: str, host_taxid: str, host_order: str,
                        is_lab_host: bool, isolation_source: str) -> str:
    """Classify the host as reservoir, spillover or laboratory.

    Driven by host.lineage rank rather than substring matching. The previous version
    classified any record whose serialised JSON merely contained 'serum', 'sero' or
    'igg' as a spillover host, which caught unrelated submitter and method text.
    """
    if is_lab_host:
        return "laboratory_host"
    if host_taxid == "9606" or host_name == "Homo sapiens":
        return "spillover_host_human"
    if host_order:
        return HOST_ORDER_ROLE.get(host_order) or f"other_host_{host_order.lower()}"
    if host_name:
        return "host_unranked"
    if isolation_source:
        return "environmental_or_unspecified"
    return "unknown"


# --------------------------------------------------------------------------------------
# Grouping
# --------------------------------------------------------------------------------------


def build_group_key(
    species_name: str,
    biosample: str,
    isolate: str,
    strain: str,
    collection_date: str,
    location: str,
    accession: str,
) -> Tuple[str, str, float]:
    """Group accessions that represent the same physical virus isolate.

    BioSample is present on only 87 of 15,991 records here (0.54%), so a BioSample-first
    cascade almost never fires and the dataset shatters into per-accession groups.
    ``isolate.name`` is present on 15,491 records (96.9%) and is the field submitters
    actually use to tie an isolate's S, M and L segments together -- grouping on it
    recovers 1,484 complete S/M/L triplets against 214 previously.
    """
    label = isolate or strain
    if species_name and biosample and label:
        return (f"t1|{species_name}|{biosample}|{normalize_key_text(label)}",
                "tier_1_biosample_isolate", 1.00)
    if species_name and isolate:
        # Normalized, not raw. Tier 2 claims 15,404 of 15,991 records, so it short-circuits
        # the normalized tier below for every record that has an isolate name -- which is
        # to say normalize_key_text never reached a grouping decision at all. 52 isolates
        # were split across spelling variants of one label ('Aa19-152/MInION' vs
        # 'Aa19-152/MinION', 'PUUV/Kazan/human/RT500/2019' vs 'PUUV/Kazan/human_RT500/2019').
        return f"t2|{species_name}|{normalize_key_text(isolate)}", "tier_2_species_isolate", 0.95
    if species_name and label:
        return f"t3|{species_name}|{normalize_key_text(label)}", "tier_3_species_isolate_normalized", 0.90
    if species_name and strain and collection_date and location:
        return f"t4|{species_name}|{strain}|{collection_date}|{location}", "tier_4_strain_date_loc", 0.75
    if species_name and strain and location:
        return f"t5|{species_name}|{strain}|{location}", "tier_5_strain_loc", 0.65
    return f"acc|{accession}", "fallback_accession_only", 0.10


# --------------------------------------------------------------------------------------
# Accession-level transform
# --------------------------------------------------------------------------------------


def transform_one_record(
    record: Dict[str, Any],
    annotation_index: Dict[str, Dict[str, Any]],
    biosample_index: Dict[str, BioSampleInfo],
    genomic_index: Dict[str, FastaRecord],
    protein_index: Dict[str, FastaRecord],
) -> Tuple[Optional[AccessionRecord], Dict[str, Any]]:
    accession = clean_text(record.get("accession")) or ""
    try:
        acc_base = accession_base(accession)
        nucleotide = record.get("nucleotide") or {}
        sequence_hash = normalize_sequence_hash(nucleotide.get("sequenceHash"))

        species = resolve_species(record)

        annotation = annotation_index.get(accession) or annotation_index.get(acc_base)
        products, protein_accessions, protein_hits = extract_annotation_details(annotation, protein_index)

        genomic = genomic_index.get(accession) or genomic_index.get(acc_base)
        nt_len = int(record.get("length") or 0) or (genomic.length if genomic else 0)
        defline = genomic.header if genomic else ""

        segment, segment_source, segment_confidence, segment_warnings = choose_segment(
            record.get("segment"), products, defline, nt_len)
        completeness, completeness_source, is_coding_complete = determine_completeness(record, segment, nt_len)

        protein_aa, protein_product, primary_protein_acc = select_primary_protein(protein_hits, segment)
        has_nss = any(hit.role == "NSs" for hit in protein_hits)

        isolate_block = record.get("isolate") or {}
        isolate = clean_text(isolate_block.get("name")) or ""
        isolation_source = clean_text(isolate_block.get("source")) or ""
        collection_date = clean_text(isolate_block.get("collectionDate")) or ""

        location_block = record.get("location") or {}
        geographic = clean_text(location_block.get("geographicLocation")) or ""
        geographic_region = clean_text(location_block.get("geographicRegion")) or ""

        biosample_acc = clean_text(record.get("biosample")) or ""
        biosample_info = biosample_index.get(biosample_acc.upper()) if biosample_acc else None

        strain = ""
        if biosample_info:
            strain = biosample_info.strain or ""
            isolate = isolate or (biosample_info.isolate or "")
            collection_date = collection_date or (biosample_info.collection_date or "")
            geographic = geographic or (biosample_info.geographic_location or "")
            isolation_source = isolation_source or (biosample_info.isolation_source or "")

        country, sub_region = split_geo(geographic)
        region = sub_region or geographic_region

        host_name, host_taxid, host_order, host_family = resolve_host_taxonomy(record)
        host_source = "data_report" if host_name else "missing"
        if not host_name and biosample_info and biosample_info.host:
            host_name = biosample_info.host
            host_taxid = biosample_info.host_taxid or ""
            host_source = "biosample"

        is_lab_host = bool(record.get("isLabHost"))
        lab_host = clean_text(record.get("labHost")) or ""
        is_vaccine_strain = bool(record.get("isVaccineStrain"))
        host_role = determine_host_role(host_name, host_taxid, host_order, is_lab_host, isolation_source)

        tissue = clean_text((biosample_info.tissue if biosample_info else None)) or ""
        bioprojects = record.get("bioprojects") or []
        bioproject = clean_text(bioprojects[0]) if isinstance(bioprojects, list) and bioprojects else ""
        sra_accessions = [str(x) for x in (record.get("sraAccessions") or []) if x]

        group_key, grouping_tier, group_confidence = build_group_key(
            species_name=species.canonical_name,
            biosample=biosample_acc,
            isolate=isolate,
            strain=strain,
            collection_date=collection_date,
            location=geographic,
            accession=accession,
        )

        anomaly = length_anomaly(segment, nt_len, completeness)
        warning_flags = uniq([
            *segment_warnings,
            *(["missing_host"] if not host_name else []),
            *(["missing_isolate_and_strain"] if not isolate and not strain else []),
            *(["lab_host"] if is_lab_host else []),
            *(["vaccine_strain"] if is_vaccine_strain else []),
            *(["missing_genomic_sequence"] if genomic is None else []),
            *(["no_annotation"] if annotation is None else []),
            *(["missing_protein_link"] if not protein_accessions else []),
            *([anomaly] if anomaly else []),
            *(["species_unclassified"] if species.status != "classified" else []),
            *(["missing_collection_date"] if not collection_date else []),
            *(["missing_location"] if not geographic else []),
        ])

        result = AccessionRecord(
            accession=accession,
            sequence_hash=sequence_hash,
            species_name=species.canonical_name,
            species_taxid=species.canonical_taxid,
            species_status=species.status,
            observed_organism=species.observed_name,
            observed_taxid=species.observed_taxid,
            ranks=species.ranks,
            segment=segment,
            segment_source=segment_source,
            segment_confidence=segment_confidence,
            completeness=completeness,
            completeness_source=completeness_source,
            is_coding_complete=is_coding_complete,
            nt_length=nt_len,
            biosample=biosample_acc,
            bioproject=bioproject or "",
            sra_accessions=sra_accessions,
            isolate=isolate,
            strain=strain,
            collection_date=collection_date,
            year=extract_year(collection_date),
            location=geographic,
            country=country,
            region=region,
            host_name=host_name,
            host_taxid=host_taxid,
            host_order=host_order,
            host_family=host_family,
            host_source=host_source,
            host_role=host_role,
            isolation_source=isolation_source,
            tissue_source=tissue,
            is_lab_host=is_lab_host,
            lab_host=lab_host,
            is_vaccine_strain=is_vaccine_strain,
            protein_accessions=uniq(protein_accessions + ([primary_protein_acc] if primary_protein_acc else [])),
            annotation_products=products,
            protein_aa=protein_aa,
            protein_product=protein_product,
            has_nss_orf=has_nss,
            group_key=group_key,
            grouping_tier=grouping_tier,
            group_confidence=group_confidence,
            is_refseq=(clean_text(record.get("sourceDatabase")) or "").upper() == "REFSEQ"
            or accession.startswith("NC_"),
            source_database=clean_text(record.get("sourceDatabase")) or "",
            warning_flags=warning_flags,
        )

        qc_row = {
            "accession": accession,
            "group_key": group_key,
            "grouping_tier": grouping_tier,
            "group_confidence": group_confidence,
            "species_name": species.canonical_name,
            "segment": segment,
            "segment_source": segment_source,
            "completeness": completeness,
            "nt_length": nt_len,
            "transform_ok": True,
            "warnings": ";".join(warning_flags),
        }
        return result, qc_row
    except Exception as exc:  # noqa: BLE001 - one bad record must not kill a 16k-record run
        logging.exception("Failed to transform accession %s", accession or "<unknown>")
        return None, {
            "accession": accession,
            "group_key": "",
            "grouping_tier": "",
            "group_confidence": 0.0,
            "species_name": "",
            "segment": "",
            "segment_source": "",
            "completeness": "",
            "nt_length": 0,
            "transform_ok": False,
            "warnings": f"transform_error:{type(exc).__name__}:{exc}",
            "traceback": traceback.format_exc(limit=5),
        }


# --------------------------------------------------------------------------------------
# Isolate reconstruction
# --------------------------------------------------------------------------------------


def choose_best_segment_record(records: Sequence[AccessionRecord]) -> Optional[AccessionRecord]:
    """Prefer coding-complete, then RefSeq, then longest, then explicit segment evidence."""
    if not records:
        return None
    return max(records, key=lambda r: (
        r.is_coding_complete,
        r.completeness == "complete",
        r.is_refseq,
        r.nt_length,
        r.segment_confidence,
    ))


def deduplicate_by_sequence_hash(records: Sequence[AccessionRecord]) -> List[AccessionRecord]:
    """Collapse records whose NCBI sequenceHash is identical.

    12,598 distinct hashes cover 15,991 records here, so 3,393 are exact resubmissions of
    a sequence already present. Counting them separately inflates group sizes and would
    duplicate rows into the ML tables.
    """
    seen: Dict[str, AccessionRecord] = {}
    passthrough: List[AccessionRecord] = []
    for record in records:
        if not record.sequence_hash:
            passthrough.append(record)
            continue
        existing = seen.get(record.sequence_hash)
        if existing is None:
            seen[record.sequence_hash] = record
        else:
            seen[record.sequence_hash] = choose_best_segment_record([existing, record]) or existing
    return list(seen.values()) + passthrough


def build_isolate_records(accessions: Sequence[AccessionRecord]) -> List[IsolateRecord]:
    groups: DefaultDict[str, List[AccessionRecord]] = defaultdict(list)
    for record in accessions:
        groups[record.group_key].append(record)

    isolates: List[IsolateRecord] = []
    for key, members in groups.items():
        members = deduplicate_by_sequence_hash(members)
        per_segment = {
            seg: choose_best_segment_record([r for r in members if r.segment == seg])
            for seg in ("S", "M", "L")
        }
        if not any(per_segment.values()):
            continue

        present = [seg for seg in ("S", "M", "L") if per_segment[seg]]
        coding_complete = [seg for seg in ("S", "M", "L")
                           if per_segment[seg] and per_segment[seg].is_coding_complete]

        if len(coding_complete) == 3:
            triplet_status = "coding_complete_triplet"
        elif len(present) == 3:
            triplet_status = "all_segments_present_not_all_complete"
        elif len(coding_complete) >= 1:
            triplet_status = "partial_triplet_with_complete_segment"
        else:
            triplet_status = "incomplete_triplet"

        ordered = [per_segment[seg] for seg in ("S", "M", "L") if per_segment[seg]]

        def pick(attr: str) -> Any:
            for record in ordered:
                value = getattr(record, attr)
                if value:
                    return value
            return "" if isinstance(getattr(ordered[0], attr), str) else getattr(ordered[0], attr)

        taxid_counts = Counter(r.species_taxid for r in ordered if r.species_taxid)
        species_names = {r.species_name for r in ordered if r.species_name}

        # Different collection dates or locations across a group's segments mean the
        # grouping key merged records that are not obviously the same physical sample.
        dates = {r.collection_date for r in ordered if r.collection_date}
        locations = {r.location for r in ordered if r.location}
        metadata_conflict = len(dates) > 1 or len(locations) > 1

        s_rec, m_rec, l_rec = per_segment["S"], per_segment["M"], per_segment["L"]
        protein_parts = [r.protein_aa for r in (s_rec, m_rec) if r and r.protein_aa]

        isolates.append(IsolateRecord(
            isolate_key=key,
            species_name=pick("species_name"),
            species_taxid=mode(taxid_counts),
            species_status=pick("species_status"),
            biosample=pick("biosample"),
            isolate=pick("isolate"),
            strain=pick("strain"),
            host_name=pick("host_name"),
            host_taxid=pick("host_taxid"),
            host_order=pick("host_order"),
            host_family=pick("host_family"),
            host_role=pick("host_role"),
            collection_date=pick("collection_date"),
            year=pick("year"),
            location=pick("location"),
            country=pick("country"),
            region=pick("region"),
            grouping_tier=ordered[0].grouping_tier,
            grouping_confidence=max(r.group_confidence for r in ordered),
            s_accession=s_rec.accession if s_rec else "",
            m_accession=m_rec.accession if m_rec else "",
            l_accession=l_rec.accession if l_rec else "",
            s_completeness=s_rec.completeness if s_rec else "absent",
            m_completeness=m_rec.completeness if m_rec else "absent",
            l_completeness=l_rec.completeness if l_rec else "absent",
            s_nt_length=s_rec.nt_length if s_rec else 0,
            m_nt_length=m_rec.nt_length if m_rec else 0,
            l_nt_length=l_rec.nt_length if l_rec else 0,
            s_protein_aa=s_rec.protein_aa if s_rec else "",
            m_protein_aa=m_rec.protein_aa if m_rec else "",
            l_protein_aa=l_rec.protein_aa if l_rec else "",
            concatenated_protein_aa="".join(protein_parts),
            segments_present=present,
            coding_complete_segments=coding_complete,
            triplet_status=triplet_status,
            has_complete_triplet=len(coding_complete) == 3,
            is_lab_derived=any(r.is_lab_host for r in ordered),
            is_vaccine_strain=any(r.is_vaccine_strain for r in ordered),
            has_nss_orf=any(r.has_nss_orf for r in ordered),
            metadata_conflict=metadata_conflict,
            # Within a group this is always False by construction (the group key contains
            # the species), so it is filled in afterwards by flag_cross_species_labels.
            taxonomic_incongruence=len(species_names) > 1,
            cross_species_partners=[],
            kmer_cluster_id="",
            s_protein_accessions=s_rec.protein_accessions if s_rec else [],
            m_protein_accessions=m_rec.protein_accessions if m_rec else [],
            l_protein_accessions=l_rec.protein_accessions if l_rec else [],
            all_accessions=sorted({r.accession for r in members}),
            warnings=uniq([w for r in ordered for w in r.warning_flags]),
        ))

    flag_cross_species_labels(isolates)
    return isolates


def flag_cross_species_labels(isolates: Sequence[IsolateRecord]) -> None:
    """Flag one isolate label whose segments were assigned to more than one species.

    This is the trustworthy cross-species signal in the package, and the closest thing
    here to the documented heterotypic reassortment case (1988 rat isolates carrying an
    HTNV S segment with SEOV M and L segments -- https://pmc.ncbi.nlm.nih.gov/articles/PMC6153690/).

    It has to be computed across groups rather than within one: the grouping key contains
    the canonical species, so every member of a group necessarily agrees on species and
    a within-group check can never fire.

    Only labels distinctive enough to mean the same sample in two different submissions
    are compared -- see is_cross_study_label.
    """
    by_label: DefaultDict[str, List[IsolateRecord]] = defaultdict(list)
    for isolate in isolates:
        label = isolate.isolate or isolate.strain
        if is_cross_study_label(label):
            by_label[normalize_key_text(label)].append(isolate)

    for members in by_label.values():
        species = {m.species_name for m in members if m.species_name}
        if len(species) < 2:
            continue
        for isolate in members:
            partners = sorted(species - {isolate.species_name})
            isolate.taxonomic_incongruence = True
            isolate.cross_species_partners = partners
            isolate.warnings = uniq(isolate.warnings + ["same_label_assigned_to_multiple_species"])


# --------------------------------------------------------------------------------------
# Alignment-free sequence comparison
# --------------------------------------------------------------------------------------
#
# NOTE ON RIGOUR: the ICTV species criterion for Hantaviridae is DEmARC on concatenated
# S/M/L expression products with a pairwise evolutionary distance cutoff of 0.1 under the
# WAG substitution model (https://ictv.global/report/chapter/hantaviridae/hantaviridae).
# That requires multiple sequence alignment and maximum-likelihood inference, which this
# script deliberately does not attempt -- it has no alignment dependency.
#
# What follows is an alignment-free k-mer (Mash-style) proxy. It is good enough to rank
# neighbours and triage reassortment candidates, and it is NOT a PED value. Every column
# it produces is named kmer_* rather than ped_* so the distinction survives into the
# output, and ortho_summary.json restates the limitation.


def kmer_profile(sequence: str, k: int = KMER_K) -> Dict[str, int]:
    if not sequence or len(sequence) < k:
        return {}
    return Counter(sequence[i:i + k] for i in range(len(sequence) - k + 1))


def profiles_to_matrix(profiles: Sequence[Dict[str, int]]) -> Tuple[Any, List[str]]:
    """Stack k-mer profiles into an L1-normalised numpy matrix over the observed vocabulary."""
    vocabulary = sorted({kmer for profile in profiles for kmer in profile})
    position = {kmer: i for i, kmer in enumerate(vocabulary)}

    # Fill via one fancy-indexed assignment rather than a scalar write per cell. These
    # profiles are sparse (~6% dense for the full ML table), so the per-cell loop was
    # doing millions of individual numpy writes to place a few million values.
    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []
    for row, profile in enumerate(profiles):
        rows.extend([row] * len(profile))
        cols.extend(position[kmer] for kmer in profile)
        vals.extend(profile.values())

    matrix = np.zeros((len(profiles), len(vocabulary)), dtype=np.float64)
    if vals:
        matrix[np.fromiter(rows, dtype=np.intp, count=len(rows)),
               np.fromiter(cols, dtype=np.intp, count=len(cols))] = np.fromiter(
                   vals, dtype=np.float64, count=len(vals))
    totals = matrix.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    return matrix / totals, vocabulary


def cosine_distance_matrix(matrix: Any) -> Any:
    """Pairwise cosine distance in [0, 1]."""
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    unit = matrix / norms
    similarity = np.clip(unit @ unit.T, -1.0, 1.0)
    return 1.0 - similarity


def single_linkage_clusters(distances: Any, threshold: float) -> List[int]:
    """Union-find single-linkage clustering on a dense distance matrix."""
    size = distances.shape[0]
    parent = list(range(size))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    rows, cols = np.where(distances <= threshold)
    for a, b in zip(rows.tolist(), cols.tolist()):
        if a < b:
            union(a, b)

    labels: Dict[int, int] = {}
    out: List[int] = []
    for index in range(size):
        root = find(index)
        if root not in labels:
            labels[root] = len(labels)
        out.append(labels[root])
    return out


def assign_kmer_clusters(isolates: Sequence[IsolateRecord], threshold: float) -> None:
    """Cluster isolates within each species by concatenated N+GPC k-mer distance.

    Replaces the previous cluster_species_by_ped, which computed a Hamming distance
    between *unaligned* protein strings padded by their length difference. That is not a
    p-distance and the clusters it produced were not interpretable.
    """
    by_species: DefaultDict[str, List[IsolateRecord]] = defaultdict(list)
    for isolate in isolates:
        if isolate.concatenated_protein_aa:
            by_species[isolate.species_name].append(isolate)

    for species, members in by_species.items():
        if len(members) == 1:
            members[0].kmer_cluster_id = f"{species}|c0"
            continue
        profiles = [kmer_profile(m.concatenated_protein_aa) for m in members]
        matrix, _vocab = profiles_to_matrix(profiles)
        distances = cosine_distance_matrix(matrix)
        for member, label in zip(members, single_linkage_clusters(distances, threshold)):
            member.kmer_cluster_id = f"{species}|c{label}"


# --------------------------------------------------------------------------------------
# Reassortment detection
# --------------------------------------------------------------------------------------


def _rank_positions(order: Dict[str, Any]) -> Dict[str, Any]:
    """Rank-lookup tables from per-segment neighbour orderings: position[segment][i, j]
    is j's rank in i's ordering. Shared by the k-mer (compute_segment_divergence) and
    tree-distance (tree_divergence_for_species) diagnostics, which differ only in how
    `order` was derived.
    """
    position: Dict[str, Any] = {}
    for segment in ("S", "M", "L"):
        ordering = order[segment]
        table = np.empty_like(ordering)
        rows = np.arange(ordering.shape[0])[:, None]
        table[rows, ordering] = np.arange(ordering.shape[1])[None, :]
        position[segment] = table
    return position


def _sl_m_rank_metrics(
    position: Dict[str, Any], span: int, index: int, s_nn: int, m_nn: int, l_nn: int
) -> Tuple[float, float]:
    """sl_rank_separation and m_rank_offset: continuous measurements of how far apart the
    per-segment nearest-neighbour orderings sit. NOT evidence of reassortment on their own
    -- see compute_segment_divergence's docstring for why the k-mer version of this signal
    was withdrawn as a classifier.
    """
    s_pos_l = float(position["L"][index, s_nn]) / span
    l_pos_s = float(position["S"][index, l_nn]) / span
    s_pos_m = float(position["S"][index, m_nn]) / span
    l_pos_m = float(position["L"][index, m_nn]) / span
    return max(s_pos_l, l_pos_s), min(s_pos_m, l_pos_m)




def compute_segment_divergence(
    isolates: Sequence[IsolateRecord],
    min_species_size: int = 8,
) -> List[SegmentDivergenceProfile]:
    """Per-segment nearest-neighbour diagnostics within each species.

    WITHDRAWN CLASSIFIER -- READ BEFORE USING THESE COLUMNS
    ------------------------------------------------------
    An earlier version of this function emitted reassortment calls
    (m_segment_reassortment_candidate / sl_cosegregation_violation) derived from
    discordance between the per-segment nearest-neighbour orderings computed here.
    Those calls were withdrawn because they failed validation against the literature
    (https://pmc.ncbi.nlm.nih.gov/articles/PMC6153690/), which reports reassortment as
    frequent in PUUV (19.1-32%), SNV, DOBV and TULV, and absent or sparse in ANDV.

    Measured on this dataset, the classifier produced the inverse ranking:

        species          n     flagged      literature
        ANDV            68      72.1%       absent or sparse
        CHOV             8     100.0%       not reported
        DOBV            62      38.7%       frequent
        PUUV           416       7.7%       frequent (19.1-32%)
        TULV           313       3.5%       frequent

    Sweeping the discordance threshold from 0.25 to 0.85 scaled all rates down but never
    changed the ordering. The cause is intrinsic to the metric: cosine distance on 3-mer
    amino-acid composition tracks within-species sequence diversity and sample size, not
    tree topology. Species represented by many near-identical outbreak sequences get
    stable neighbour orderings and few flags; species represented by a handful of
    divergent, geographically scattered sequences get noisy orderings and many flags.
    The 100% rate for CHOV at n=8 is the clearest symptom.

    Genuine reassortment detection requires testing phylogenetic incongruence between
    per-segment trees built from aligned sequences (e.g. MAFFT + IQ-TREE), which this
    script has no dependency for. So this function now emits only the underlying
    distances and rank separations as continuous diagnostics, and marks every assessable
    isolate 'requires_phylogenetic_confirmation'.

    The one reassortment-adjacent signal that IS trustworthy here is
    taxonomic_incongruence: NCBI assigning an isolate's segments to different species.
    That is a metadata fact, not an inference, and it is reported separately.
    """
    profiles: List[SegmentDivergenceProfile] = []
    by_species: DefaultDict[str, List[IsolateRecord]] = defaultdict(list)

    def unassessable(isolate: IsolateRecord, species: str, reason: str) -> SegmentDivergenceProfile:
        return SegmentDivergenceProfile(
            isolate_key=isolate.isolate_key, species_name=species, assessable=False,
            reason=reason, s_neighbour="", m_neighbour="", l_neighbour="",
            s_distance=0.0, m_distance=0.0, l_distance=0.0,
            sl_rank_separation=0.0, m_rank_offset=0.0,
            taxonomic_incongruence=isolate.taxonomic_incongruence,
            reassortment_status="not_assessable",
            reported_species_prior=REPORTED_REASSORTMENT.get(species, "not reported"))

    for isolate in isolates:
        if not (isolate.s_protein_aa and isolate.m_protein_aa and isolate.l_protein_aa):
            profiles.append(unassessable(isolate, isolate.species_name,
                                         "requires S, M and L protein sequences"))
            continue
        by_species[isolate.species_name].append(isolate)

    for species, members in by_species.items():
        prior = REPORTED_REASSORTMENT.get(species, "not reported")
        if len(members) < min_species_size:
            for isolate in members:
                profiles.append(unassessable(
                    isolate, species,
                    f"only {len(members)} assessable isolate(s) in this species "
                    f"(need {min_species_size})"))
            continue

        neighbour: Dict[str, List[int]] = {}
        distance: Dict[str, List[float]] = {}
        ranking: Dict[str, Any] = {}
        for segment, attribute in (("S", "s_protein_aa"), ("M", "m_protein_aa"), ("L", "l_protein_aa")):
            matrix, _vocab = profiles_to_matrix([kmer_profile(getattr(m, attribute)) for m in members])
            distances = cosine_distance_matrix(matrix)
            np.fill_diagonal(distances, np.inf)
            order = np.argsort(distances, axis=1)
            ranking[segment] = order
            neighbour[segment] = order[:, 0].tolist()
            distance[segment] = [float(distances[i, order[i, 0]]) for i in range(len(members))]

        position = _rank_positions(ranking)
        span = max(1, len(members) - 1)

        for index, isolate in enumerate(members):
            s_nn, m_nn, l_nn = neighbour["S"][index], neighbour["M"][index], neighbour["L"][index]
            sl_separation, m_offset = _sl_m_rank_metrics(position, span, index, s_nn, m_nn, l_nn)

            profiles.append(SegmentDivergenceProfile(
                isolate_key=isolate.isolate_key,
                species_name=species,
                assessable=True,
                reason="",
                s_neighbour=members[s_nn].isolate_key,
                m_neighbour=members[m_nn].isolate_key,
                l_neighbour=members[l_nn].isolate_key,
                s_distance=round(distance["S"][index], 6),
                m_distance=round(distance["M"][index], 6),
                l_distance=round(distance["L"][index], 6),
                sl_rank_separation=round(sl_separation, 4),
                m_rank_offset=round(m_offset, 4),
                taxonomic_incongruence=isolate.taxonomic_incongruence,
                reassortment_status="requires_phylogenetic_confirmation",
                reported_species_prior=prior,
            ))
    return profiles


# --------------------------------------------------------------------------------------
# Pathogenicity and surveillance risk
# --------------------------------------------------------------------------------------


@lru_cache(maxsize=None)
def species_risk_profile(species_name: str) -> Dict[str, Any]:
    """Look up curated risk metadata. Unlisted species are 'unknown', never 'safe'.

    Cached: this is a pure function over ~39 distinct species but is called once per row
    in every table builder, ~91k times per run. Callers must treat the returned dict as
    read-only.
    """
    entry = SPECIES_RISK.get(species_name)
    if entry is None:
        return {
            "virus_abbrev": "",
            "disease_syndrome": "unknown",
            "cfr_band": -1,
            "cfr_band_label": CFR_BAND_LABEL[-1],
            "cfr_reported_range": "",
            "reservoir_host": "",
            "endemic_region": "",
            "pathogenicity_class": "unknown",
            "risk_citation": "",
            "reservoir_order": "",
        }
    band = entry["cfr_band"]
    syndrome = entry["syndrome"]
    pathogenicity = SYNDROME_TO_PATHOGENICITY.get(syndrome, "unknown")
    return {
        "virus_abbrev": entry["abbrev"],
        "disease_syndrome": syndrome,
        "cfr_band": band,
        "cfr_band_label": CFR_BAND_LABEL.get(band, CFR_BAND_LABEL[-1]),
        "cfr_reported_range": entry["cfr_reported"],
        "reservoir_host": entry["reservoir"],
        "endemic_region": entry["region"],
        "pathogenicity_class": pathogenicity,
        "risk_citation": entry["citation"],
        "reservoir_order": RESERVOIR_GENUS_ORDER.get(entry["reservoir"].split()[0], ""),
    }


def host_jump_status(observed_order: str, reservoir_order: str) -> str:
    """Compare the sampled host's order against the species' expected reservoir order.

    A pathogenic rodent-borne virus recovered from a shrew or bat is a documented and
    genuinely interesting observation (SEOV in Chinese shrews and bats, ANDV-clade
    sequences in Brazilian bats both appear in this package), not a data error.
    """
    if not observed_order or not reservoir_order:
        return "unknown"
    if observed_order == reservoir_order:
        return "expected_reservoir_order"
    if observed_order == "Primates":
        return "human_or_primate_spillover"
    return "host_order_mismatch"


def surveillance_priority(
    isolate: IsolateRecord,
    risk: Dict[str, Any],
    profile: Optional[SegmentDivergenceProfile],
    phylo_call: Optional[str] = None,
) -> Dict[str, Any]:
    """Transparent additive priority score.

    Every term is emitted as its own column so the total can be audited, disagreed with
    or re-weighted downstream. This is a triage aid for surveillance, not a validated
    epidemiological risk model.

    Two different kinds of reassortment evidence exist in this pipeline, and only one of
    them contributes here:

    - The k-mer segment-divergence diagnostics (compute_segment_divergence, `profile`
      here) do NOT contribute a term, and never will under this function. An earlier
      version weighted them as reassortment evidence; that classifier failed validation
      against the literature and feeding it into a risk score would have propagated an
      artefact of sample size into a ranking meant to guide surveillance.
    - `phylo_call`, the bootstrap-gated tree-topology call from --phylogenetics
      (classify_reassortment_call), is a different kind of evidence -- real per-segment
      alignments and maximum-likelihood trees, gated on ultrafast-bootstrap support
      specifically so an unsupported pattern is never reported as a finding. That
      reasoning doesn't disqualify it the way it disqualified the k-mer diagnostics; it's
      the reason to trust it. Only 'supported_reassortment_candidate' contributes --
      never 'weak_signal_low_bootstrap_support' or 'no_conflict', and never at all when
      --phylogenetics wasn't run (phylo_call is None, the default).

    taxonomic_incongruence -- a metadata fact, not an inference -- has always contributed
    and still does, independent of either kind of reassortment evidence above.
    """
    weights = SURVEILLANCE_WEIGHTS
    cfr_term = weights["cfr_by_band"].get(int(risk["cfr_band"]), 0.5)
    spillover_term = weights["human_spillover"] if isolate.host_role == "spillover_host_human" else 0.0
    incongruence_term = weights["taxonomic_incongruence"] if isolate.taxonomic_incongruence else 0.0
    jump_term = weights["host_jump"] if host_jump_status(
        isolate.host_order, risk["reservoir_order"]) == "host_order_mismatch" else 0.0
    completeness_term = weights["complete_triplet"] if isolate.has_complete_triplet else 0.0
    nss_term = weights["nss_orf"] if isolate.has_nss_orf else 0.0
    reassortment_term = weights["supported_reassortment"] if phylo_call == "supported_reassortment_candidate" else 0.0
    total = (cfr_term + spillover_term + incongruence_term + jump_term
             + completeness_term + nss_term + reassortment_term)
    return {
        "priority_cfr_term": cfr_term,
        "priority_human_spillover_term": spillover_term,
        "priority_taxonomic_incongruence_term": incongruence_term,
        "priority_host_jump_term": jump_term,
        "priority_triplet_term": completeness_term,
        "priority_nss_term": nss_term,
        "priority_reassortment_term": reassortment_term,
        "surveillance_priority_score": round(total, 3),
    }




# --------------------------------------------------------------------------------------
# Output tables
# --------------------------------------------------------------------------------------


def build_accession_rows(accessions: Sequence[AccessionRecord]) -> List[Dict[str, Any]]:
    rows = []
    for record in accessions:
        risk = species_risk_profile(record.species_name)
        rows.append({
            "accession": record.accession,
            "sequence_hash": record.sequence_hash,
            "species_name": record.species_name,
            "species_taxid": record.species_taxid,
            "species_status": record.species_status,
            "observed_organism": record.observed_organism,
            "observed_taxid": record.observed_taxid,
            "genus": record.ranks.get("genus", ""),
            "subfamily": record.ranks.get("subfamily", ""),
            "family": record.ranks.get("family", ""),
            "order": record.ranks.get("order", ""),
            "class": record.ranks.get("class", ""),
            "phylum": record.ranks.get("phylum", ""),
            "virus_abbrev": risk["virus_abbrev"],
            "pathogenicity_class": risk["pathogenicity_class"],
            "expected_reservoir_order": risk["reservoir_order"],
            "host_jump_status": host_jump_status(record.host_order, risk["reservoir_order"]),
            "segment": record.segment,
            "segment_source": record.segment_source,
            "segment_confidence": record.segment_confidence,
            "completeness": record.completeness,
            "completeness_source": record.completeness_source,
            "is_coding_complete": int(record.is_coding_complete),
            "nt_length": record.nt_length,
            "biosample": record.biosample,
            "bioproject": record.bioproject,
            "sra_accessions": ";".join(record.sra_accessions),
            "isolate": record.isolate,
            "strain": record.strain,
            "collection_date": record.collection_date,
            "year": record.year,
            "location": record.location,
            "country": record.country,
            "region": record.region,
            "host_name": record.host_name,
            "host_taxid": record.host_taxid,
            "host_order": record.host_order,
            "host_family": record.host_family,
            "host_source": record.host_source,
            "host_role": record.host_role,
            "isolation_source": record.isolation_source,
            "tissue_source": record.tissue_source,
            "is_lab_host": int(record.is_lab_host),
            "lab_host": record.lab_host,
            "is_vaccine_strain": int(record.is_vaccine_strain),
            "protein_accessions": ";".join(record.protein_accessions),
            "annotation_products": ";".join(record.annotation_products),
            "protein_product": record.protein_product,
            "has_nss_orf": int(record.has_nss_orf),
            "group_key": record.group_key,
            "grouping_tier": record.grouping_tier,
            "group_confidence": record.group_confidence,
            "is_refseq": int(record.is_refseq),
            "source_database": record.source_database,
            "warning_flags": ";".join(record.warning_flags),
        })
    return rows


def build_species_rows(accessions: Sequence[AccessionRecord],
                       isolates: Sequence[IsolateRecord]) -> List[Dict[str, Any]]:
    """One row per canonical ICTV species -- the top of the aim-2 grouping hierarchy."""
    by_species: DefaultDict[str, List[AccessionRecord]] = defaultdict(list)
    for record in accessions:
        by_species[record.species_name].append(record)
    isolates_by_species: DefaultDict[str, List[IsolateRecord]] = defaultdict(list)
    for isolate in isolates:
        isolates_by_species[isolate.species_name].append(isolate)

    rows = []
    for species, members in sorted(by_species.items()):
        species_isolates = isolates_by_species.get(species, [])
        risk = species_risk_profile(species)
        segments = Counter(r.segment for r in members if r.segment)
        hosts = Counter(r.host_name for r in members if r.host_name)
        countries = Counter(r.country for r in members if r.country)
        years = sorted({r.year for r in members if r.year})
        rows.append({
            "species_name": species,
            "species_status": members[0].species_status,
            "observed_taxids": ";".join(sorted({r.observed_taxid for r in members if r.observed_taxid})),
            "canonical_taxids": ";".join(sorted({r.species_taxid for r in members if r.species_taxid})),
            "legacy_synonyms": ";".join(sorted({r.observed_organism for r in members
                                                if r.observed_organism and r.observed_organism != species})),
            "subfamily": members[0].ranks.get("subfamily", ""),
            "family": members[0].ranks.get("family", ""),
            "order": members[0].ranks.get("order", ""),
            **{k: v for k, v in risk.items()},
            "reported_reassortment": REPORTED_REASSORTMENT.get(species, "not reported"),
            "accession_count": len(members),
            "isolate_count": len(species_isolates),
            "complete_triplet_count": sum(1 for i in species_isolates if i.has_complete_triplet),
            "s_count": segments.get("S", 0),
            "m_count": segments.get("M", 0),
            "l_count": segments.get("L", 0),
            "distinct_hosts": len(hosts),
            "top_host": mode(hosts),
            "human_derived_count": sum(1 for r in members if r.host_taxid == "9606"),
            "distinct_countries": len(countries),
            "top_country": mode(countries),
            "year_min": years[0] if years else "",
            "year_max": years[-1] if years else "",
        })
    return rows


def build_isolate_rows(isolates: Sequence[IsolateRecord],
                       divergence: Dict[str, SegmentDivergenceProfile]) -> List[Dict[str, Any]]:
    rows = []
    for isolate in isolates:
        risk = species_risk_profile(isolate.species_name)
        profile = divergence.get(isolate.isolate_key)
        rows.append({
            "isolate_key": isolate.isolate_key,
            "species_name": isolate.species_name,
            "species_taxid": isolate.species_taxid,
            "species_status": isolate.species_status,
            "virus_abbrev": risk["virus_abbrev"],
            "pathogenicity_class": risk["pathogenicity_class"],
            "kmer_cluster_id": isolate.kmer_cluster_id,
            "biosample": isolate.biosample,
            "isolate": isolate.isolate,
            "strain": isolate.strain,
            "host_name": isolate.host_name,
            "host_taxid": isolate.host_taxid,
            "host_order": isolate.host_order,
            "host_family": isolate.host_family,
            "host_role": isolate.host_role,
            "collection_date": isolate.collection_date,
            "year": isolate.year,
            "location": isolate.location,
            "country": isolate.country,
            "region": isolate.region,
            "grouping_tier": isolate.grouping_tier,
            "grouping_confidence": isolate.grouping_confidence,
            "s_accession": isolate.s_accession,
            "m_accession": isolate.m_accession,
            "l_accession": isolate.l_accession,
            "s_completeness": isolate.s_completeness,
            "m_completeness": isolate.m_completeness,
            "l_completeness": isolate.l_completeness,
            "s_nt_length": isolate.s_nt_length,
            "m_nt_length": isolate.m_nt_length,
            "l_nt_length": isolate.l_nt_length,
            "segments_present": ";".join(isolate.segments_present),
            "coding_complete_segments": ";".join(isolate.coding_complete_segments),
            "triplet_status": isolate.triplet_status,
            "has_complete_triplet": int(isolate.has_complete_triplet),
            "concatenated_protein_aa": isolate.concatenated_protein_aa,
            "concatenated_protein_length": len(isolate.concatenated_protein_aa),
            "is_lab_derived": int(isolate.is_lab_derived),
            "is_vaccine_strain": int(isolate.is_vaccine_strain),
            "has_nss_orf": int(isolate.has_nss_orf),
            "metadata_conflict": int(isolate.metadata_conflict),
            "taxonomic_incongruence": int(isolate.taxonomic_incongruence),
            "cross_species_partners": ";".join(isolate.cross_species_partners),
            "divergence_assessable": int(bool(profile and profile.assessable)),
            "reassortment_status": profile.reassortment_status if profile else "not_assessable",
            "sl_rank_separation": profile.sl_rank_separation if profile else 0.0,
            "m_rank_offset": profile.m_rank_offset if profile else 0.0,
            "s_protein_accessions": ";".join(isolate.s_protein_accessions),
            "m_protein_accessions": ";".join(isolate.m_protein_accessions),
            "l_protein_accessions": ";".join(isolate.l_protein_accessions),
            "all_accessions": ";".join(isolate.all_accessions),
            "warnings": ";".join(isolate.warnings),
        })
    return rows


def build_segment_divergence_rows(profiles: Sequence[SegmentDivergenceProfile]) -> List[Dict[str, Any]]:
    return [{
        "isolate_key": profile.isolate_key,
        "species_name": profile.species_name,
        "assessable": int(profile.assessable),
        "not_assessable_reason": profile.reason,
        "reassortment_status": profile.reassortment_status,
        "taxonomic_incongruence": int(profile.taxonomic_incongruence),
        "sl_rank_separation": profile.sl_rank_separation,
        "m_rank_offset": profile.m_rank_offset,
        "s_nearest_neighbour": profile.s_neighbour,
        "m_nearest_neighbour": profile.m_neighbour,
        "l_nearest_neighbour": profile.l_neighbour,
        "s_kmer_distance": profile.s_distance,
        "m_kmer_distance": profile.m_distance,
        "l_kmer_distance": profile.l_distance,
        "reported_species_prior": profile.reported_species_prior,
        "distance_metric": "cosine distance on 3-mer amino-acid profiles (alignment-free; NOT ICTV PED)",
        "interpretation_warning": "These are divergence diagnostics, NOT reassortment calls. "
                                 "A classifier built on them failed validation against the "
                                 "literature (inverted PUUV vs ANDV ranking) because the metric "
                                 "tracks within-species diversity and sample size rather than "
                                 "tree topology. Confirm any candidate with aligned per-segment "
                                 "phylogenies (MAFFT + IQ-TREE) before drawing conclusions.",
    } for profile in profiles]


def build_risk_rows(
    isolates: Sequence[IsolateRecord],
    divergence: Dict[str, SegmentDivergenceProfile],
) -> List[Dict[str, Any]]:
    """ortho_risk_table.csv is always phylo-blind now (Stage 2's phylogenetics analysis is a
    separate script/run; combine on isolate_key -- see README's Outputs table). This is
    byte-identical to the pre-split phylogenetics-free run's output.
    """
    rows = []
    for isolate in isolates:
        risk = species_risk_profile(isolate.species_name)
        profile = divergence.get(isolate.isolate_key)
        terms = surveillance_priority(isolate, risk, profile, None)
        rows.append({
            "isolate_key": isolate.isolate_key,
            "species_name": isolate.species_name,
            "virus_abbrev": risk["virus_abbrev"],
            "pathogenicity_class": risk["pathogenicity_class"],
            "disease_syndrome": risk["disease_syndrome"],
            "cfr_band": risk["cfr_band"],
            "cfr_band_label": risk["cfr_band_label"],
            "cfr_reported_range": risk["cfr_reported_range"],
            "reservoir_host": risk["reservoir_host"],
            "endemic_region": risk["endemic_region"],
            "risk_citation": risk["risk_citation"],
            "host_role": isolate.host_role,
            "host_order": isolate.host_order,
            "expected_reservoir_order": risk["reservoir_order"],
            "host_jump_status": host_jump_status(isolate.host_order, risk["reservoir_order"]),
            "observed_host": isolate.host_name,
            "country": isolate.country,
            "year": isolate.year,
            "triplet_status": isolate.triplet_status,
            "has_complete_triplet": int(isolate.has_complete_triplet),
            "has_nss_orf": int(isolate.has_nss_orf),
            "reassortment_status": profile.reassortment_status if profile else "not_assessable",
            "taxonomic_incongruence": int(isolate.taxonomic_incongruence),
            "cross_species_partners": ";".join(isolate.cross_species_partners),
            # "" (not None) when --phylogenetics wasn't run or didn't reach this isolate --
            # matches the string-typed convention every other column in this table uses;
            # surveillance_priority is the one place None vs "not analyzed" is worth
            # distinguishing (see ortho_reassortment_table.csv), not this summary table.
            "phylo_reassortment_call": "",
            **terms,
        })
    rows.sort(key=lambda r: r["surveillance_priority_score"], reverse=True)
    return rows


def build_host_rows(accessions: Sequence[AccessionRecord],
                    isolates: Sequence[IsolateRecord]) -> List[Dict[str, Any]]:
    grouped: DefaultDict[Tuple[str, str], List[AccessionRecord]] = defaultdict(list)
    for record in accessions:
        if record.host_name or record.host_taxid:
            grouped[(record.host_name, record.host_taxid)].append(record)
    isolate_counts: Counter = Counter((i.host_name, i.host_taxid) for i in isolates)

    rows = []
    for (name, taxid), members in grouped.items():
        species = Counter(r.species_name for r in members if r.species_name)
        rows.append({
            "host_name": name,
            "host_taxid": taxid,
            "host_order": members[0].host_order,
            "host_family": members[0].host_family,
            "host_role": mode(Counter(r.host_role for r in members)),
            "host_risk_prior": HOST_ORDER_RISK_PRIOR.get(members[0].host_order, "unknown"),
            "accession_count": len(members),
            "isolate_count": isolate_counts.get((name, taxid), 0),
            "distinct_virus_species": len(species),
            "top_virus_species": mode(species),
            "all_virus_species": ";".join(sorted(species)),
            "representative_biosample": next((r.biosample for r in members if r.biosample), ""),
            "representative_isolate": next((r.isolate for r in members if r.isolate), ""),
        })
    rows.sort(key=lambda r: r["accession_count"], reverse=True)
    return rows


def _species_effort_stats(
    accessions: Sequence[AccessionRecord], isolates: Sequence[IsolateRecord],
) -> Dict[str, Dict[str, Any]]:
    """Per-species surveillance-effort covariates for build_ml_tables -- see that
    function's docstring for why these matter as confound controls, not biology.

    Counts are over ALL accessions/isolates of a species, not just the ML-usable subset
    build_ml_tables filters down to: the point is to measure how much this species has
    been sampled overall, which is a property of the species' surveillance history, not
    of which isolates happened to have a usable concatenated protein.

    effort_years_since_first is the SPAN between the earliest and latest recorded year
    for that species (year_max - year_min), not years-since-first-to-today -- today's
    date isn't a meaningful reference point for a virus nobody has sequenced since 2004,
    and would conflate "how long ago sampling started" with "how long has it continued."
    """
    acc_by_species: DefaultDict[str, List[AccessionRecord]] = defaultdict(list)
    for record in accessions:
        acc_by_species[record.species_name].append(record)
    iso_by_species: DefaultDict[str, List[IsolateRecord]] = defaultdict(list)
    for isolate in isolates:
        iso_by_species[isolate.species_name].append(isolate)

    stats: Dict[str, Dict[str, Any]] = {}
    for species in set(acc_by_species) | set(iso_by_species):
        acc_members = acc_by_species.get(species, [])
        iso_members = iso_by_species.get(species, [])
        years = sorted(int(r.year) for r in acc_members if r.year.isdigit())
        countries = {r.country for r in acc_members if r.country}
        stats[species] = {
            "effort_n_accessions": len(acc_members),
            "effort_n_isolates": len(iso_members),
            "effort_n_countries": len(countries),
            "effort_n_years": len(set(years)),
            "effort_first_year": years[0] if years else "",
            "effort_years_since_first": (years[-1] - years[0]) if years else 0,
            "effort_n_complete_triplets": sum(1 for i in iso_members if i.has_complete_triplet),
        }
    return stats


_EFFORT_STATS_DEFAULT: Dict[str, Any] = {
    "effort_n_accessions": 0, "effort_n_isolates": 0, "effort_n_countries": 0,
    "effort_n_years": 0, "effort_first_year": "", "effort_years_since_first": 0,
    "effort_n_complete_triplets": 0,
}


def build_ml_tables(
    accessions: Sequence[AccessionRecord],
    isolates: Sequence[IsolateRecord],
    divergence: Dict[str, SegmentDivergenceProfile],
    max_features: int,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Build the ML feature matrix and the scaling artifact.

    Four deliberate changes from the previous implementation:

    - Features are emitted RAW. The previous version z-scored every column across the
      whole dataset before writing, which leaks test-fold statistics into training. The
      per-column mean/std are written to ortho_feature_scaling.json instead, so scaling
      can be fitted inside the caller's cross-validation fold.
    - Feature selection is unsupervised (top-variance k-mers). The previous selector
      scored features as abs(between)/variance against label indices produced by sorting
      taxids alphabetically, which imposes an arbitrary ordinal scale on unordered
      categories.
    - group_id and year are included so downstream CV can split by species and by time.
      Random splits over this data overstate accuracy badly, because segments from the
      same outbreak are near-identical.
    - effort_* columns (n_accessions, n_isolates, n_countries, n_years, first_year,
      years_since_first, n_complete_triplets -- see _species_effort_stats) are added per
      species. label_pathogenicity_class, label_cfr_band and label_disease_syndrome are
      pure functions of species_name, so "no known human disease" and "nobody has
      sequenced this reservoir" are literally the same label value. That confound
      correlates with country, year and research funding: a species judged harmless
      because it has 6 accessions from one 2004 paper is not evidence of anything. These
      are CONFOUND CONTROLS, not biology -- they exist so a model can be trained or
      ablated with sampling intensity explicitly accounted for, not so it can learn to
      predict pathogenicity from effort_n_accessions. A model whose held-out performance
      collapses when the effort_* block is dropped was substantially reading sampling
      intensity, not virology, and that result would otherwise have gone unnoticed.

    is_lab_derived is now FLAGGED, not filtered. An earlier version excluded every
    is_lab_derived isolate from this matrix, on the same "not a naturally occurring
    sample" reasoning --phylogenetics still uses. But NCBI's isLabHost field marks the
    immediate recorded host, not sequence provenance: of 164 is_lab_derived isolates, 21
    have host_taxid 9606 -- these are lab-acquired/occupational human infections in
    exactly the highest-CFR species (HTNV, ANDV, PUUV, DOBV) -- and most of the rest are
    samples from colony-maintained reservoir animals (Mus musculus, Rattus norvegicus,
    Apodemus agrarius), not attenuated or culture-adapted sequences. Blanket-filtering
    them discarded real positive evidence of human/reservoir infection for a
    pathogenicity-prediction task, exactly the risk this row's own citation warns about
    (https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-025-02059-0 --
    human-derived sequence similarity to reservoir sequences is itself the spillover-risk
    signal; excluding it removes the baseline). is_vaccine_strain isolates remain
    excluded: a vaccine strain is deliberately attenuated, so labelling its sequence with
    its parent species' wild-type pathogenicity class would inject a real, not
    hypothetical, false signal. --phylogenetics is unaffected by this: it keeps excluding
    both is_lab_derived and is_vaccine_strain, because a passage/attenuation history is a
    tree-topology risk (branch-length distortion) regardless of whether the sequence
    itself is informative for a similarity/classification feature.
    """
    usable = [i for i in isolates if i.concatenated_protein_aa and not i.is_vaccine_strain]
    if not usable:
        return [], {"features": [], "note": "no usable isolates"}

    # Computed over ALL accessions/isolates (not just `usable`) -- see
    # _species_effort_stats's docstring for why.
    effort_by_species = _species_effort_stats(accessions, isolates)

    profiles = [kmer_profile(i.concatenated_protein_aa) for i in usable]
    matrix, vocabulary = profiles_to_matrix(profiles)

    variances = matrix.var(axis=0)
    keep = np.argsort(variances)[::-1][:max_features]
    keep = sorted(keep.tolist())
    selected = [vocabulary[i] for i in keep]
    reduced = matrix[:, keep]

    rows: List[Dict[str, Any]] = []
    for position, isolate in enumerate(usable):
        risk = species_risk_profile(isolate.species_name)
        profile = divergence.get(isolate.isolate_key)
        row: Dict[str, Any] = {
            "isolate_key": isolate.isolate_key,
            "group_id": isolate.species_name,
            "kmer_cluster_id": isolate.kmer_cluster_id,
            "year": isolate.year,
            "species_name": isolate.species_name,
            "species_taxid": isolate.species_taxid,
            "host_name": isolate.host_name,
            "host_taxid": isolate.host_taxid,
            "host_order": isolate.host_order,
            "host_role": isolate.host_role,
            "country": isolate.country,
            "s_nt_length": isolate.s_nt_length,
            "m_nt_length": isolate.m_nt_length,
            "l_nt_length": isolate.l_nt_length,
            "n_segments_present": len(isolate.segments_present),
            "n_coding_complete": len(isolate.coding_complete_segments),
            "has_complete_triplet": int(isolate.has_complete_triplet),
            "has_nss_orf": int(isolate.has_nss_orf),
            "protein_length": len(isolate.concatenated_protein_aa),
            # Continuous divergence diagnostics only. The former binary reassortment
            # features were removed: they failed literature validation (see
            # compute_segment_divergence) and would have trained a model on sample size.
            "divergence_assessable": int(bool(profile and profile.assessable)),
            "sl_rank_separation": profile.sl_rank_separation if profile else 0.0,
            "m_rank_offset": profile.m_rank_offset if profile else 0.0,
            "taxonomic_incongruence": int(isolate.taxonomic_incongruence),
            "host_jump_status": host_jump_status(isolate.host_order, risk["reservoir_order"]),
            # Retained, not filtered -- see this function's docstring. Lets a caller
            # ablate/stratify on lab-vs-wild provenance instead of losing the isolates.
            "is_lab_derived": int(isolate.is_lab_derived),
        }
        # Confound controls, not biology -- see this function's docstring. Prefixed so a
        # caller can select/drop the whole block (df.filter(like="effort_")) in one shot.
        row.update(effort_by_species.get(isolate.species_name, _EFFORT_STATS_DEFAULT))
        row.update({
            # Labels
            "label_pathogenicity_class": risk["pathogenicity_class"],
            "label_disease_syndrome": risk["disease_syndrome"],
            "label_cfr_band": risk["cfr_band"],
            "label_is_human_derived": int(isolate.host_taxid == "9606"),
        })
        for column, kmer in enumerate(selected):
            row[f"kmer_{kmer}"] = float(reduced[position, column])
        rows.append(row)

    feature_columns = [f"kmer_{kmer}" for kmer in selected]
    n_label_bearing_groups = len({row["group_id"] for row in rows if row["label_cfr_band"] != -1})
    scaling = {
        "note": "Fit scaling inside your CV fold. These statistics are provided for "
                "reference and reproducibility only; they were computed on the full table "
                "and must not be applied before splitting.",
        "k": KMER_K,
        "alphabet": "amino acid",
        "selection": "top-variance k-mers (unsupervised)",
        "n_features": len(feature_columns),
        "features": feature_columns,
        "mean": {name: float(reduced[:, i].mean()) for i, name in enumerate(feature_columns)},
        "std": {name: float(reduced[:, i].std()) for i, name in enumerate(feature_columns)},
        "recommended_cv": "GroupKFold on group_id (species), or a temporal split on year",
        # Distinct species (group_id) with a real cfr_band, i.e. not the -1 'unknown'
        # placeholder -- how many groups actually carry a label at all, separate from how
        # many groups exist in total.
        "n_label_bearing_groups": n_label_bearing_groups,
        "caveat": (
            "Labels (label_pathogenicity_class, label_disease_syndrome, label_cfr_band) "
            "are constant within group_id: every isolate of one species carries the same "
            "label. Under GroupKFold on group_id, the effective sample size is therefore "
            f"the number of label-bearing SPECIES ({n_label_bearing_groups}), not the "
            "number of isolates in this table. Any split that is NOT grouped by group_id "
            "leaks completely, because 3-mer amino-acid composition separates species "
            "trivially -- a model can identify the species (and therefore its label) from "
            "the kmer_* columns alone, with no biology learned at all."
        ),
    }
    return rows, scaling


def build_ml_segment_rows(accessions: Sequence[AccessionRecord]) -> List[Dict[str, Any]]:
    rows = []
    for record in accessions:
        if record.segment not in {"S", "M", "L"}:
            continue
        risk = species_risk_profile(record.species_name)
        rows.append({
            "accession": record.accession,
            "sequence_hash": record.sequence_hash,
            "segment": record.segment,
            "segment_source": record.segment_source,
            "species_name": record.species_name,
            "species_taxid": record.species_taxid,
            "pathogenicity_class": risk["pathogenicity_class"],
            "cfr_band": risk["cfr_band"],
            "host_name": record.host_name,
            "host_taxid": record.host_taxid,
            "host_order": record.host_order,
            "host_role": record.host_role,
            "country": record.country,
            "year": record.year,
            "nt_length": record.nt_length,
            "completeness": record.completeness,
            "is_coding_complete": int(record.is_coding_complete),
            "is_lab_host": int(record.is_lab_host),
            "group_key": record.group_key,
        })
    return rows


def build_ml_sequence_rows(isolates: Sequence[IsolateRecord],
                           genomic_index: Dict[str, FastaRecord]) -> List[Dict[str, Any]]:
    """Wide, sequence-level ML table: one row per isolate with its raw S/M/L genomic
    sequences side by side, for tasks that want the sequences themselves rather than
    derived k-mer features (see build_ml_tables for that).

    Three filters, all isolate-level:
    - has_complete_triplet: coding-complete S, M and L all present (see
      build_isolate_records). Weaker "all three segments present" isolates are dropped,
      since a partial/non-coding segment here would just be a gap in the sequence column.
    - is_lab_derived excluded: unlike build_ml_tables (which keeps these, flagged, because
      k-mer composition doesn't care about passage history -- see its docstring), a raw
      sequence handed to a model as "wild-type" input is exactly the case that history
      would mislead.
    - is_vaccine_strain excluded: an attenuated strain's sequence is not representative of
      wild-type pathogenicity, same reasoning build_ml_tables already applies.

    An isolate can still be dropped here even after passing those filters if its genomic
    FASTA sequence is missing for one of the three accessions (coding-complete is
    determined from annotation length, not genomic-FASTA presence, so the two can
    disagree in rare cases) -- these are logged and skipped rather than emitting a row
    with a blank sequence column.
    """
    rows: List[Dict[str, Any]] = []
    skipped_missing_sequence = 0
    for isolate in isolates:
        if not isolate.has_complete_triplet or isolate.is_lab_derived or isolate.is_vaccine_strain:
            continue

        s_rec = genomic_index.get(isolate.s_accession)
        m_rec = genomic_index.get(isolate.m_accession)
        l_rec = genomic_index.get(isolate.l_accession)
        if not (s_rec and m_rec and l_rec):
            skipped_missing_sequence += 1
            continue

        risk = species_risk_profile(isolate.species_name)
        rows.append({
            "species_name": isolate.species_name,
            "isolate_name": isolate.isolate or isolate.strain or isolate.isolate_key,
            "s_sequence": s_rec.sequence,
            "m_sequence": m_rec.sequence,
            "l_sequence": l_rec.sequence,
            "cfr_reported_range": risk["cfr_reported_range"],
        })

    if skipped_missing_sequence:
        logging.warning(
            "ML sequence table: %s complete-triplet isolate(s) skipped, missing genomic "
            "FASTA sequence for at least one segment", skipped_missing_sequence)
    return rows


def build_suppressed_rows(accessions: Sequence[AccessionRecord],
                          isolates: Sequence[IsolateRecord]) -> List[Dict[str, Any]]:
    representative = {a for i in isolates for a in (i.s_accession, i.m_accession, i.l_accession) if a}
    hash_owner: Dict[str, str] = {}
    for record in accessions:
        if record.accession in representative and record.sequence_hash:
            hash_owner.setdefault(record.sequence_hash, record.accession)

    rows = []
    for record in accessions:
        if record.accession in representative:
            continue
        if record.sequence_hash and record.sequence_hash in hash_owner:
            reason = f"duplicate_sequence_of:{hash_owner[record.sequence_hash]}"
        elif not record.segment:
            reason = "segment_unresolved"
        elif record.warning_flags:
            reason = ";".join(record.warning_flags)
        else:
            reason = "not_best_representative_for_group"
        rows.append({
            "accession": record.accession,
            "sequence_hash": record.sequence_hash,
            "species_name": record.species_name,
            "species_taxid": record.species_taxid,
            "segment": record.segment,
            "completeness": record.completeness,
            "nt_length": record.nt_length,
            "biosample": record.biosample,
            "isolate": record.isolate,
            "strain": record.strain,
            "collection_date": record.collection_date,
            "location": record.location,
            "host_name": record.host_name,
            "host_taxid": record.host_taxid,
            "group_key": record.group_key,
            "suppression_reason": reason,
        })
    return rows


# --------------------------------------------------------------------------------------
# Input resolution and provenance
# --------------------------------------------------------------------------------------


def resolve_inputs(root: Path, output: Optional[Path] = None) -> InputPaths:
    root = root.expanduser().resolve()

    # Accept either the repo root (containing all_hanta/) or the package dir itself.
    candidates = [root / "all_hanta" / "ncbi_dataset" / "data", root / "ncbi_dataset" / "data", root]
    data_dir = next((c for c in candidates if (c / DATA_REPORT).exists()), None)
    if data_dir is None:
        raise FileNotFoundError(
            f"Could not locate {DATA_REPORT}. Looked in:\n  " + "\n  ".join(str(c) for c in candidates))

    package_dir = data_dir.parent.parent
    md5sum_file = package_dir / MD5SUM_FILE

    required = {
        "data_report": data_dir / DATA_REPORT,
        "annotation_report": data_dir / ANNOTATION_REPORT,
        "biosample_report": data_dir / BIOSAMPLE_REPORT,
        "dataset_catalog": data_dir / "dataset_catalog.json",
    }
    for name, path in required.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing required file ({name}): {path}")

    output_dir = (output.expanduser().resolve() if output else root / "output")
    output_dir.mkdir(parents=True, exist_ok=True)

    def optional(name: str) -> Optional[Path]:
        path = data_dir / name
        return path if path.exists() else None

    return InputPaths(
        root=root,
        package_dir=package_dir,
        data_dir=data_dir,
        data_report=required["data_report"],
        annotation_report=required["annotation_report"],
        biosample_report=required["biosample_report"],
        dataset_catalog=required["dataset_catalog"],
        md5sum_file=md5sum_file if md5sum_file.exists() else None,
        genomic_fasta=optional(GENOMIC_FASTA),
        protein_fasta=optional(PROTEIN_FASTA),
        cds_fasta=optional(CDS_FASTA),
        output_dir=output_dir,
        log_file=output_dir / "ortho_pipeline.log",
    )


def resolve_datasets_cli_version() -> str:
    """Version of NCBI's `datasets` CLI on PATH, if present. Purely informational -- this
    pipeline never shells out to `datasets` for anything else, and never requires it; the
    point is recording what tool a user could re-fetch an equivalent package with, for
    someone trying to reproduce this exact run later. Best-effort: any failure (not
    installed, no --version support, unexpected output) is recorded as the reason rather
    than raised, since provenance metadata must never be able to fail the run.
    """
    binary = shutil.which("datasets")
    if not binary:
        return "not found on PATH"
    try:
        result = subprocess.run([binary, "--version"], capture_output=True, text=True,
                                check=False, timeout=10)
        output = (result.stdout or result.stderr or "").strip()
        return output or "unknown (empty --version output)"
    except Exception as exc:  # noqa: BLE001 - provenance is best-effort, never fatal
        return f"unknown ({type(exc).__name__})"


def summarize(
    accessions: Sequence[AccessionRecord],
    isolates: Sequence[IsolateRecord],
    divergence: Sequence[SegmentDivergenceProfile],
    qc_rows: Sequence[Dict[str, Any]],
    checksums: Sequence[ChecksumResult],
    curation_issues: Sequence[Dict[str, str]],
    paths: InputPaths,
    cluster_threshold: float,
    min_species_size: int,
    data_records: Sequence[Dict[str, Any]] = (),
) -> Dict[str, Any]:
    catalog: Dict[str, Any] = {}
    if paths.dataset_catalog.exists():
        catalog = json.loads(paths.dataset_catalog.read_text(encoding="utf-8"))

    species_names = {a.species_name for a in accessions if a.species_name}
    status_counts = Counter(d.reassortment_status for d in divergence)

    # GenBank's own per-record releaseDate, carried straight through from the raw
    # data_report.jsonl records (not stored on AccessionRecord -- nothing downstream of
    # ingestion needs a per-record value, only this package-level min/max). This is what
    # identifies WHICH NCBI snapshot the package is, which the md5 integrity check below
    # does not: md5 proves these files are byte-identical to what NCBI published, not
    # which query/date produced them -- two packages can each pass their OWN md5sum.txt
    # while covering different accession sets, if NCBI was queried on different dates.
    release_dates = sorted(r["releaseDate"] for r in data_records if r.get("releaseDate"))

    return {
        "pipeline_version": VERSION,
        "generated_at_utc": datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "dataset_catalog_api_version": catalog.get("apiVersion", "unknown"),

        "provenance": {
            "taxon_id": ORTHOHANTAVIRUS_TAXID,
            "datasets_cli_version": resolve_datasets_cli_version(),
            "package_release_date_min": release_dates[0] if release_dates else None,
            "package_release_date_max": release_dates[-1] if release_dates else None,
            "note": "release_date_min/max identify WHICH NCBI snapshot this package is; "
                    "the md5 check under integrity below proves the files are intact, not "
                    "which query/date produced them. datasets_cli_version records what a "
                    "user could re-fetch an equivalent package with; it is never used by "
                    "this pipeline itself.",
        },

        "integrity": {
            "model": "verified against NCBI's published md5sum.txt (file level) and "
                     "nucleotide.sequenceHash / CRC32 (content level)",
            "md5sum_file_present": paths.md5sum_file is not None,
            "all_verified": all(r.verified for r in checksums) if checksums else False,
            "files": [asdict(r) for r in checksums],
            "sequence_hash_note": "NCBI strips leading zeros from sequenceHash; values are "
                                  "zero-padded to 8 hex characters before use.",
            "records_with_sequence_hash": sum(1 for a in accessions if a.sequence_hash),
            "distinct_sequence_hashes": len({a.sequence_hash for a in accessions if a.sequence_hash}),
        },

        "counts": {
            "accession_records": len(accessions),
            "failed_transforms": sum(1 for row in qc_rows if not row.get("transform_ok", True)),
            "canonical_species": len(species_names),
            "distinct_observed_organisms": len({a.observed_organism for a in accessions if a.observed_organism}),
            "distinct_observed_taxids": len({a.observed_taxid for a in accessions if a.observed_taxid}),
            "isolates": len(isolates),
            "complete_triplets": sum(1 for i in isolates if i.has_complete_triplet),
            "all_segments_present": sum(1 for i in isolates if len(i.segments_present) == 3),
            "records_with_host": sum(1 for a in accessions if a.host_name),
            "records_with_biosample": sum(1 for a in accessions if a.biosample),
            "human_derived_records": sum(1 for a in accessions if a.host_taxid == "9606"),
            "lab_host_records": sum(1 for a in accessions if a.is_lab_host),
            "vaccine_strain_records": sum(1 for a in accessions if a.is_vaccine_strain),
        },

        "segment_counts": dict(Counter(a.segment or "unresolved" for a in accessions)),
        "segment_source_counts": dict(Counter(a.segment_source for a in accessions)),
        "completeness_counts": dict(Counter(a.completeness for a in accessions)),
        "grouping_tier_counts": dict(Counter(a.grouping_tier for a in accessions)),
        "triplet_status_counts": dict(Counter(i.triplet_status for i in isolates)),
        "host_order_counts": dict(Counter(a.host_order or "unranked" for a in accessions)),
        "host_role_counts": dict(Counter(a.host_role for a in accessions)),
        "pathogenicity_counts": dict(Counter(
            species_risk_profile(a.species_name)["pathogenicity_class"] for a in accessions)),
        "host_jump_counts": dict(Counter(
            host_jump_status(a.host_order, species_risk_profile(a.species_name)["reservoir_order"])
            for a in accessions)),
        "reassortment_status_counts": dict(status_counts),
        "taxonomic_incongruence_count": sum(1 for i in isolates if i.taxonomic_incongruence),
        "qc_warning_counts": dict(Counter(
            w for row in qc_rows for w in str(row.get("warnings", "")).split(";") if w)),

        "methods": {
            "species_resolution": "canonical ICTV species taken from the node below genus "
                                  "Orthohantavirus (taxId 1980442) in virus.lineage, canonicalised "
                                  "by binomial name so legacy and current taxids for the same "
                                  "species merge.",
            "grouping": "isolate.name-first tier cascade; BioSample is a confirmation signal, not "
                        "the primary key, because it is present on <1% of records.",
            "sequence_distance": "cosine distance on 3-mer amino-acid composition (alignment-free).",
            "kmer_cluster_threshold": cluster_threshold,
            "reassortment": {
                "status": "NOT CALLED -- classifier withdrawn after failing validation",
                "what_is_emitted": "per-segment nearest neighbours, k-mer distances, and "
                                   "normalized rank separations (sl_rank_separation, "
                                   "m_rank_offset) as continuous diagnostics only",
                "min_species_size": min_species_size,
                "assessable_isolates": sum(1 for d in divergence if d.assessable),
                "why_withdrawn": "A classifier built on these rank separations produced the "
                                 "inverse of the reported literature ranking: ANDV 72% flagged "
                                 "(reported absent or sparse) vs PUUV 7.7% (reported 19.1-32%), "
                                 "and CHOV 100% at n=8. Sweeping the threshold from 0.25 to "
                                 "0.85 scaled all rates down without changing the ordering. "
                                 "Cosine distance on 3-mer amino-acid composition tracks "
                                 "within-species diversity and sample size, not tree topology.",
                "to_detect_reassortment_properly": "Build per-segment alignments (MAFFT) and "
                                                   "trees (IQ-TREE), then test topological "
                                                   "incongruence between the S, M and L trees.",
                "trustworthy_signal": "taxonomic_incongruence -- NCBI assigning an isolate's "
                                      "segments to different species. This is metadata, not "
                                      "inference, and it does feed the risk score.",
            },
            "reservoir_curation_audit": {
                "curated_species": len(SPECIES_RISK),
                "disagreements_with_data": len(curation_issues),
                "details": list(curation_issues),
                "note": "Curated reservoirs are checked against the modal observed host order "
                        "on every run. Six of the original entries disagreed with the package's "
                        "own host annotations and were corrected against the data.",
            },
            "LIMITATION": "Distances here are alignment-free proxies. The ICTV species criterion "
                          "for Hantaviridae is DEmARC with a pairwise evolutionary distance cutoff "
                          "of 0.1 under the WAG model, which requires multiple sequence alignment "
                          "and maximum-likelihood inference. Columns are named kmer_* and never "
                          "ped_* for this reason. Any published claim about a specific reassortant "
                          "needs confirmation with a real alignment and phylogeny (e.g. MAFFT + "
                          "IQ-TREE).",
            "citations": {
                "ictv_taxonomy": CIT_ICTV,
                "pathogenicity_and_cfr": CIT_ZOONOTIC,
                "ecdc_factsheet": CIT_ECDC,
                "reassortment": CIT_REASSORT,
            },
        },

        "inputs": {
            "package_dir": paths.package_dir.name,
            "files": sorted(p.name for p in paths.data_dir.iterdir() if p.is_file()),
        },
    }


# --------------------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------------------


def run_pipeline(paths: InputPaths, args: argparse.Namespace) -> None:
    logging.info("Verifying package integrity against NCBI md5sum.txt")
    checksums = verify_inputs(paths)
    report_checksums(checksums, allow_mismatch=args.allow_checksum_mismatch)

    logging.info("Loading reports")
    data_records = list(iter_jsonl(paths.data_report))
    logging.info("  data_report: %s records", len(data_records))
    annotation_index = build_annotation_index(iter_jsonl(paths.annotation_report))
    logging.info("  annotation_report: %s records (%s index keys)",
                 indexed_record_count(annotation_index), len(annotation_index))
    biosample_index = build_biosample_index(iter_jsonl(paths.biosample_report))
    logging.info("  biosample_report: %s samples", len(biosample_index))

    logging.info("Parsing FASTA")
    genomic_index = parse_fasta(paths.genomic_fasta)
    protein_index = parse_fasta(paths.protein_fasta)
    logging.info("  genomic: %s sequences, protein: %s sequences",
                 indexed_record_count(genomic_index), indexed_record_count(protein_index))

    # Serial by design: this work is pure-Python and GIL-bound, so the previous
    # ThreadPoolExecutor added 15,991 futures of overhead for no parallelism.
    logging.info("Transforming accession records")
    accessions: List[AccessionRecord] = []
    qc_rows: List[Dict[str, Any]] = []
    for record in data_records:
        transformed, qc_row = transform_one_record(
            record, annotation_index, biosample_index, genomic_index, protein_index)
        if transformed is not None:
            accessions.append(transformed)
        qc_rows.append(qc_row)
    logging.info("  %s transformed, %s failed", len(accessions), len(data_records) - len(accessions))

    curation_issues = audit_reservoir_curation(accessions)
    for issue in curation_issues:
        logging.warning("Reservoir curation disagrees with data for %s: curated %s, observed %s",
                        issue["species"], issue["curated_reservoir"], issue["modal_observed_host"])
    if not curation_issues:
        logging.info("Reservoir curation audit: all curated reservoirs agree with observed hosts")

    logging.info("Reconstructing isolates")
    isolates = build_isolate_records(accessions)
    logging.info("  %s isolates, %s complete triplets",
                 len(isolates), sum(1 for i in isolates if i.has_complete_triplet))

    logging.info("Clustering by k-mer distance (threshold %.3f)", args.cluster_threshold)
    assign_kmer_clusters(isolates, args.cluster_threshold)

    logging.info("Computing per-segment divergence diagnostics")
    divergence = compute_segment_divergence(isolates, min_species_size=args.min_species_size)
    divergence_by_key = {d.isolate_key: d for d in divergence}
    assessable = sum(1 for d in divergence if d.assessable)
    incongruent = sum(1 for i in isolates if i.taxonomic_incongruence)
    logging.info("  %s assessable, %s with taxonomic incongruence", assessable, incongruent)
    logging.warning(
        "Reassortment is NOT called by this pipeline. The k-mer classifier that used to do "
        "so failed validation against the literature (it ranked ANDV above PUUV, the exact "
        "inverse of reported frequencies) because the metric tracks within-species "
        "diversity and sample size, not tree topology. Assessable isolates are marked "
        "'requires_phylogenetic_confirmation'; confirm candidates with MAFFT + IQ-TREE.")

    logging.info("Building ML tables")
    ml_rows, scaling = build_ml_tables(accessions, isolates, divergence_by_key, args.max_features)

    out = paths.output_dir

    logging.info("Writing outputs")
    output_digests: Dict[str, str] = {}

    def write_csv(name: str, rows: Sequence[Dict[str, Any]]) -> None:
        output_digests[name] = safe_write_csv(out / name, rows)

    def write_json(name: str, payload: Any) -> None:
        output_digests[name] = safe_write_json(out / name, payload)

    write_csv("ortho_accession_table.csv", build_accession_rows(accessions))
    write_csv("ortho_species_table.csv", build_species_rows(accessions, isolates))
    write_csv("ortho_isolate_table.csv", build_isolate_rows(isolates, divergence_by_key))
    write_csv("ortho_host_table.csv", build_host_rows(accessions, isolates))
    write_csv("ortho_ml_feature_matrix.csv", ml_rows)
    write_csv("ortho_ml_sequence_table.csv", build_ml_sequence_rows(isolates, genomic_index))
    write_csv("ortho_suppressed.csv", build_suppressed_rows(accessions, isolates))
    write_csv("ortho_qc_report.csv", qc_rows)
    write_json("ortho_feature_scaling.json", scaling)

    summary = summarize(accessions, isolates, divergence, qc_rows, checksums,
                        curation_issues, paths, args.cluster_threshold,
                        args.min_species_size, data_records)
    write_json("ortho_summary.json", summary)

    # ortho_risk_table.csv is always phylo-blind now -- Stage 2 (phylogeny.py) is a separate
    # script/run; combine its ortho_reassortment_table.csv on isolate_key if you need both
    # (see README's Outputs table). This is byte-identical to the pre-split behavior when
    # --phylogenetics wasn't passed, which was already the common case.
    write_csv("ortho_risk_table.csv", build_risk_rows(isolates, divergence_by_key))

    # Written once, here, after every output that might exist -- including the optional
    # reassortment table and the risk table above -- is already on disk, so it always
    # covers exactly what's present in `out`. Run `sha256sum -c ortho_output_checksums.txt`
    # from the output directory to confirm nothing was truncated or altered afterward (by
    # an interrupted write, a bad copy, git-lfs mishandling, etc).
    write_output_checksums(out / "ortho_output_checksums.txt", output_digests)

    logging.info("Done. %s accessions -> %s species -> %s isolates (%s complete triplets)",
                 summary["counts"]["accession_records"],
                 summary["counts"]["canonical_species"],
                 summary["counts"]["isolates"],
                 summary["counts"]["complete_triplets"])


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Orthohantavirus NCBI Datasets production pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", default=".",
                        help="Package root containing all_hanta/ncbi_dataset/data")
    parser.add_argument("--output", help="Output directory; defaults to <input>/output")
    parser.add_argument("--allow-checksum-mismatch", action="store_true",
                        help="Continue even if inputs do not match NCBI's md5sum.txt "
                             "(the mismatch is still recorded in ortho_summary.json)")
    parser.add_argument("--cluster-threshold", type=float, default=0.02,
                        help="Cosine k-mer distance for single-linkage clustering within a species")
    parser.add_argument("--max-features", type=int, default=200,
                        help="Number of top-variance k-mer features in the ML matrix")
    parser.add_argument("--min-species-size", type=int, default=8,
                        help="Minimum assessable isolates in a species before reassortment "
                             "is called at all")
    # --concordance-tolerance / --discordance-threshold were removed along with the
    # withdrawn reassortment classifier; no threshold made its ranking agree with the
    # literature, so exposing one invited false confidence.
    parser.add_argument("--verbose", action="store_true", help="Enable debug logging")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    paths = resolve_inputs(Path(args.input), output=Path(args.output) if args.output else None)
    configure_logging(paths.log_file, verbose=args.verbose)
    logging.info("Orthohantavirus pipeline v%s", VERSION)
    logging.info("Package: %s", paths.data_dir)
    logging.info("Output : %s", paths.output_dir)
    run_pipeline(paths, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


