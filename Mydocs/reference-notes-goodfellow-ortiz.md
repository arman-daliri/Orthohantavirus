# Reference notes: Goodfellow 2024 and Ortiz 2025

Extracted before the source PDFs were removed from the repo. **Both papers are Open
Access (CC BY)** — the PDFs were removed to keep the repo small, not for licensing
reasons, and either can be re-downloaded from its DOI at any time.

---

## Goodfellow et al. (2024) — "Limestone Canyon virus" is Montaño virus

Goodfellow SM, Nofchissey RA, Morley VJ, Coan KE, Schwalm KC, Cook JA, Dunnum JL,
Hanfelt-Goade D, Dinwiddie DL, Domman DB, Dragoo JW, Kuhn JH, Bradfute SB.
*Genome sequencing identifies "Limestone Canyon virus" as Montaño virus
(Hantaviridae: Orthohantavirus montanoense) circulating in brush deermice in New Mexico.*
**npj Viruses 2:11 (2024).** doi:[10.1038/s44298-024-00016-6](https://doi.org/10.1038/s44298-024-00016-6)
Licence: **CC BY 4.0** (Open Access).

### Core finding
"Limestone Canyon virus" (LSCV) is **not a distinct virus**. Its coding-complete genome
groups with **Montaño virus (MTNV, *Orthohantavirus montanoense*)**. Because the MTNV
coding-complete genome was deposited in 2012 — before LSCV was characterised — MTNV has
priority, and **"Limestone Canyon virus"/"LSCV" are junior synonyms that should be
abandoned**.

### Specimen and accessions
- Host: brush deermouse (*Peromyscus boylii*), MSB:Mamm:332771, Red Canyon, Manzano
  Mountains, central New Mexico, September 2019.
- Coding-complete S/M/L deposited as **OR148902–OR148904**.
- Eight 2006 specimens from Hidalgo and Grant Counties (SW New Mexico):
  **OR552601–OR552625**.
- Museum metadata available via Arctos (https://arctosdb.org).

### Divergence of "LSCV" from MTNV
| Segment | Nucleotide divergence | Amino-acid divergence |
|---|---|---|
| S | 29.9 % | 6.5 % |
| M | 22.1 % | 11.8 % |
| L | 6.5 % | 9.6 % |

### Why this matters to this project
1. **DEmARC.** Species are established by the ICTV *Hantaviridae* Study Group only when a
   coding-complete exemplar (S + M + L) exists and differs sufficiently under **DEmARC**
   (DivErsity pArtitioning by hieRarchical Clustering), per Laenen et al. 2019. This is
   the same demarcation framework the RPI design spec §7 refers to when it calls for a
   "DEmARC-style analysis" — note the spec means the *approach*, not the tool.
2. **Reassortment.** The paper explicitly motivates reassortment surveillance: multiple
   *Peromyscus* species are syntopic in New Mexico with overlapping distributions, giving
   "high potential for multiple orthohantaviruses to co-occur in the same area and even
   in the same individual rodent," and calls for assessing "intra- and inter-species
   genomic segment reassortment and hence the possibility for novel orthohantavirus
   emergence."
3. **Naming hygiene.** *Orthohantavirus montanoense* (MTNV) is present in our dataset
   (1 complete triplet). Any output naming "Limestone Canyon virus" should be treated as
   a synonym of MTNV.

### Read composition (BLAST hits, MSB:Mamm:332771 WGS)
MTNV 24.8 %, KHAV 17.6 %, SNV 10.4 %, CADV 3.3 %, FUSV 3.3 %, RKPV 2.4 %, HTNV 2.4 %.

### Methods worth mirroring
CLUSTAL Omega alignment; MEGA-X maximum likelihood, **General Time Reversible**,
**1,000 bootstraps**; trees rendered with iTOL. (Our pipeline uses IQ-TREE, but the
GTR + 1000-bootstrap convention matches.)

---

## Ortiz et al. (2025) — rodent hosts and genotypes in Southern South America

Ortiz N, Pinotti JD, Andreo V, González-Ittig RE, Gardenal CN.
*Orthohantavirus rodent hosts and genotypes in Southern South America: A narrative review.*
**PLOS Negl Trop Dis 19(9):e0013489 (2025).**
doi:[10.1371/journal.pntd.0013489](https://doi.org/10.1371/journal.pntd.0013489)
Licence: **CC BY** (Open Access).

### Scope
Narrative review covering Argentina, Chile, Paraguay and Uruguay: **14 rodent reservoir
species hosting 15 distinct orthohantavirus genotypes**. HPS case-fatality in South
America **12–40 %**; HFRS **5–10 %**.

### Table 1 — reservoir host → genotype (condensed)

| Reservoir host | Orthohantavirus / genotype | Distribution | Human disease |
|---|---|---|---|
| *Akodon azarae* | *O. andesense* (Pergamino) | Uruguay, E Argentina, Paraguay, S Brazil | Unknown |
| *Akodon montensis* | Jaborá (Apé Aimé synonymised under Jaborá, Rivera 2015) | NE Argentina, Paraguay, Brazil | Unknown |
| *Calomys callosus* | *O. negraense* | Argentina, Paraguay, E Bolivia, C-W Brazil | HPS |
| *Calomys fecundus* | *O. negraense* | NW Argentina, S Bolivia | HPS (Argentina) |
| *Calomys laucha* | *O. negraense* | N Paraguay, SE Bolivia, C Argentina, Brazil | HPS (Paraguay) |
| *Holochilus chacarius* | Alto Paraguay | C-E/N-E Argentina, W Paraguay, Brazil | HPS |
| *Necromys lasiurus* | *O. andesense* (Maciel, Araraquara) | N Argentina, Bolivia, Paraguay, Brazil | Araraquara → HPS in Brazil; Maciel not linked to human cases |
| *Oligoryzomys chacoensis* | *O. andesense* (Orán) | SE Bolivia, W Paraguay, SW Brazil, N Argentina | HPS |
| *Oligoryzomys flavescens* s.s. | HU39694, *O. andesense* (Lechiguanas, Buenos Aires, Central Plata) | C/NE Argentina, Uruguay | HPS |
| *Oligoryzomys longicaudatus* | *O. andesense* (Andes) | S Argentina, C-S Chile | HPS; **person-to-person transmission** |
| *Oligoryzomys nigripes* | Juquitiba, *O. andesense* (Lechiguanas) | C/SE Brazil, E Paraguay, Uruguay, NE Argentina | HPS |
| *Oligoryzomys occidentalis* | *O. andesense* (Bermejo) | C-W/NW Argentina, S/C-W Bolivia | HPS |
| *Scapteromys aquaticus* | Leyes | C-E/N Argentina, S Paraguay | Unknown |
| *Rattus norvegicus* | Seoul virus (SEOV) | Worldwide | HFRS (Asia, Europe) |

### Why this matters to this project
1. **Factor E host mapping.** This is a curated reservoir-host → genotype table for a
   region our dataset covers heavily (ANDV has 45 complete triplets). It is a source for
   validating `host_rank_proximity` and the ecological-opportunity term against published
   host assignments rather than inferring them from record metadata alone.
2. **The intra-species blind spot, again.** Note how many distinct *named genotypes*
   (Pergamino, Maciel, Araraquara, Orán, Lechiguanas, Buenos Aires, Central Plata, Andes,
   Bermejo) sit **inside** *O. andesense*. Reassortment among these is exactly the
   intra-species class gate 1 cannot see, since all of them assign to one reference
   species. This is direct literature support for prioritising gate 4.
3. **Sympatry drives co-infection.** The review notes spillover is more common for
   cricetid-borne orthohantaviruses than other host groups, and that sympatric,
   ecologically similar or phylogenetically related hosts facilitate horizontal
   transmission — the ecological precondition Factor E is meant to score.
