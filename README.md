# Orthohantavirus
*************************************************************************
*************************************************************************
# Step1 Bash: download all Orthohantavirus virus genomes + BioSample + CDS + proteins
datasets download virus genome taxon 1980442 \
  --include genome,biosample,cds,protein \
  --filename ortho.zip

# Bash: unzip into a working directory
unzip ortho.zip -d ortho

dataformat tsv virus-genome \
  --package ortho.zip \
  --fields accession,segment,completeness,length,biosample-acc,'virus-infraspecific-isolate','virus-infraspecific-strain','isolate-collection-date',geo-location,host-name \
  > metadata.tsv

**************************************************************************
**************************************************************************

# Step2: run all_hanta46.py code to classify the groups
## Code Report: `all_hanta46.py`

## What This Script Does
This pipeline processes raw NCBI orthohantavirus packages into ML-ready datasets. It cross-checks metadata across multiple files, matches nucleotide and protein records, applies segment and completeness rules, removes likely lab-passaged records, and flags possible reassortment. 

## Input Files and Checks
The script requires six files in a hard-coded local folder: `data_report.jsonl`, `annotation_report.jsonl`, `biosample_report.jsonl`, `genomic.fna`, `cds.fna`, and `protein.faa`. The function `ensure_inputs()` checks that all six are present and stops immediately if any are missing. This makes the pipeline strict but not portable, because the file paths are fixed in the code.

## Core Methods
The script first loads the three JSONL files line by line with `iter_jsonl()` and raises an error if a line is not valid JSON. It then builds lookup tables with `build_annotation_index()` and `build_biosample_index()`, and parses the FASTA files with `parse_fasta()` so sequences can be found by accession or base accession. For each record, it flattens nested JSON, finds the accession, links matching annotation and BioSample data, and merges everything into one context for later searches. 

The script extracts organism name, taxonomy ID, isolate, strain, collection date, country, region, isolation source, tissue source, BioSample accession, and BioProject accession from these contexts. Host candidates are found from host fields first, then from isolation-source fields if needed, and host labels are cleaned with `normalize_host_label()`. Protein records are matched by accession, and if no direct match is found, isolate or strain text is used to search protein FASTA headers. 

## Segment and Completeness
`choose_segment()` assigns S, M, or L by searching metadata text and protein headers for key words. Nucleocapsid or nucleoprotein points to S, glycoprotein-related words point to M, and polymerase-related words point to L. If no text clue is found, it falls back to length: 5000 bp or more suggests L, 2500 bp or more suggests M, and 1000 bp or more suggests S. 

`determine_completeness()` gives one of six labels. `complete` is used when the record says complete or the text says “complete genome” or “complete sequence”. `partial` is used when the structured field says partial, and `coding_complete_only` is used when the text says “complete CDS” without partial. `likely_complete` is assigned when the nucleotide length meets the minimum expected length for the detected segment: 5000 for L, 2500 for M, and 1000 for S. In simple English, this means the sequence is probably complete because it is long enough for that segment type. `partial_or_unknown` is used when there is a length but no other rule fits, and `unknown` is used when there is no length. For grouping, `complete`, `coding_complete_only`, and `likely_complete` count as passing statuses. 

## Grouping and Filtering
The script removes records that look like lab-passaged or lab-adapted sequences by searching for words such as `vero`, `passage`, `cell line`, `adapted`, and `in vitro` in the combined metadata text. It then builds a group key. If a BioSample accession exists, it uses that first; otherwise it uses isolate or strain plus collection date, location, and host; if nothing useful exists, it falls back to the accession itself. The confidence score for the compound key is the fraction of the four fields that are known, rounded to three decimals. 

After grouping, only groups with complete S, M, and L segments are kept. For each kept group, the script picks the best segment representative by preferring complete records, then RefSeq accessions, then longer sequences. 

It predicts a host for each segment using the most common host among complete members of that segment.

Labels reassortment as `0`, `1`, or `"unknown"` depending on whether the three segment host predictions match.

## Outputs
The script writes three output files. `orthohantavirus_ml_ready.json` stores one enriched record per retained accession and includes accession, group ID, segment, completeness, host data, protein accessions, sequence length, group members, and predicted hosts for S, M, and L. `orthohantavirus_host_prediction_master.csv` stores one row per retained group with the group ID, primary host label, representative protein accessions, reassortment status, group size, observed hosts, country, and year. `orthohantavirus_summary.json` stores counts for records, groups, segment types, completeness types, reassortment, parse errors, lab-dropped records, and output file paths. 

## Observations
The report is strongest when it highlights the layered fallback logic and the separate group-level versus record-level pattern counts. The hard-coded paths and keyword-based lab filter are the main weaknesses. One important correction is that `complete_members` in the CSV is hard-coded to `3`, while the actual number of records in a group is tracked separately in `group_size`. 
