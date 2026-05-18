# Orthohantavirus

# Step1: Download data_report.tsv from NCBI

# once conda works
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli

# 1. Download (metadata‑only package for Orthohantavirus)
datasets download virus genome taxon Orthohantavirus \
  --include none --filename orthohanta.zip

# 2. Unpack the ZIP
unzip orthohanta.zip -d orthohanta

# 3. Create TSV from the internal JSONL
dataformat tsv virus-genome \
  --inputfile orthohanta/ncbi_dataset/data/data_report.jsonl \
  --fields accession,length,completeness,biosample-acc,virus-infraspecific-isolate,virus-infraspecific-strain,collection-date,geo-location,host-name,virus-name,virus-sci-name \
  > data_report.tsv


dataformat tsv virus-genome   --inputfile orthohanta/ncbi_dataset/data/data_report.jsonl   --fields accession,bioprojects,biosample-acc,completeness,gene-count,length,segment,release-date,update-date,sourcedb,geo-location,geo-region,geo-state,host-common-name,host-name,host-tax-id,host-infraspecific-isolate,host-infraspecific-strain,host-infraspecific-breed,host-infraspecific-cultivar,host-infraspecific-ecotype,host-infraspecific-sex,virus-common-name,virus-name,virus-tax-id,virus-infraspecific-isolate,virus-infraspecific-strain,virus-infraspecific-breed,virus-infraspecific-cultivar,virus-infraspecific-ecotype,virus-infraspecific-sex,isolate-collection-date,isolate-lineage,isolate-lineage-source,purpose-of-sampling,submitter-country,submitter-affiliation,submitter-names,sra-accs,protein-count,matpeptide-count,mol-type,is-annotated,is-lab-host,is-vaccine-strain,lab-host,virus-pangolin,host-pangolin   > data_report.tsv

head -n 5 data_report.tsv

# Step2: run hanta.ipynb to get orthohanta_cleaned.tsv
