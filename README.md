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

# Step2: run code to classify the groups
