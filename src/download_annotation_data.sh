#!/bin/bash

mature_file=mature.fa
if [ $1 ] && [ -f "$2/mirnas.fa" ]; then
    echo "Loaded $mature_file from cache."
else 
    echo "Downloading $mature_file from miRBase..."
    # mirbase no longer supplies a fa.gz and the fa is littered with html tags, remove them if they're there
    curl -fL "https://www.mirbase.org/download/CURRENT/$mature_file" | sed -e 's/<br>/\n/g' -e 's/&gt;/>/g' -e 's/<\/\?p>//g' > "$2/mirnas.fa"
fi

gtf_file=Homo_sapiens.GRCh38.$3.chr.gtf
if [ $1 ] && [ -f "$2/annotation.gtf" ]; then
    echo "Loaded $gtf_file from cache."
else 
    echo "Downloading $gtf_file.gz from Ensembl..."
    curl "ftp://ftp.ensembl.org/pub/release-$3/gtf/homo_sapiens/$gtf_file.gz" | gunzip -c > "$2/annotation.gtf"
fi

mane_file=MANE.GRCh38.v$4.ensembl_genomic.gtf
if [ $1 ] && [ -f "$2/mane.gtf" ]; then
    echo "Loaded $mane_file from cache."
else 
    echo "Downloading $mane_file.gz from NCBI..."
    echo "$mane_file"
    curl "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_$4/$mane_file.gz" | gunzip -c > "$2/mane.gtf"
fi