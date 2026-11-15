#!/bin/bash

mature_file=mature.fa
if [ $1 ] && [ -f "$2/mirnas.fa" ]; then
    echo "Loaded $mature_file from cache."
else 
    echo "Downloading $mature_file.gz from miRBase..."
    curl "https://www.mirbase.org/ftp/CURRENT/$mature_file.gz" | gunzip -c > "$2/mirnas.fa"
fi

gtf_file=Homo_sapiens.GRCh38.$3.chr.gtf
if [ $1 ] && [ -f "$2/annotation.gtf" ]; then
    echo "Loaded $gtf_file from cache."
else 
    echo "Downloading $gtf_file.gz from Ensembl..."
    curl "ftp://ftp.ensembl.org/pub/release-$3/gtf/homo_sapiens/$gtf_file.gz" | gunzip -c > "$2/annotation.gtf"
fi

mane_file=MANE.GRCh38.v$4.select_ensembl_genomic.gtf
if [ $1 ] && [ -f "$2/mane.gtf" ]; then
    echo "Loaded $mane_file from cache."
else 
    echo "Downloading $mane_file.gz from NCBI..."
    curl "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_$4/$mane_file.gz" | gunzip -c > "$2/mane.gtf"
fi