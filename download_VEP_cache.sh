#!/usr/bin/bash

set -euo pipefail

cache_dir="${HOME}/.vep"

mkdir -p ${cache_dir}

cd ${cache_dir}

curl -O https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh38.tar.gz
tar xzf homo_sapiens_vep_111_GRCh38.tar.gz
rm -f homo_sapiens_vep_111_GRCh38.tar.gz

curl -O https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh37.tar.gz
tar xzf homo_sapiens_vep_111_GRCh37.tar.gz
rm -f homo_sapiens_vep_111_GRCh37.tar.gz

curl -O https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
bgzip -i Homo_sapiens.GRCh38.dna.primary_assembly.fa

curl -O https://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
bgzip -i Homo_sapiens.GRCh37.75.dna.primary_assembly.fa

