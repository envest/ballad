#!/usr/bin/bash

input_vcf=$1
output_vcf=$2
assembly=$3

dir_cache="${HOME}/.vep"

if [ ${assembly} == "GRCh38" ]; then

	fasta="${dir_cache}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

elif [ ${assembly} == "GRCh37" ]; then

	fasta="${dir_cache}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"

else

	echo "Assembly must be GRCh37 or GRCh38"
	exit 1

fi

vep --input_file ${input_vcf} --format vcf \
	--output_file ${output_vcf} --vcf \
	--force_overwrite \
	--no_stats \
	--no_check_variants_order \
	--cache \
	--offline \
	--dir_cache ${dir_cache} \
	--fasta ${fasta} \
	--assembly ${assembly} \
	--pick \
	--hgvs \
	--symbol \
	--canonical \
	--check_existing \
	--exclude_null_alleles \
	--af

