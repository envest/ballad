#!/usr/bin/bash

dir_cache="${HOME}/.vep"
assembly="GRCh38"
fasta="${dir_cache}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

vep --input_file test.VEP_input.hg38.vcf --format vcf \
	--output_file output.txt --vcf \
	--force_overwrite \
	--no_stats \
	--no_check_variants_order \
	--cache \
	--offline \
	--dir_cache ${dir_cache} \
	--fasta ${fasta} \
	--assembly ${assembly} \
	--variant_class \
	--nearest symbol \
	--gene_phenotype \
	--regulatory \
	--show_ref_allele \
	--hgvs \
	--symbol \
	--canonical \
	--domains \
	--check_existing \
	--exclude_null_alleles \
	--af \
	--max_af \
	--af_gnomadg
