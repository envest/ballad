# vep2loom: convert a vartrix mtx file to loom
# inputs: vep.vcf output.loom

import sys
import numpy as np
import loompy

vcf_filename = sys.argv[1]
output_loom_filename = sys.argv[2]

vcf_file = open(vcf_filename, "r")

annotation_keys = ["Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", 
  "Feature_type", "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", 
  "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", 
  "Existing_variation", "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", 
  "HGNC_ID", "CANONICAL", "HGVS_OFFSET", "AF", "CLIN_SIG", "SOMATIC", "PHENO"]
n_annotation_keys = len(annotation_keys)

# convert vcf_file to dict of row_attr and layers
chr_list = []
pos_list = []
ref_list = []
alt_list = []
variant_list = []
annotations_dict = {}
n_records = 0
for k in annotation_keys:
  annotations_dict[k] = []
  
for line in vcf_file:
    if line.startswith("#"):
        continue
    else:
        n_records += 1
        record = line.strip().split("\t")
        CHROM = record[0]
        POS = int(record[1])
        REF = record[3]
        ALT = record[4]
        presplit_annotations = record[7]
        if presplit_annotations.startswith("CSQ="):
            annotations = record[7].split("=")[1].split("|")
            annotations_converted = [None if x == '' else x for x in annotations]
        else:
            annotations_converted = [None]*n_annotation_keys
        chr_list.append(CHROM)
        pos_list.append(POS)
        ref_list.append(REF)
        alt_list.append(ALT)
        variant_list.append(":".join([CHROM, str(POS), REF, ALT]))
        for ki in range(n_annotation_keys):
          annotations_dict[annotation_keys[ki]].append(annotations_converted[ki])

vcf_file.close()

# set up loom components
row_attrs = {"CHROM" : np.array(chr_list), "POS" : np.array(pos_list), "REF" : np.array(ref_list), "ALT" : np.array(alt_list), "variant" : np.array(variant_list)} | {k : np.array(v) for k,v in annotations_dict.items()}

col_attrs = {"VEP" : [0]}

layers = {"" : np.zeros((n_records, 1))}

# create loom output
loompy.create(output_loom_filename,
  layers,
  row_attrs,
  col_attrs)

