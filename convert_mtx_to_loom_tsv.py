# convert_mtx_to_loom_tsv: convert a vartrix mtx file to loom and tsv
# inputs: consensus.mtx variants.vcf barcodes.tsv output.loom output.csv

import sys
import numpy as np
import loompy

mtx_filename = sys.argv[1]
var_filename = sys.argv[2]
bc_filename = sys.argv[3]
output_loom_filename = sys.argv[4]
output_tsv_filename = sys.argv[5]

mtx_file = open(mtx_filename, "r")
var_file = open(var_filename, "r")
bc_file = open(bc_filename, "r")

# get from n_var and n_bc from third row of mtx file
mtx_file.readline()
mtx_file.readline()
n_var, n_bc, n_nonzero = [int(x) for x in mtx_file.readline().strip().split()]

if n_var == 0:
    print("No variants found in " + mtx_filename)
    sys.exit()

# set up consensus array with all missing values
consensus_array = np.zeros((n_var, n_bc)) + 3

# convert var_file to dict of row_attr
# match names to tapestri loom
chr_list = []
pos_list = []
ref_list = []
alt_list = []
variant_list = []
for line in var_file:
    if line.startswith("#"):
        continue
    else:
        record = line.strip().split("\t")
        CHROM = record[0]
        POS = int(record[1])
        REF = record[3]
        ALT = record[4]
        chr_list.append(CHROM)
        pos_list.append(POS)
        ref_list.append(REF)
        alt_list.append(ALT)
        variant_list.append(":".join([CHROM, str(POS), REF, ALT]))

row_attrs = {"CHROM" : np.array(chr_list), 
        "POS" : np.array(pos_list),
        "REF" : np.array(ref_list),
        "ALT" : np.array(alt_list),
        "variant" : np.array(variant_list)}

col_attrs = {"barcode": [bc.strip() for bc in bc_file.readlines()]}

mtx_file.close()
var_file.close()
bc_file.close()

# use index values from mtx to update loom
mtx_file.close()
mtx_file = open(mtx_filename, "r")
mtx_file.readline()
mtx_file.readline()
mtx_file.readline()

for line in mtx_file:
    var, bc, gt = [int(x) for x in line.strip().split()]
    if gt == 1: # REF only
        consensus_array[var - 1, bc - 1] = 0
    elif gt in [2,3]: # ALT only or ALT+REF
        consensus_array[var - 1, bc - 1] = 1
    else:
        exit("sparse consensus.mtx contains value other than 1,2,3")

# create loom output
loompy.create(output_loom_filename, 
        consensus_array, 
        {k:np.array(v) for k,v in row_attrs.items()}, 
        {k:np.array(v) for k,v in col_attrs.items()})

# create tsv output
output_tsv = open(output_tsv_filename, "w")

header_line = "\t".join(["variant"] + list(col_attrs["barcode"])) + "\n"
output_tsv.write(header_line)

for row_index in range(n_var):
    var = ":".join([chr_list[row_index], str(pos_list[row_index]), str(ref_list[row_index]), str(alt_list[row_index])])
    write_line = "\t".join([var] + [str(int(x)) for x in consensus_array[row_index]]) + "\n"
    output_tsv.write(write_line)

output_tsv.close()

