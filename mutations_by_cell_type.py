# Analyze mutations by cell type
# python mutations_by_cell_type.py sample.loom mutations.tsv cell_barcodes.tsv output_dir output_prefix

import os
import sys
import loompy
import numpy as np
from scipy.stats import binomtest
from scipy.stats import false_discovery_control

loom_filename = sys.argv[1]
mutations_filename = sys.argv[2]
cell_barcodes_filename = sys.argv[3]
output_dir = sys.argv[4]
output_prefix = sys.argv[5]

output_filename = os.path.join(output_dir, output_prefix + ".mutations_by_cell_type.tsv")

if not os.path.isdir(output_dir):
  sys.exit(output_dir + " directory does not exist")

cell_type_dict = {}
with open(cell_barcodes_filename, "r") as f:
    for line in f.readlines():
        bc, tumor_normal, level_one, level_two, umap1, umap2 = line.strip().split()

        compartment_0 = ":".join([str(x) for x in ["L0", tumor_normal, "NA", "NA"]])
        compartment_1 = ":".join([str(x) for x in ["L1", tumor_normal, level_one, "NA"]])
        compartment_2 = ":".join([str(x) for x in ["L2", tumor_normal, level_one, level_two]])
  
        for compartment in [compartment_0, compartment_1, compartment_2]:
            if compartment in cell_type_dict:
                cell_type_dict[compartment].append(bc)
            else:
                cell_type_dict[compartment] = [bc]

ds = loompy.connect(loom_filename)

mutation_dict = {}
mutation_prop = {}
loom_mutation_list = []

for i in range(ds.shape[0]):
  loom_mutation_list.append(":".join([ds.ra["CHROM"][i], str(ds.ra["POS"][i]), ds.ra["REF"][i], ds.ra["ALT"][i]]))

with open(mutations_filename, "r") as f:
    for line in f.readlines():
        m = line.strip()
        if m in loom_mutation_list:
            m_index = np.where([x == m for x in loom_mutation_list])[0][0]
            mutation_dict[m] = m_index
            n_00 = sum(ds[m_index,] == 0)
            n_01 = sum(ds[m_index,] == 1)
            n_11 = sum(ds[m_index,] == 2)
            mutation_prop[m] = (n_01 + n_11)/(n_00 + n_01 + n_11)
        else:
            print("Mutation " + m + " not in loom file")
  
output_header_list = ["mutation", "annotation_level", "tumor_normal", "level_one", "level_two", "cell_type_total", "cell_type_with_coverage", "n_00", "n_01", "n_11", "n_na", "cell_type_total_mutated", "cell_type_proportion_mutated", "overall_proportion_mutated", "binomial_pvalue"]
output_list_dict = {k : [] for k in output_header_list}

for m,i in mutation_dict.items():

    print(m)

    cells_00 = ds.ca["barcode"][np.where(ds[i,] == 0)]
    cells_01 = ds.ca["barcode"][np.where(ds[i,] == 1)]
    cells_11 = ds.ca["barcode"][np.where(ds[i,] == 2)]
    cells_na = ds.ca["barcode"][np.where(ds[i,] == 3)]
  
    for compartment in cell_type_dict.keys():

        n_00 = sum([x in cells_00 for x in cell_type_dict[compartment]])
        n_01 = sum([x in cells_01 for x in cell_type_dict[compartment]])
        n_11 = sum([x in cells_11 for x in cell_type_dict[compartment]])
        n_na = sum([x in cells_na for x in cell_type_dict[compartment]])
    
        cell_type_n = n_00 + n_01 + n_11 + n_na
        cell_type_cov = n_00 + n_01 + n_11
        cell_type_mut = n_01 + n_11
    
        binom_result = binomtest(k = cell_type_mut, n = cell_type_cov, p = mutation_prop[m], alternative = "greater")
        binom_p = binom_result.pvalue
  
        output_line_list = [m] + compartment.split(":") + [x for x in [cell_type_n, cell_type_cov, n_00, n_01, n_11, n_na, cell_type_mut, cell_type_mut/cell_type_cov, mutation_prop[m], binom_p]]

        for k,v in zip(output_header_list, output_line_list):
            output_list_dict[k].append(v)

# BH adjusted p-value
adjusted_pvalues = false_discovery_control(output_list_dict["binomial_pvalue"], method = "bh")
output_list_dict["adjusted_pvalue"] = adjusted_pvalues

output_file = open(output_filename, "w")
output_file.write("\t".join(output_header_list + ["adjusted_pvalue"]) + "\n")

n_lines = len(adjusted_pvalues)
for i in range(n_lines):
    output_line_list = [str(output_list_dict[x][i]) for x in output_header_list + ["adjusted_pvalue"]]
    output_file.write("\t".join(output_line_list) + "\n")

ds.close()
output_file.close()

