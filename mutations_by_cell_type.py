# Analyze mutations by cell type
# python mutations_by_cell_type.py sample.loom mutations.tsv cell_barcodes.tsv output_dir output_prefix

import os
import sys
import loompy
import numpy as np
from scipy.stats import binomtest

loom_filename = sys.argv[1]
mutations_filename = sys.argv[2]
cell_barcodes_filename = sys.argv[3]
output_dir = sys.argv[4]
output_prefix = sys.argv[5]

output_filename = os.path(output_dir, output_prefix + ".mutations_by_cell_type.tsv")

if not output_dir.isdir():
  sys.exit(output_dir + " directory does not exist")

cell_type_dict = {}
with f as open(cell_barcodes_filename, "r"):
  line = f.readline()
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
  loom_mutation_list.append(":".join([ds.ra["CHROM"][i], str(ds.ra["POS"][i]), ds.ra["REF"][i], ds.ra["ALT"][i])

with f as open(mutations_filename, "r"):
  m = f.readline()
  if m in loom_mutation_list:
    m_index = np.where(loom_mutation_list == m)[0][0]
    mutation_dict[m] = m_index
    n_00 = sum(ds[m_index,] == 0)
    n_01 = sum(ds[m_index,] == 1)
    n_11 = sum(ds[m_index,] == 2)
    mutation_prop[m] = (n_01 + n_11)/(n_00 + n_01 + n_11)
  else:
    sys.message("Mutation " + m + " not in loom file")
   
  
output_header_list = ["annotation_level", "tumor_normal", "level_one", "level_two", "cell_type_total", "cell_type_with_coverage", "cell_type_missing", "n_00", "n_01", "n_11", "cell_type_total_mutated", "cell_type_proportion_mutated", "binomial_pvalue", "binomial_statistic"]
output_file = open(output_filename, "w")
output_file.write("\t".join(output_header_list) + "\n")

for m,i in mutation_dict.items():
  
  cells_00 = ds.ca["barcode"][np.where(ds[i,] == 0)]
  cells_01 = ds.ca["barcode"][np.where(ds[i,] == 1)]
  cells_11 = ds.ca["barcode"][np.where(ds[i,] == 2)]
  cells_na = ds.ca["barcode"][np.where(ds[i,] == 3)]
  
  for compartment in cell_type_dict.keys():
    
    n_00 = sum(cell_type_dict[compartment] in cells_00)
    n_01 = sum(cell_type_dict[compartment] in cells_01)
    n_11 = sum(cell_type_dict[compartment] in cells_11)
    n_na = sum(cell_type_dict[compartment] in cells_na)
    
    cell_type_n = n_00 + n_01 + n_11 + n_na
    cell_type_cov = n_00 + n_01 + n_11
    cell_type_mut = n_01 + n_11
    
    binom_result = binomtest(k = cell_type_mut, n = cell_type_cov, p = mutation_prop[m], alternative = "greater")
    binom_p = binom_result.pvalue()
    binom_stat = binom_result.statistic()
  
    output_line_list = compartment.split(":") + [str(x) for x in [cell_type_n, cell_type_cov, cell_type_na, n_00, n_01, n_11, cell_type_mut, cell_type_mut/cell_type_cov, binom_p, binom_stat]]

    output_file.write("\t".join(output_line_list) + "\t")

ds.close()
output_file.close()

