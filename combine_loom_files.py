import sys
import loompy
import numpy as np
import pandas as pd

# Usage: python combine_loom_files.py loom1,loom2 prefix1,prefix2 output.loom
# Combines two or more loom files by finding overlapping variants and pasting columns together
# Each barcode is prefixed with an identifier (prefix)

# Read command line options

loom_file_arg = sys.argv[1]
loom_file_list = loom_file_arg.split(",")

prefix_arg = sys.argv[2]
prefix_list = prefix_arg.split(",")

output_file = sys.argv[3]

# Pre-process each loom file separately

n_loom_files = len(loom_file_list)
new_barcode_list = []
sorted_id_list = []
sorted_id_list_dict = {}
barcode_indexes = []

for i in range(n_loom_files):
    loom_file_path = loom_file_list[i]
    barcode_prefix = prefix_list[i]
    ds = loompy.connect(loom_file_path)
    new_barcode = list(barcode_prefix + "_" + ds.ca["barcode"])
    new_barcode_list.extend(new_barcode)
    sorted_id_list.append(sorted(ds.ra["id"]))
    sorted_id_list_dict[i] = {x:'' for x in sorted_id_list[i]}
    barcode_indexes.append(len(ds.ca["barcode"]))
    ds.close()

# Remove duplicate rows

keep_rows_rm_dup = []
for i in range(n_loom_files):
    keep_rows = []
    for index in range(len(sorted_id_list[i])):
        n_matches = 0
        for j in range(n_loom_files):
            if sorted_id_list[i][index] in sorted_id_list_dict[j]:
                n_matches += 1
        if n_matches == n_loom_files:
            keep_rows.append(index)
    ds_ids_df = pd.DataFrame(sorted_id_list[i])
    indexes_of_duplicates = ds_ids_df[ds_ids_df.duplicated(keep = "first")].index.values.tolist()
    for k in indexes_of_duplicates:
        if k in keep_rows:
            keep_rows.remove(k)
    keep_rows_rm_dup.append(np.array(keep_rows))

# Write filtered looms to file, separately

### from https://linnarssonlab.org/loompy/cookbook/index.html#combining-data-using-scan-and-new
### and https://github.com/linnarsson-lab/loompy/issues/59#issuecomment-406060518
### Note: when sorting by key, boolean values to keep rows must also be sorted by same key

with loompy.new(output_file) as dsout:  # Create a new, empty, loom file
    for i in range(n_loom_files):
        f = loom_file_list[i]
        with loompy.connect(f) as ds:
            for (ix, selection, view) in ds.scan(axis = 1, key = "id"):
                dsout.add_columns(layers = view.layers[keep_rows_rm_dup[i], :],
                        col_attrs = view.ca,
                        row_attrs = view.ra[keep_rows_rm_dup[i]])
    dsout.ca["barcode"] = new_barcode_list
    #del dsout.ra["CLNDEF"]
    #del dsout.ra["COMMON"]
    #del dsout.ra["QUAL"]
    #del dsout.ra["RSID"]
    #del dsout.ra["nonref_cells"]

# Check that everything worked as expected

ds = loompy.connect(output_file)

for i in range(n_loom_files):
    loom_file_path = loom_file_list[i]
    lf = loompy.connect(loom_file_path)
    if np.mean(np.array(sorted(lf.ra["id"]))[keep_rows_rm_dup[i]] == ds.ra["id"]) != 1:
        exit("Combined loom file does not match input file " + loom_file_path)
    lf.close()

ds.close()

