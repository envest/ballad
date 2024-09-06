import sys
import loompy
import numpy as np
import random

# Usage: python downsample_loom_cells.py proportion seed input.loom downsampled_cells.loom
# Downsample cells from a loom file
# Proportion indicates what proportion of cells to keep
# Set seed for reproducibility

# Read command line options

loom_file_path = sys.argv[1]
keep_proportion = float(sys.argv[2])
random_seed = sys.argv[3]
output_file = sys.argv[4]

# check keep_proportion is (0,1]
if keep_proportion <= 0 or keep_proportion > 1:
    sys.exit("Keep proportion should be (0,1]: greater than 0 and less than or equal to 1.")

# read in loom file
ds = loompy.connect(loom_file_path)

# cells to keep
random.seed(random_seed)
cells_list = ds.ca["barcode"].tolist()
n_cells_total = len(cells_list) 
n_cells_keep = round(n_cells_total*keep_proportion)
cells_keep = random.sample(cells_list, n_cells_keep)
cell_indexes_keep = np.where([cell in cells_keep for cell in cells_list])[0] # Select the cells in keep list

### from https://linnarssonlab.org/loompy/cookbook/index.html#combining-data-using-scan-and-new
### and https://github.com/linnarsson-lab/loompy/issues/59#issuecomment-406060518

with loompy.new(output_file) as dsout:  # Create a new, empty, loom file
    for (ix, selection, view) in ds.scan(items = cell_indexes_keep, axis = 1):
        dsout.add_columns(view.layers, col_attrs = view.ca, row_attrs = view.ra)

# close input loom file
ds.close()

