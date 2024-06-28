import sys
import loompy
import numpy as np

loom_filename = sys.argv[1]
bc_filename = sys.argv[2]
selected_cells_loom_filename = sys.argv[3]

# barcodes
bc_file = open(bc_filename, "r")
bc_list = bc_file.readlines().strip().split()

# create new loom file and select cells to fill it
with loompy.new(selected_cells_loom_filename) as selected_cells_ds:
  with loompy.connect(loom_filename) as ds:
    keep_barcodes_index = np.where([x in bc_list for x in ds.ca["barcode"]])[0]
    for (ix, selection, view) in ds.scan(items = keep_barcodes_index, axis = 1):
      selected_cells_ds.add_columns(view.layers, col_attrs = view.ca, row_attrs = view.ra)
