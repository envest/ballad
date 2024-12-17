import sys
import loompy

loom_filename = sys.argv[1]

ds = loompy.connect(loom_filename)

n_cells = len(ds.ca["barcode"])

print(n_cells)

ds.close()
