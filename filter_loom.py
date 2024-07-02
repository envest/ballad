import sys
import loompy
import numpy as np

loom_filename = sys.argv[1]
output_file_prefix = sys.argv[2]
output_dir = sys.argv[3]

liftover_file = sys.argv[4]
annotation_loom_file = sys.argv[5]
rna_loom_file = sys.argv[6]
#atac_loom_file = sys.argv[6]

gq_min = 30
dp_min = 10
scVAF_min = 0.2
min_prop_cell_with_GT = 0.5
min_prop_var_with_GT = 0.5
min_prop_var_mutated = 0.01

var_cb_file = output_dir + "/" + output_file_prefix + ".variant_cell_barcode.tsv"
var_sum_file = output_dir + "/" + output_file_prefix + ".variant_summary.tsv"
filtered_loom_file = output_dir + "/" + output_file_prefix + ".loom"

liftover_dict = {"GRCh37" : {}, "GRCh38" : {}}
with open(liftover_file, "r") as lo:
  for line in lo.readlines():
    lo_list = line.strip().split()
    if lo_list[2] != "NA" and len(lo_list) == 6:
      hg19_chr, hg19_pos, hg38_chr, hg38_pos, REF, ALT = lo_list
      GRCh37_variant = ":".join(["chr" + hg19_chr, str(hg19_pos), REF, ALT])
      GRCh38_variant = ":".join(["chr" + hg38_chr, str(hg38_pos), REF, ALT])
      liftover_dict["GRCh37"][GRCh37_variant] = GRCh38_variant
      liftover_dict["GRCh38"][GRCh38_variant] = GRCh37_variant

# set up dictionaries to modify/filter
layers_dict = {}
row_dict = {}
col_dict = {}

# connect to pre-filtered loom file
with loompy.connect(loom_filename) as ds:
  for k,v in ds.layers.items():
    layers_dict[k] = v[:,:]
  for k,v in ds.ra.items():
    row_dict[k] = v
  for k,v in ds.ca.items():
    col_dict[k] = v

row_dict["variant_GRCh37"] = [":".join(["chr" + str(c), str(p), r, a]) for c,p,r,a in 
  zip(row_dict["CHROM"], 
    row_dict["POS"], 
    row_dict["REF"], 
    row_dict["ALT"])]
    
row_dict["variant_GRCh38"] = [liftover_dict["GRCh37"][x] if x in liftover_dict["GRCh37"] else None for x in row_dict["variant_GRCh37"]]

# set up annotation, RNA, and ATAC dictionaries
with loompy.connect(annotation_loom_file) as ds_annot:
  annot_dict = {}
  for k in ds_annot.ra.keys():
    annot_dict[k] = []
  for v in row_dict["variant_GRCh38"]:
    annot_index = None
    if v is not None and v in ds_annot.ra["variant"]:
      annot_index = np.where(ds_annot.ra["variant"] == v)[0][0]
    for k in ds_annot.ra.keys():
      if annot_index is None:
        annot_dict[k].append(None)
      else:
        annot_dict[k].append(ds_annot.ra[k][annot_index])

# how many cells are REF/ALT from RNA
with loompy.connect(rna_loom_file) as ds_rna:
  rna_n_cells_ref = np.zeros(ds_rna.shape[0])
  rna_n_cells_alt = np.zeros(ds_rna.shape[0])
  for(ix, selection, view) in ds_rna.scan(axis = 0):
    rna_n_cells_ref[selection.min():selection.max() + 1] = np.apply_along_axis(sum, 1, view[:,:] == 1)
    rna_n_cells_alt[selection.min():selection.max() + 1] = np.apply_along_axis(sum, 1, view[:,:] > 1)
  rna_n_cells_missing = ds_rna.shape[1] - rna_n_cells_ref - rna_n_cells_alt

# how many cells are REF/ALT from ATAC  
with loompy.connect(atac_loom_file) as ds_atac:
  atac_n_cells_ref = np.zeros(ds_atac.shape[0])
  atac_n_cells_alt = np.zeros(ds_atac.shape[0])
  for(ix, selection, view) in ds_atac.scan(axis = 0):
    atac_n_cells_ref[selection.min():selection.max() + 1] = np.apply_along_axis(sum, 1, view[:,:] == 1)
    atac_n_cells_alt[selection.min():selection.max() + 1] = np.apply_along_axis(sum, 1, view[:,:] > 1)
  atac_n_cells_missing = ds_atac.shape[1] - atac_n_cells_ref - atac_n_cells_alt

# gather RNA information about each variant
with loompy.connect(rna_loom_file) as ds_rna:  
  rna_dict = {"n_cells_REF" : np.zeros(), "n_cells_ALT" : [], "n_cells_missing" : []}
  for k in ds_rna.ra.keys():
    rna_dict[k] = []
  for v in row_dict["variant_GRCh38"]:
    rna_index = None
    if v is not None and v in ds_rna.ra["variant"]:
      rna_index = np.where(ds_rna.ra["variant"] == v)[0][0]
      n_cells_REF = rna_n_cells_ref[rna_index]
      n_cells_ALT = rna_n_cells_alt[rna_index]
      n_cells_missing = rna_n_cells_missing[rna_index]
    else:
      n_cells_REF = -1
      n_cells_ALT = -1
      n_cells_missing = -1
    rna_dict["n_cells_REF"].append(n_cells_REF)
    rna_dict["n_cells_ALT"].append(n_cells_ALT)
    rna_dict["n_cells_missing"].append(n_cells_missing)
    for k in ds_rna.ra.keys():
      if rna_index is None:
        rna_dict[k].append(None)
      else:
        rna_dict[k].append(ds_rna.ra[k][rna_index])
        
rna_dict["n_cells_ALT"] = np.array(rna_dict["n_cells_ALT"])

# gather ATAC information about each variant
with loompy.connect(atac_loom_file) as ds_atac:  
  atac_dict = {"n_cells_REF" : np.zeros(), "n_cells_ALT" : [], "n_cells_missing" : []}
  for k in ds_atac.ra.keys():
    atac_dict[k] = []
  for v in row_dict["variant_GRCh38"]:
    atac_index = None
    if v is not None and v in ds_atac.ra["variant"]:
      atac_index = np.where(ds_atac.ra["variant"] == v)[0][0]
      n_cells_REF = atac_n_cells_ref[atac_index]
      n_cells_ALT = atac_n_cells_alt[atac_index]
      n_cells_missing = atac_n_cells_missing[atac_index]
    else:
      n_cells_REF = -1
      n_cells_ALT = -1
      n_cells_missing = -1
    atac_dict["n_cells_REF"].append(n_cells_REF)
    atac_dict["n_cells_ALT"].append(n_cells_ALT)
    atac_dict["n_cells_missing"].append(n_cells_missing)
    for k in ds_atac.ra.keys():
      if atac_index is None:
        atac_dict[k].append(None)
      else:
        atac_dict[k].append(ds_atac.ra[k][atac_index])
        
atac_dict["n_cells_ALT"] = np.array(atac_dict["n_cells_ALT"])

# Filters 1-3: set genotype to missing if GQ < gq_min, DP < dp_min, or scVAF < scVAF_min
# can be applied simultaneously because each genotype is independent

def calculate_scVAF(GT, AD, DP):
    scVAF = np.array([[1.0]*GT.shape[1]]*GT.shape[0])
    where_GT_is_12 = np.where(np.logical_and(list(GT == 1) or list(GT == 2), list(DP > 0)))
    i = where_GT_is_12[0]
    j = where_GT_is_12[1]
    scVAF[i,j] = np.divide(AD[i,j], DP[i,j])
    return(scVAF)

scVAF = calculate_scVAF(layers_dict[""], layers_dict["AD"], layers_dict["DP"])

set_missing_index_array = np.where((layers_dict["GQ"] < gq_min) | (layers_dict["DP"] < dp_min) | (scVAF < scVAF_min))

i = set_missing_index_array[0]
j = set_missing_index_array[1]

# set genotypes to missing
layers_dict[""][i,j] = 3

# Filter 4: remove variants with genotype in < X% of cells
prop_var_with_GT = 1 - np.apply_along_axis(np.mean, 1, layers_dict[""] == 3)
keep_var_index_F4 = np.where(prop_var_with_GT >= min_prop_var_with_GT)[0]

for k,v in layers_dict.items():
  layers_dict[k] = v[keep_var_index_F4,:]

for k,v in row_dict.items():
  row_dict[k] = v[keep_var_index_F4]
  
for k,v in rna_dict.items():
  rna_dict[k] = v[keep_var_index_F4]
  
#for k,v in atac_dict.items():
#  atac_dict[k] = v[keep_var_index_F4]

# Filter 5: remove cells with genotype in < X% of variants
prop_cell_with_GT = 1 - np.apply_along_axis(np.mean, 0, layers_dict[""] == 3)
keep_cell_index_F5 = np.where(prop_cell_with_GT >= min_prop_cell_with_GT)[0]

for k,v in layers_dict.items():
  layers_dict[k] = v[:,keep_cell_index_F5]

for k,v in col_dict.items():
  col_dict[k] = v[keep_cell_index_F5]

# Filter 6: remove variants mutated in < X% of cells
def prop_GT_mutated(n_GT_00, n_GT_not_missing):
    prop_GT_00 = np.array([0.0]*len(n_GT_00))
    has_GT = np.where(n_GT_not_missing > 0)
    prop_GT_00[has_GT] = n_GT_00[has_GT]/n_GT_not_missing[has_GT]
    return(1 - prop_GT_00)

n_GT_00 = np.apply_along_axis(sum, 1, layers_dict[""] == 0)
n_GT_not_missing = np.apply_along_axis(sum, 1, layers_dict[""] != 3)
prop_var_mutated = prop_GT_mutated(n_GT_00, n_GT_not_missing)

multiome_evidence = np.logical_or(rna_dict["n_cells_ALT"] > 0, atac_dict["n_cells_ALT"] > 0)
keep_var_index_F6 = np.where(prop_var_mutated >= min_prop_var_mutated or multiome_evidence)[0]

for k,v in layers_dict.items():
  layers_dict[k] = v[keep_var_index_F6,:]

for k,v in row_dict.items():
  row_dict[k] = v[keep_var_index_F6]

# write filtered loom file
loompy.create(filtered_loom_file, layers_dict, row_dict, col_dict)

# write other output files

var_cb = open(var_cb_file, "w")
var_cb_header_list = ["chr", "pos", "ref", "alt", "cell_barcode", "GT", "DP", "AD_REF", "AD_ALT", "GQ"]
var_cb_header_line = "\t".join(var_cb_header_list) + "\n"
var_cb.write(var_cb_header_line)

var_sum = open(var_sum_file, "w")
var_sum_header_list = ["chr", "pos", "ref", "alt", "n_cells", "n_00", "n_01", "n_11", "n_na", "prop_cells_mutated", "alt_allele_freq"]
var_sum_header_line = "\t".join(var_sum_header_list) + "\n"
var_sum.write(var_sum_header_line)

n_var = layers_dict[""].shape[0]
n_cells = layers_dict[""].shape[1]

for var_index in range(n_var):
    chr = row_dict["CHROM"][var_index]
    pos = row_dict["POS"][var_index]
    ref = row_dict["REF"][var_index]
    alt = row_dict["ALT"][var_index]
    for cb_index in range(n_cells):
        cb = col_dict["barcode"][cb_index]
        gt = layers_dict[""][var_index, cb_index]
        dp = layers_dict["DP"][var_index, cb_index]
        ad_ref = layers_dict["RO"][var_index, cb_index]
        ad_alt = layers_dict["AD"][var_index, cb_index]
        gq = layers_dict["GQ"][var_index, cb_index]
        var_cb_list = [str(x) for x in [chr, pos, ref, alt, cb, gt, dp, ad_ref, ad_alt, gq]]
        var_line = "\t".join(var_cb_list) + "\n"
        var_cb.write(var_line)
    var = layers_dict[""][var_index,:]
    n_na = (var == 3).sum()
    n_not_na = n_cells - n_na
    n_00 = (var == 0).sum()
    n_01 = (var == 1).sum()
    n_11 = (var == 2).sum()
    prop_cells_mutated = (n_01 + n_11)/(n_not_na)
    alt_allele_freq = (n_01 + 2*n_11)/(2*n_not_na)
    sum_line = [str(x) for x in [chr, pos, ref, alt, str(n_cells), str(n_00), str(n_01), str(n_11), str(n_na), str(prop_cells_mutated), str(alt_allele_freq)]]
    var_sum.write("\t".join(sum_line) + "\n")

var_cb.close()
var_sum.close()
