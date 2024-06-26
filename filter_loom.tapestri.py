import sys
import shutil
import loompy
import numpy as np

loom_filename = sys.argv[1]
output_file_prefix = sys.argv[2]
output_dir = sys.argv[3]

gq_min = 30
dp_min = 10
scVAF_min = 0.2
min_prop_cell_with_GT = 0.5
min_prop_var_with_GT = 0.5
min_prop_var_mutated = 0.01

var_cb_file = output_dir + "/" + output_file_prefix + ".variant_cell_barcode.tsv"
var_sum_file = output_dir + "/" + output_file_prefix + ".variant_summary.tsv"
filtered_loom_file = output_dir + "/" + output_file_prefix + ".loom"

# copy existing loom file to filtered loom file
shutil.copy2(loom_filename, filtered_loom_file)
ds = loompy.connect(filtered_loom_file)

keep_var_index = []
keep_cb_index = []

n_variants = ds.shape[0]
n_cells = ds.shape[1]

# Filters 1-3: set genotype to missing if GQ < gq_min, DP < dp_min, or scVAF < scVAF_min
# applied simultaneously because each genotype is independent

def calculate_scVAF(GT, AD, DP):
    scVAF = np.array([[None]*GT.shape[1]]*GT.shape[0])
    where_GT_is_01 = np.where((GT[:,:] == 1 | GT[:,:] == 2) & (DP[:,:] > 0))
    i = where_GT_is_01[0]
    j = where_GT_is_01[1]
    scVAF[i,j] = np.divide(AD[:,:][i,j], DP[:,:][i,j])
    return(scVAF)

scVAF = calculate_scVAF(ds, ds["AD"], ds["DP"])

set_missing_index_array = np.where((ds["GQ"][:,:] < gq_min) | (ds["DP"][:,:] < dp_min) | (scVAF[:,:] is not None & scVAF[:,:] < scVAF_min))
i = set_missing_index_array[0]
j = set_missing_index_array[1]

ds[i,j] = 3

# Filter 4: remove variants with genotype in < X% of cells
prop_var_with_GT = 1 - np.apply_along_axis(np.mean, 1, ds[:,:] == 3)
keep_var_index_F4 = np.where(prop_var_with_GT >= min_prop_var_with_GT)
ds = ds[:,keep_var_index_F4]

# Filter 5: remove cells with genotype in < X% of variants
prop_cell_with_GT = 1 - np.apply_along_axis(np.mean, 0, ds[:,:] == 3)
keep_cell_index_F5 = np.where(prop_cell_with_GT >= min_prop_cell_with_GT)
ds = ds[keep_cell_index_F5,:]

# Filter 6: remove variants mutated in < X% of cells
def prop_GT_mutated(n_GT_00, n_GT_not_missing):
    prop_GT_00 = np.array([0.0]*len(n_GT_00))
    has_GT = np.where(n_GT_not_missing > 0)
    prop_GT_00[has_GT] = n_GT_00[has_GT]/n_GT_not_missing[has_GT]
    return(1 - prop_GT_00)

n_GT_00 = np.apply_along_axis(sum, 1, ds[:,:] == 0)
n_GT_not_missing = np.apply_along_axis(sum, 1, ds[:,:] != 3)
prop_var_mutated = prop_GT_mutated(n_GT_00, n_GT_not_missing)

keep_var_index_F6 = np.where(prop_var_mutated >= min_prop_var_mutated)
ds = ds[:,keep_var_index_F6]

# write output files

var_cb = open(var_cb_file, "w")
var_cb_header_list = ["chr", "pos", "ref", "alt", "cell_barcode", "GT", "DP", "AD_REF", "AD_ALT", "GQ"]
var_cb_header_line = "\t".join(var_cb_header_list) + "\n"
var_cb.write(var_cb_header_line)

var_sum = open(var_sum_file, "w")
var_sum_header_list = ["chr", "pos", "ref", "alt", "n_cells", "n_00", "n_01", "n_11", "n_na", "prop_cells_mutated", "alt_allele_freq"]
var_sum_header_line = "\t".join(var_sum_header_list) + "\n"
var_sum.write(var_sum_header_line)

n_var = ds.shape[0]
n_cells = ds.shape[1]

for var_index in range(n_var):
    
    chr = ds.ra["CHROM"][var_index]
    pos = ds.ra["POS"][var_index]
    ref = ds.ra["REF"][var_index]
    alt = ds.ra["ALT"][var_index]

    for cb_index in range(n_cells):
        
        cb = ds.ca["barcode"][cb_index]
        gt = ds_GT_copy[var_index, cb_index]
        dp = ds["DP"][var_index, cb_index]
        ad_ref = ds["RO"][var_index, cb_index]
        ad_alt = ds["AD"][var_index, cb_index]
        gq = ds["GQ"][var_index, cb_index]

        var_cb_list = [str(x) for x in [chr, pos, ref, alt, cb, gt, dp, ad_ref, ad_alt, gq]]
        var_line = "\t".join(var_cb_list) + "\n"
        var_cb.write(var_line)

    var = ds[var_index,:]
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
ds.close()
