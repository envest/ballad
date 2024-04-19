import sys
import loompy
import numpy as np

loom_filename = sys.argv[1]
output_file_prefix = sys.argv[2]
output_dir = sys.argv[3]

gq_min = 30
dp_min = 10
scVAF_min = 0.2
max_cell_NA = 0.5
max_var_NA = 0.5
min_pct_mut = 0.01

var_cb_file = output_dir + "/" + output_file_prefix + ".variant_cell_barcode.tsv"
var_sum_file = output_dir + "/" + output_file_prefix + ".variant_summary.tsv"

ds = loompy.connect(loom_filename)
ds_GT_copy = ds[""][:,:]

keep_var_index = []
keep_cb_index = []

n_variants = ds.shape[0]
n_cells = ds.shape[1]

# genotype filter

def calculate_scVAF(GT, AD, DP):
    scVAF = np.array([[1.0]*GT.shape[1]]*GT.shape[0])
    where_GT_is_01 = np.where((GT[:,:] == 1) & (DP[:,:] > 0))
    i = where_GT_is_01[0]
    j = where_GT_is_01[1]
    scVAF[i,j] = np.divide(AD[:,:][i,j], DP[:,:][i,j])
    return(scVAF)

scVAF = calculate_scVAF(ds, ds["AD"], ds["DP"])

set_missing_index_array = np.where((ds["GQ"][:,:] < gq_min) | (ds["DP"][:,:] < dp_min) | ((ds[:,:] == 1) & (scVAF[:,:] < scVAF_min)) | ((ds[:,:] == 1) & (scVAF[:,:] > 1 - scVAF_min)))
i = set_missing_index_array[0]
j = set_missing_index_array[1]

ds_GT_copy[i,j] = 3

# cell barcode filter

keep_cb_index = np.where(np.apply_along_axis(sum, 0, ds_GT_copy == 3)/n_variants < max_cell_NA)[0]

# variant filter

def prop_GT_mutated(n_GT_00, n_GT_not_missing):
    prop_GT_00 = np.array([0.0]*len(n_GT_00))
    has_GT = np.where(n_GT_not_missing > 0)
    prop_GT_00[has_GT] = n_GT_00[has_GT]/n_GT_not_missing[has_GT]
    return(1 - prop_GT_00)

n_GT_00 = np.apply_along_axis(sum, 1, ds_GT_copy[:, keep_cb_index] == 0)
n_GT_not_missing = np.apply_along_axis(sum, 1, ds_GT_copy[:, keep_cb_index] != 3)
prop_var_mutated = prop_GT_mutated(n_GT_00, n_GT_not_missing)
prop_var_NA = np.apply_along_axis(np.mean, 1, ds_GT_copy[:, keep_cb_index] == 3)

keep_var_index = np.where((prop_var_mutated > min_pct_mut) & (prop_var_NA < max_var_NA))[0]

# write output files

var_cb = open(var_cb_file, "w")
var_cb_header_list = ["chr", "pos", "ref", "alt", "cell_barcode", "GT", "DP", "AD_REF", "AD_ALT", "GQ"]
var_cb_header_line = "\t".join(var_cb_header_list) + "\n"
var_cb.write(var_cb_header_line)

var_sum = open(var_sum_file, "w")
var_sum_header_list = ["chr", "pos", "ref", "alt", "n_cells", "n_00", "n_01", "n_11", "n_na", "prop_cells_mutated", "alt_allele_freq"]
var_sum_header_line = "\t".join(var_sum_header_list) + "\n"
var_sum.write(var_sum_header_line)

for var_index in keep_var_index:
    
    chr = ds.ra["CHROM"][var_index]
    pos = ds.ra["POS"][var_index]
    ref = ds.ra["REF"][var_index]
    alt = ds.ra["ALT"][var_index]

    for cb_index in keep_cb_index:
        
        cb = ds.ca["barcode"][cb_index]
        gt = ds_GT_copy[var_index, cb_index]
        dp = ds["DP"][var_index, cb_index]
        ad_ref = ds["RO"][var_index, cb_index]
        ad_alt = ds["AD"][var_index, cb_index]
        gq = ds["GQ"][var_index, cb_index]

        var_cb_list = [str(x) for x in [chr, pos, ref, alt, cb, gt, dp, ad_ref, ad_alt, gq]]
        var_line = "\t".join(var_cb_list) + "\n"
        var_cb.write(var_line)

    var = ds_GT_copy[var_index, keep_cb_index]
    n_cells = len(var)
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

