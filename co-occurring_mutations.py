import sys
import loompy
import numpy as np

loom_filename = sys.argv[1]
mutations = sys.argv[2]
output_filename = sys.argv[3]

ds = loompy.connect(loom_filename)

n_variants = ds.shape[0]
n_barcodes = ds.shape[1]

variant_list = []
for i in range(n_variants):
    new_variant = ":".join(["chr" + ds.ra["CHROM"][i], str(ds.ra["POS"][i]), ds.ra["REF"][i], ds.ra["ALT"][i]])
    variant_list.append(new_variant)

mutation_list = mutations.split(",")

for m in mutation_list:
    if m not in variant_list:
        print("Mutation " + m + " not in variant list.")

mutation_dict = {m : None for m in mutation_list}
allele_dict = {}

for m in mutation_list:
    m_index = np.where([x == m for x in variant_list])[0][0]
    mutation_dict[m] = ds[m_index,]

output_file = open(output_filename, "w")
output_file.write("\t".join(["barcode"] + mutation_list) + "\n")

for barcode_index in range(n_barcodes):
    bc = ds.ca["barcode"][barcode_index]
    print_list = [bc]
    allele_list = []
    for m in mutation_list:
        allele_list.append(str(mutation_dict[m][barcode_index]))
        alleles = ",".join([str(x) for x in allele_list])
    if alleles in allele_dict:
        allele_dict[alleles] += 1
    else:
        allele_dict[alleles] = 1
    print_list.extend(allele_list)
    output_file.write("\t".join(print_list) + "\n")

output_file.close()

print("\t".join(["barcode"] + mutation_list))

for allele_index in reversed(np.argsort(list(allele_dict.values()))):
    k = list(allele_dict.keys())[allele_index]
    print("\t".join([str(allele_dict[k])] + k.split(",")))



