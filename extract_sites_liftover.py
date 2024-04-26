import sys
import loompy
from liftover import get_lifter

converter = get_lifter('hg19', 'hg38', one_based=True)

loom_filename = sys.argv[1]
output_file_prefix = sys.argv[2]
output_dir = sys.argv[3]

output_filename = output_dir + "/" + output_file_prefix + ".coverage_sites.tsv"
output_file = open(output_filename, "w")
output_header_list = ["hg19_chr", "hg19_pos", "hg38_chr", "hg38_pos", "REF", "ALT", "strand"]
output_file.write("\t".join(output_header_list) + "\n")

ds = loompy.connect(loom_filename)

print(ds.shape)

loom_hg19_coordinates = zip(ds.ra["CHROM"], ds.ra["POS"], ds.ra["REF"], ds.ra["ALT"])

for hg19_chr, hg19_pos, ref, alt in loom_hg19_coordinates:

    try:
        hg38_chr, hg38_pos, strand = converter[hg19_chr][hg19_pos][0]
    except:
        hg38_chr = "NA"
        hg38_pos = "NA"

    output_line_list = [str(x) for x in ["chr"+hg19_chr, hg19_pos, hg38_chr, hg38_pos, ref, alt, "+"]]

    output_file.write("\t".join(output_line_list) + "\n")

output_file.close()

