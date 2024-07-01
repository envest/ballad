import sys
import loompy
from liftover import get_lifter

converter = get_lifter('hg19', 'hg38', one_based=True)

loom_filename = sys.argv[1]
output_file_prefix = sys.argv[2]
output_dir = sys.argv[3]

output_filename = output_dir + "/" + output_file_prefix + ".coverage_sites_liftover.tsv"
vep_input_hg19_filename = output_dir + "/" + output_file_prefix + ".VEP_input.hg19.vcf"
vep_input_hg38_filename = output_dir + "/" + output_file_prefix + ".VEP_input.hg38.vcf"

output_header_list = ["hg19_chr", "hg19_pos", "hg38_chr", "hg38_pos", "REF", "ALT"]
vcf_header_list = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]

output_file = open(output_filename, "w")
output_file.write("\t".join(output_header_list) + "\n")

vep_input_hg19_file = open(vep_input_hg19_filename, "w")
vep_input_hg19_file.write("##fileformat=VCFv4.2\n")
for c in [str(x) for x in range(1, 23)] + ["X", "Y"]:
    vep_input_hg19_file.write("##contig=<ID=chr" + c + ">\n")
vep_input_hg19_file.write("\t".join(vcf_header_list) + "\n")

vep_input_hg38_file = open(vep_input_hg38_filename, "w")
vep_input_hg38_file.write("##fileformat=VCFv4.2\n")
for c in [str(x) for x in range(1, 23)] + ["X", "Y"]:
    vep_input_hg38_file.write("##contig=<ID=chr" + c + ">\n")
vep_input_hg38_file.write("\t".join(vcf_header_list) + "\n")

ds = loompy.connect(loom_filename)

loom_hg19_coordinates = zip(ds.ra["CHROM"], ds.ra["POS"], ds.ra["REF"], ds.ra["ALT"])

for hg19_chr, hg19_pos, ref, alt in loom_hg19_coordinates:

    try:
        hg38_chr, hg38_pos, strand = converter[hg19_chr][hg19_pos][0]
        hg38_chr_rm = hg38_chr.strip("chr")
    except:
        hg38_chr_rm = "NA"
        hg38_pos = "NA"

    output_line_list = [str(x) for x in [hg19_chr, hg19_pos, hg38_chr_rm, hg38_pos, ref, alt]]
    vep_input_hg19_line_list = [str(x) for x in [hg19_chr, hg19_pos, ".", ref, alt, ".", ".", ".", "."]]
    vep_input_hg38_line_list = [str(x) for x in [hg38_chr, hg38_pos, ".", ref, alt, ".", ".", ".", "."]]

    output_file.write("\t".join(output_line_list) + "\n")
    vep_input_hg19_file.write("\t".join(vep_input_hg19_line_list) + "\n")
    
    if hg38_chr_rm != "NA":
        vep_input_hg38_file.write("\t".join(vep_input_hg38_line_list) + "\n")

output_file.close()
vep_input_hg19_file.close()
vep_input_hg38_file.close()

