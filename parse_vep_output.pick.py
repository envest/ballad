import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]

input_file = open(input_filename, "r")
output_file = open(output_filename, "w")

for line in input_file:

    if line.startswith("##INFO"):
        INFO_header_list = line.strip('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: ').strip('">\n').split("|")
        header_list = ["CHROM", "POS", "REF", "ALT"]
        header_list.extend(INFO_header_list)
        output_file.write("\t".join(header_list) + "\n")
        next

    if not line.startswith("#"):
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT = line.strip().split("\t")

        if INFO == ".":

            INFO_list = ["NA"]*len(INFO_header_list)

        else:

            INFO_list = ["NA" if a == '' else a for a in INFO.split("=")[1].split("|")]
        
        output_list = [CHROM, POS, REF, ALT]
        output_list.extend(INFO_list)
        output_list_replace_na = [x.replace(".", "NA") for x in output_list]
        output_file.write("\t".join(output_list) + "\n")

input_file.close()
output_file.close()
