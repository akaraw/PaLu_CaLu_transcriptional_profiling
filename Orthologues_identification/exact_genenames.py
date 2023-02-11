def extract_genenames(input_file, bat_gtf_file, human_gtf_file):
    gene_names = {}
    attributes_bat_dict = {}
    attributes_human_dict = {}
    with open(bat_gtf_file) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "CDS":
                continue
            attributes = fields[-1].split(";")
            for attribute in attributes:
                parts = attribute.strip().split(" ",1)
                if len(parts) !=2:
                    continue
                key, value = parts
                key = key.strip()
                value = value.strip("\"")
                if key == "gene_id":
                    gene_id = value
                if key == "protein_id":
                    protein_id = value
            attributes_bat_dict[protein_id] = gene_id
    with open(human_gtf_file) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "CDS":
                continue
            attributes = fields[-1].split(";")
            for attribute in attributes:
                parts = attribute.strip().split(" ",1)
                if len(parts) !=2:
                    continue
                key, value = parts
                key = key.strip()
                value = value.strip("\"")
                if key == "gene_id":
                    gene_id = value
                if key == "protein_id":
                    protein_id = value
            attributes_human_dict[protein_id] = gene_id

    with open(input_file) as f:
        for line in f:
            bid, hid, evalu,pident = line.strip().split("\t")
            if bid in attributes_bat_dict:
                bat_gene_id = attributes_bat_dict[bid]
                human_gene_id = attributes_human_dict[hid]
                print(bid, bat_gene_id, hid, human_gene_id, evalu, pident)
            else:
                print(bid, "not found", hid, "NA", evalu, pident)

import sys

input_file = sys.argv[1]
bat_gtf_file = sys.argv[2]
human_gtf_file = sys.argv[3]
extract_genenames(input_file, bat_gtf_file, human_gtf_file)
