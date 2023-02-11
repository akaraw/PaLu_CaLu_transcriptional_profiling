import sys
from Bio import SeqIO

fasta_file = sys.argv[1]

def extract_longest_isoform(fasta_file):
    isoform_dict = {}

    # Parse the fasta file using SeqIO
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        sequence = str(record.seq)

        if "isoform" in header:
            # Extract the isoform name from the header
            line = header.split(" isoform ")[0]
            split_line = line.split()
            isoform_name = " ".join(split_line[1:])
            # If this is the first time we have seen this isoform, add it to the dictionary
            if isoform_name not in isoform_dict:
                isoform_dict[isoform_name] = {"length": len(sequence), "header": header, "sequence": sequence}
            else:
                # If this is not the first time we have seen this isoform, compare the length with the previous isoform
                if len(sequence) > isoform_dict[isoform_name]["length"]:
                    isoform_dict[isoform_name] = {"length": len(sequence), "header": header, "sequence": sequence}
        else:
            isoform_name = header
            isoform_dict[isoform_name] = {"length": len(sequence), "header": header, "sequence": sequence}
    # Create a list of the longest isoforms
    longest_isoforms = [isoform_dict[key] for key in isoform_dict.keys()]

    return longest_isoforms

longest_isoforms = extract_longest_isoform(fasta_file)

# Print the longest isoforms
for isoform in longest_isoforms:
    print(">" + isoform["header"])
    print(isoform["sequence"])
