import os
def replace_values(htseq_dir, orthotsv_bat):
    # Create a dictionary from the second file, with the 4th column as the key
    # and the 1st column as the value
    replace_dict = {}
    ortholist = orthotsv_bat
    with open(ortholist, "r") as f:
        for line in f:
            if line.startswith("_"):
                continue
            values = line.strip().split()
            replace_dict[values[0]] = values[0]

    folder_path =htseq_dir

    for filename in os.listdir(folder_path):
        if filename.endswith(".tsv"):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, "r") as f1:
                new_filename = filename[:-4] + "ortho.tsv"
                new_file_path = os.path.join(folder_path, new_filename)
                with open(new_file_path, "w") as f2:
                    for line in f1:
                        if line.startswith("_"):
                            line = line.strip()
                            f2.write(line + "\n")
                        else:
                            values = line.strip().split()
                            if values[0] in replace_dict:
                                values[0] = replace_dict[values[0]]
                            else:
                                continue
                            f2.write("\t".join(values) + "\n")

import sys

htseq_dir = sys.argv[1]
ortholist = sys.argv[2]


replace_values(htseq_dir, ortholist)
