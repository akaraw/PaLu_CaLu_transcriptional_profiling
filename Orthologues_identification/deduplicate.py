import sys

def deduplicate(data):
    deduplicated = []
    seen = set()
    for line in data:
        if line[0] not in seen:
            seen.add(line[0])
            deduplicated.append(line)
        else:
            for d in deduplicated:
                if d[0] == line[0]:
                    if float(line[3]) > float(d[3]):
                        deduplicated.remove(d)
                        deduplicated.append(line)
    return deduplicated

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit()

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    with open(input_file, 'r') as f:
        data = [line.strip().split() for line in f.readlines()]

    deduplicated = deduplicate(data)

    with open(output_file, 'w') as f:
        for line in deduplicated:
            f.write("\t".join(line) + "\n")
