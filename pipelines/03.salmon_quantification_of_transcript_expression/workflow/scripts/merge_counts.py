
import argparse
from operator import itemgetter

# Open count files as dictionaries (key=geneid, value=count)
def open_count_files(input_file, delimiter="\t"):
    with open(input_file, "r") as f:
        counts = {}
        lines = f.readlines()

        for line in lines:
            key, values = line.strip().split(delimiter, 1)
            counts[key] = values.split('\t')[4] ## NumReads is the fifth and last column of quant.sf (Salmon output)

    return counts

def merge_counts(input_list):
    dictionaries = []

    for file_path in open(input_list, "r"):
        counts = open_count_files(file_path.strip('\n'))
        dictionaries.append(counts)

    if len(dictionaries) > 1:
        common_keys = set(dictionaries[0].keys())

        for counts in dictionaries[1:]:
            common_keys &= set(counts.keys())

        combined_values = {}
        for key in common_keys:
            values_list = [counts[key] for counts in dictionaries]
            combined_values[key] = '\t'.join(map(str, values_list))

        # Print common keys and values
        for key, value in sorted(combined_values.items(), key=itemgetter(0)):
            yield f"{key}\t{value}"

        # Print keys and values for dictionaries that have unique keys
        for counts in dictionaries:
            unique_keys = set(counts.keys()) - common_keys
            for key in sorted(unique_keys):
                yield f"{key}\t{counts[key]}"

    else:
        print("Not enough files for comparison.")

        
## Implementation ##
def parse_arguments():
    parser = argparse.ArgumentParser(description='Merge counts from multiple samples based on geneid.')
    parser.add_argument('-c','--counts', required=True, type=str, help='List of counts files you want to merge')
    parser.add_argument('-o','--out', required=False , type=str, help='Output file')
    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
        print("Error: No arguments provided.")

    return parser.parse_args()


def main():
    
    parser = argparse.ArgumentParser(description='Merge count files')
    args = parse_arguments()

    with open(args.out if args.out else "merged_counts.tsv", "w") as out:
        for c in (merge_counts(args.counts)):
            print (c, file=out)

    
if __name__ == "__main__":
    main()
    
