import argparse

def merge_files(file1_path, file2_path, output_path):
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2, open(output_path, 'w') as outfile:
        for line1, line2 in zip(file1, file2):
            parts1 = line1.strip().split()
            parts2 = line2.strip().split()[1:]  # Skip the first element of line2
            merged_line = parts1 + parts2
            outfile.write(' '.join(merged_line) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Merge two files line by line, skipping first element of second file\'s lines.')
    parser.add_argument('file1', help='Path to the first input file')
    parser.add_argument('file2', help='Path to the second input file')
    parser.add_argument('output', help='Path to the output file')

    args = parser.parse_args()
    merge_files(args.file1, args.file2, args.output)

if __name__ == '__main__':
    main()

