#!/usr/bin/env python3

import sys
import numpy as np
import tempfile
import os
import pysam

os.environ['TMPDIR'] = os.getcwd()

def process_genome(file_path):
    hs_genome = {}
    with open(file_path, 'r') as file:
        tag = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                tag = line[1:].split()[0]
                hs_genome[tag] = ""
            else:
                hs_genome[tag] += line
    return hs_genome

def process_target(file_path, hs_genome, temp_file_path):
    with open(file_path, 'r') as file, open(temp_file_path, 'a') as temp_file:
        for line in file:
            line = line.strip()
            #if not line.startswith('chr'):
            #    continue  # Skip lines not starting with 'chr'
            ar = line.split('\t')
            
            # Genomic coordinates
            chrom = ar[1]
            start = int(ar[5]) + 1
            end = int(ar[6]) + 1
            seq_length = end - start
            
            # Form the tag for the output
            tag = f"{ar[0]}\t0\t{chrom}\t{start}\t60\t{seq_length}M\t*\t0\t0"

            # Extract 6mA and 5mC sites
            ref6mA = ar[8]
            ref5mC = ar[10]

            try:
                ref6mA_array = np.fromstring(ref6mA[1:-1], sep=', ')  # Assuming comma-separated values without brackets
                #ref5mC_array = np.fromstring(ref5mC[1:-1], sep=', ')
            except ValueError:
                #print(f"Error parsing arrays for line: {line}")
                continue
            
            # Combine and filter valid modification sites
            refmodsite = ref6mA_array
            #refmodsite = np.concatenate([ref6mA_array, ref5mC_array])
            refmodsite = refmodsite[refmodsite != -1]

            if refmodsite.size > 0:
                refmodml = np.full(refmodsite.shape, 255)  # Create an array filled with 255
                refmodsite = np.sort(refmodsite - int(ar[5]))  # Sort the mod sites and adjust by start
                refmodsitedf = (np.diff(refmodsite)) - 1
                refmodfinal = np.insert(refmodsitedf, 0, refmodsite[0])  # Keep the first element
                refmodfinal = [int(x) for x in refmodfinal]
                mmstr = ','.join(map(str, refmodfinal))
                mlstr = ','.join(map(str, refmodml))
                
                try:
                    sequence = hs_genome[chrom][int(ar[5]):int(ar[6])]
                except KeyError:
                    print(f"Chromosome {chrom} not found in genome data.")
                    continue

                temp_file.write(f"{tag}\t{sequence}\t*\tMM:Z:N+n?,{mmstr}\tML:B:C,{mlstr}\n")
            else:
                try:
                    sequence = hs_genome[chrom][int(ar[5]):int(ar[6])]
                except KeyError:
                    print(f"Chromosome {chrom} not found in genome data.")
                    continue

                temp_file.write(f"{tag}\t{sequence}\t*\n")

def extract_header(bam_file_path, temp_file):
    with pysam.AlignmentFile(bam_file_path, 'rb', check_sq = False) as bam_file:
        for header_line in bam_file.text.splitlines():
            temp_file.write(header_line + '\n')


def merge_and_sort(temp_file_path, output_bam_path):
    # Convert the temporary file to a BAM file
    temp_bam_path = temp_file_path.replace(".sam", ".bam")
    pysam.view("-bS", "-o", temp_bam_path, temp_file_path, catch_stdout = False)

    # Sort the BAM file
    #sorted_bam_path = output_bam_path.replace(".bam", ".sorted.bam")
    pysam.sort("-o", output_bam_path, temp_bam_path)
    pysam.index(output_bam_path)

    # Clean up temporary BAM file
    os.remove(temp_bam_path)
    #return sorted_bam_path


def main():
    if len(sys.argv) < 5:
        print("Usage: script.py <genome_file> <target_file> <aligned bam> <output.bam>")
        sys.exit(1)

    genome_file = sys.argv[1]
    target_file = sys.argv[2]
    alignbam_file = sys.argv[3]
    out_bam = sys.argv[4]
    # Process the genome file
    hs_genome = process_genome(genome_file)

    # Create a temporary file to store processed target data
    with tempfile.NamedTemporaryFile(delete=False, suffix = ".sam") as temp_file:
        temp_file_path = temp_file.name

    with open(temp_file_path, 'w') as temp_file:
        extract_header(alignbam_file, temp_file)

    # Process the target file
    process_target(target_file, hs_genome, temp_file_path)

    #header = extract_header(alignbam_file, temp_file_path)
    # Optionally, print the header if provided
    merge_and_sort(temp_file_path, out_bam)

    print(f"Sorted BAM file created: {out_bam}")
    os.remove(temp_file_path)
    #os.remove(out_bam)

if __name__ == "__main__":
    main()
