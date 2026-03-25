#!/usr/bin/env python
import pysam
import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Split bam into individual bams')
parser.add_argument('-i', '--input', type=str, required=True, help='The name of Sam/Bam files for parsing')
parser.add_argument('-b', '--barcodes', type=str, required=True, help='Barcode file for parsing')
parser.add_argument('-o', '--output', type=str, default='./splitted_shift/', help='Output directory (default: ./splitted_shift/)')
parser.add_argument('--prefix', type=str, default='', help='Prefix for output bam files (optional)')
args = parser.parse_args()

# Input files
bam_file = args.input
barcode_file_path = args.barcodes
out_dir = args.output
prefix = args.prefix

# Create output directory
os.makedirs(out_dir, exist_ok=True)

# Variables to keep record of barcode index
oldCB = ''
count = 0

# Read the valid barcode file
with open(barcode_file_path, "r") as f:
    barcodes = [line.strip() for line in f if line.strip()]
print(f"Loaded {len(barcodes)} barcodes from file")

# Read input BAM file and loop through reads
sam_file = pysam.AlignmentFile(bam_file, "rb")
sam_dict = {}
sam_list = []
no_cb = 0

print(f"Processing BAM file: {bam_file}")

for read in sam_file.fetch(until_eof=True):  # since there is no index for it, until_eof=True is required
    count += 1
    
    # Progress indicator
    if count % 1000000 == 0:
        print(f"Processed {count:,} reads...")
    
    if "CB:Z" not in read.tostring():
        no_cb += 1
    else:
        try:
            CB = read.get_tag('CB')
            if CB == oldCB:
                sam_list.append(read)
            else:
                if oldCB == '':
                    sam_list = [read]
                    oldCB = CB
                else:
                    sam_dict[oldCB] = sam_list
                    if len(sam_dict) % 100 == 0:
                        print(f"Added {oldCB} to dict (total: {len(sam_dict)})")
                    sam_list = [read]
                    oldCB = CB
        except KeyError:
            no_cb += 1

# Add the last record to the dict
if CB and sam_list:
    sam_dict[CB] = sam_list

print(f"\nProcessing complete. Found {len(sam_dict)} unique barcodes")

# Filter valid barcodes
valid_barcodes = [barcode for barcode in barcodes if barcode in sam_dict]
print(f"\n{len(valid_barcodes)} barcodes match between barcode file and BAM file")

# Write individual BAM files
for i, barcode in enumerate(valid_barcodes):
    output_file = os.path.join(out_dir, f"{prefix}{barcode}.bam")
    
    # Progress indicator
    if i % 100 == 0:
        print(f"Writing bam files: {i}/{len(valid_barcodes)}...")
    
    CB_sam_file = pysam.AlignmentFile(output_file, "wb", template=sam_file)
    for read in sam_dict[barcode]:
        CB_sam_file.write(read)
    CB_sam_file.close()
    
    # Optional: Create index for each BAM file
    # pysam.index(output_file)

print(f"\n{'='*60}")
print(f"Summary:")
print(f"{'='*60}")
print(f"Total alignments processed: {count:,}")
print(f"Alignments with no cell barcodes: {no_cb:,}")
print(f"Alignments with valid barcodes: {count - no_cb:,}")
print(f"Unique barcodes found: {len(sam_dict):,}")
print(f"Barcodes written to individual files: {len(valid_barcodes):,}")
print(f"Output directory: {os.path.abspath(out_dir)}")
print(f"{'='*60}")

# Close the main BAM file
sam_file.close()