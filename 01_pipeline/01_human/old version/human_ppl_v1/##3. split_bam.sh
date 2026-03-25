# Split bam by barcode

#!/bin/python
import pysam
import os
import subprocess

# input files
unshifted_file = "PBMC_chrM_unmapped_bwa_unshifted_realign_sorted.bam"
shifted_file = "PBMC_chrM_unmapped_bwa_shifted_realign_sorted.bam"
out_dir = "./PBMC_split_bam/"

os.makedirs(out_dir, exist_ok=True)

def get_all_barcodes(bam_file):
    """Get all unique barcodes from BAM file"""
    barcodes = set()
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for read in samfile.fetch(until_eof=True):
            try:
                barcodes.add(read.get_tag('CB'))
            except:
                continue
    return barcodes

# Get all barcodes from both files
print("Collecting barcodes...")
all_barcodes = set()
for bam_file in [unshifted_file, shifted_file]:
    barcodes = get_all_barcodes(bam_file)
    all_barcodes.update(barcodes)
    print(f"Found {len(barcodes)} barcodes in {bam_file}")

print(f"Total unique barcodes: {len(all_barcodes)}")

# Create directories for all barcodes
for barcode in all_barcodes:
    barcode_dir = os.path.join(out_dir, barcode)
    os.makedirs(barcode_dir, exist_ok=True)

def split_bam_by_barcode(input_file, output_suffix):
    """Split BAM file by barcode into pre-created directories"""
    print(f"Splitting {input_file}...")
    
    CB_hold = 'unset'
    lit = []
    
    with pysam.AlignmentFile(input_file, "rb") as samfile:
        for read in samfile.fetch(until_eof=True):
            try:
                CB_itr = read.get_tag('CB')
            except:
                continue
            
            if CB_itr != CB_hold:
                if lit:  # process previous batch
                    barcode_dir = os.path.join(out_dir, CB_hold)
                    output_file = os.path.join(barcode_dir, f"{output_suffix}.bam")
                    with pysam.AlignmentFile(output_file, "wb", template=samfile) as out_bam:
                        for read_item in lit:
                            out_bam.write(read_item)
                
                CB_hold = CB_itr
                lit = []
            
            lit.append(read)
        
        # Process last batch
        if lit:
            barcode_dir = os.path.join(out_dir, CB_hold)
            output_file = os.path.join(barcode_dir, f"{output_suffix}.bam")
            with pysam.AlignmentFile(output_file, "wb", template=samfile) as out_bam:
                for read_item in lit:
                    out_bam.write(read_item)

# Process both files
split_bam_by_barcode(unshifted_file, "unshifted")
split_bam_by_barcode(shifted_file, "shifted")

print("All files processed successfully!")
