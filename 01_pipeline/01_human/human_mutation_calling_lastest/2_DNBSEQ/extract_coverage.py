#!/usr/bin/env python3
import os
import sys
import gzip
import glob
import argparse
from multiprocessing import Pool, cpu_count
import time
from tqdm import tqdm

def process_file(file_info):
    """处理单个文件"""
    f, barcode = file_info
    try:
        results = []
        with open(f, 'r') as file:
            for line in file:
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    chrom = parts[0]
                    if chrom == 'chrM' or chrom.isdigit() or chrom in ['X', 'Y']:
                        results.append(f"{barcode}\t{parts[1]}\t{parts[3]}\n")
        return results
    except Exception as e:
        print(f"Error processing {f}: {e}")
        return []

def parallel_method(input_dir, output_file, num_processes=8, batch_size=500):
    """并行处理方法"""
    # 查找所有文件
    files = glob.glob(os.path.join(input_dir, "*.counts"))
    
    if not files:
        print(f"No .counts files found in {input_dir}")
        sys.exit(1)
    
    print(f"Found {len(files):,} files to process")
    
    # 准备文件信息列表
    file_info_list = [(f, os.path.basename(f).replace('.counts', '')) for f in files]
    
    # 设置进程数
    num_processes = min(cpu_count(), num_processes)
    print(f"Using {num_processes} processes")
    
    start_time = time.time()
    
    # 直接写入最终文件
    with gzip.open(output_file, 'wt', encoding='utf-8') as f_out:
        # 写入表头
        f_out.write("Barcode\tPosition\tCount\n")
        
        # 分批处理以避免内存问题
        total_batches = (len(file_info_list) + batch_size - 1) // batch_size
        
        for batch_idx in range(0, len(file_info_list), batch_size):
            batch = file_info_list[batch_idx:batch_idx + batch_size]
            batch_num = batch_idx // batch_size + 1
            
            print(f"\nProcessing batch {batch_num}/{total_batches} ({len(batch)} files)")
            batch_start = time.time()
            
            # 使用进程池处理当前批次
            with Pool(processes=num_processes) as pool:
                # 使用imap_unordered，chunksize控制任务分配粒度
                results = list(tqdm(
                    pool.imap_unordered(process_file, batch, chunksize=10),
                    total=len(batch),
                    desc=f"Batch {batch_num}"
                ))
            
            # 写入当前批次的结果
            print(f"Writing batch {batch_num} results...")
            write_start = time.time()
            
            for file_results in results:
                for line in file_results:
                    f_out.write(line)
            
            write_time = time.time() - write_start
            batch_time = time.time() - batch_start
            
            # 进度统计
            processed_so_far = min(batch_idx + batch_size, len(file_info_list))
            elapsed_total = time.time() - start_time
            speed = processed_so_far / elapsed_total if elapsed_total > 0 else 0
            
            print(f"Batch {batch_num} completed in {batch_time:.1f}s (write: {write_time:.1f}s)")
            print(f"Progress: {processed_so_far:,}/{len(file_info_list):,} files ({processed_so_far/len(file_info_list)*100:.1f}%)")
            print(f"Average speed: {speed:.1f} files/sec")
    
    return files, start_time

def simple_method(input_dir, output_file):
    """简单的串行方法"""
    files = glob.glob(os.path.join(input_dir, "*.counts"))
    
    print(f"Processing {len(files):,} files (simple method)...")
    
    start_time = time.time()
    
    with gzip.open(output_file, 'wt', encoding='utf-8') as f_out:
        f_out.write("Barcode\tPosition\tCount\n")
        
        processed = 0
        for f in tqdm(files, desc="Processing"):
            barcode = os.path.basename(f).replace('.counts', '')
            
            with open(f, 'r') as f_in:
                for line in f_in:
                    if not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) >= 4:
                        chrom = parts[0]
                        if chrom == 'chrM' or chrom.isdigit() or chrom in ['X', 'Y']:
                            f_out.write(f"{barcode}\t{parts[1]}\t{parts[3]}\n")
            
            processed += 1
            
            # 每1000个文件报告一次进度
            if processed % 1000 == 0:
                elapsed = time.time() - start_time
                speed = processed / elapsed
                remaining = (len(files) - processed) / speed if speed > 0 else 0
                print(f"  {processed:,}/{len(files):,} files - {speed:.1f} files/sec - ETA: {remaining/60:.1f} min")
    
    return files, start_time

def print_summary(files, start_time, output_file, method_name):
    """打印汇总信息"""
    total_time = time.time() - start_time
    
    print(f"\n{'='*60}")
    print(f"✓ {method_name.upper()} COMPLETE!")
    print(f"{'='*60}")
    print(f"Total files processed: {len(files):,}")
    print(f"Total time: {total_time:.1f} seconds")
    print(f"Average speed: {len(files)/total_time:.1f} files/second")
    print(f"Output file: {output_file}")
    
    # 检查输出文件大小
    if os.path.exists(output_file):
        size_bytes = os.path.getsize(output_file)
        size_mb = size_bytes / (1024**2)
        size_gb = size_bytes / (1024**3)
        print(f"Output file size: {size_mb:.2f} MB ({size_gb:.3f} GB)")
        
        # 显示文件前几行
        print(f"\nFirst few lines of output:")
        try:
            with gzip.open(output_file, 'rt', encoding='utf-8') as f:
                lines_read = 0
                for line in f:
                    print(f"  {line.strip()}")
                    lines_read += 1
                    if lines_read >= 3:
                        break
        except Exception as e:
            print(f"  Could not read output file: {e}")

def main():
    parser = argparse.ArgumentParser(
        description='Extract coverage information from .counts files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i /path/to/counts -o output.tsv.gz
  %(prog)s -i /path/to/counts -o output.tsv.gz --simple
  %(prog)s -i /path/to/counts -o output.tsv.gz -p 16 -b 1000
  %(prog)s --input /path/to/counts --output output.tsv.gz --simple --verbose
        """
    )
    
    parser.add_argument('-i', '--input', type=str, required=True,
                       help='Input directory containing .counts files')
    
    parser.add_argument('-o', '--output', type=str, required=True,
                       help='Output gzipped TSV file (e.g., output.tsv.gz)')
    
    parser.add_argument('-p', '--processes', type=int, default=8,
                       help='Number of parallel processes (default: 8)')
    
    parser.add_argument('-b', '--batch-size', type=int, default=500,
                       help='Batch size for parallel processing (default: 500)')
    
    parser.add_argument('--simple', action='store_true',
                       help='Use simple serial method instead of parallel')
    
    parser.add_argument('--chunksize', type=int, default=10,
                       help='Chunk size for parallel processing (default: 10)')
    
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')
    
    parser.add_argument('--check', action='store_true',
                       help='Check a few files before processing')
    
    args = parser.parse_args()
    
    # 验证输入目录
    if not os.path.isdir(args.input):
        print(f"Error: Input directory does not exist: {args.input}")
        sys.exit(1)
    
    # 验证输出文件目录
    output_dir = os.path.dirname(args.output) or '.'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    # 检查文件（可选）
    if args.check:
        files = glob.glob(os.path.join(args.input, "*.counts"))
        if not files:
            print(f"Error: No .counts files found in {args.input}")
            sys.exit(1)
        
        print(f"Found {len(files)} .counts files")
        print("\nChecking first 3 files:")
        for i, f in enumerate(files[:3]):
            print(f"\nFile {i+1}: {os.path.basename(f)}")
            try:
                with open(f, 'r') as test_file:
                    lines = test_file.readlines()[:3]
                    for j, line in enumerate(lines):
                        print(f"  Line {j+1}: {line.strip()}")
            except Exception as e:
                print(f"  Error reading file: {e}")
        
        answer = input("\nContinue processing? (y/n): ")
        if answer.lower() != 'y':
            print("Processing cancelled.")
            sys.exit(0)
    
    try:
        if args.simple:
            print("Using simple serial method...")
            files, start_time = simple_method(args.input, args.output)
            print_summary(files, start_time, args.output, "SERIAL PROCESSING")
        else:
            print("Using parallel method...")
            files, start_time = parallel_method(
                args.input, 
                args.output, 
                args.processes, 
                args.batch_size
            )
            print_summary(files, start_time, args.output, "PARALLEL PROCESSING")
    
    except KeyboardInterrupt:
        print("\n\nProcessing interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nError during processing: {e}")
        print("\nTrying simple method instead...")
        try:
            files, start_time = simple_method(args.input, args.output)
            print_summary(files, start_time, args.output, "SERIAL FALLBACK PROCESSING")
        except Exception as e2:
            print(f"Simple method also failed: {e2}")
            sys.exit(1)

if __name__ == "__main__":
    main()