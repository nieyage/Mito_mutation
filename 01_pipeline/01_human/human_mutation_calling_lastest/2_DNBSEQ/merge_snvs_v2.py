import os
import glob
import argparse

def merge_snv_files(input_dir, output_file):
    # 1. 查找所有 .snv 文件
    search_path = os.path.join(input_dir, "*.snv")
    file_list = sorted(glob.glob(search_path))
    
    if not file_list:
        print(f"Error: No .snv files found in {input_dir}")
        return

    print(f"Found {len(file_list)} SNV files, starting merge...")

    with open(output_file, 'w') as fout:
        # 2. 从第一个文件提取表头并写入 cell_id
        with open(file_list[0], 'r') as f:
            header = f.readline().strip()
            fout.write(f"cell_id\t{header}\n")

        # 3. 逐个合并数据行
        for i, file_path in enumerate(file_list):
            # 获取文件名（不含后缀）作为细胞名
            cell_name = os.path.basename(file_path).replace(".snv", "")
            
            with open(file_path, 'r') as fin:
                # 跳过当前文件的表头行
                next(fin)
                for line in fin:
                    # 去除行末换行符并重新组合，确保 cell_id 后面是 Tab
                    line_content = line.strip()
                    if line_content:
                        fout.write(f"{cell_name}\t{line_content}\n")
            
            if (i + 1) % 500 == 0:
                print(f"Processed {i + 1} files...")

    print(f"Success! Merged file saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge SNV files and add cell_id column")
    parser.add_argument("-i", "--input", required=True, help="Input directory")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    
    args = parser.parse_args()
    merge_snv_files(args.input, args.output)