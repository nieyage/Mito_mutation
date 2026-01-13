import os
import pandas as pd
import shutil
import sys
import argparse
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import multiprocessing
from functools import partial

def process_single_file(filename, input_dir, output_dir, barcode_to_sample):
    """
    处理单个文件的函数（总是移动文件）
    
    Args:
        filename: 文件名
        input_dir: 输入目录
        output_dir: 输出目录
        barcode_to_sample: barcode到sample的映射字典
    
    Returns:
        tuple: (状态, 文件名, 目标目录, 错误信息)
    """
    try:
        # 提取barcode
        if filename.endswith('.counts'):
            barcode = filename[:-7]  # 移除 .counts
        elif filename.endswith('.snv'):
            barcode = filename[:-4]  # 移除 .snv
        else:
            return ('skip', filename, None, '非目标文件类型')
        
        # 查找对应的sample
        if barcode in barcode_to_sample:
            sample_name = barcode_to_sample[barcode]
            
            # 创建有效的目录名
            donor_dir = sample_name.replace('/', '_').replace('\\', '_')
            donor_path = os.path.join(output_dir, donor_dir)
            
            # 创建目标文件夹（如果不存在）
            if not os.path.exists(donor_path):
                os.makedirs(donor_path, exist_ok=True)
            
            # 移动文件
            src_path = os.path.join(input_dir, filename)
            dst_path = os.path.join(donor_path, filename)
            
            # 总是移动文件（而不是复制）
            shutil.move(src_path, dst_path)
            return ('moved', filename, donor_dir, None)
        else:
            return ('unmatched', filename, None, None)
            
    except Exception as e:
        return ('error', filename, None, str(e))

def organize_files_by_donor_parallel(input_dir, annotation_file, output_dir=None, 
                                    max_workers=None, method='thread'):
    """
    根据CSV文件中的barcode-sample映射，将文件整理到对应的Donor文件夹（并行版本）
    
    Args:
        input_dir: 包含输入文件的目录
        annotation_file: CSV注解文件路径
        output_dir: 输出目录（可选，默认为输入目录）
        max_workers: 最大并行工作线程/进程数
        method: 并行方法，'thread'（线程）或'process'（进程）
    """
    # 如果没有指定输出目录，则使用输入目录
    if output_dir is None:
        output_dir = input_dir
    else:
        # 确保输出目录存在
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            print(f"创建输出目录: {output_dir}")
    
    # 读取CSV文件
    try:
        df = pd.read_csv(annotation_file)
        print(f"成功读取注解文件: {annotation_file}")
        print(f"包含 {len(df)} 条barcode记录")
    except FileNotFoundError:
        print(f"错误：找不到注解文件 {annotation_file}")
        return
    except Exception as e:
        print(f"读取CSV文件时出错: {e}")
        return
    
    # 创建字典映射：barcode -> sample
    barcode_to_sample = dict(zip(df['barcode'], df['sample']))
    
    # 获取输入目录下所有的 .counts 和 .snv 文件
    try:
        files = os.listdir(input_dir)
    except FileNotFoundError:
        print(f"错误：找不到输入目录 {input_dir}")
        return
    
    count_files = [f for f in files if f.endswith('.counts')]
    snv_files = [f for f in files if f.endswith('.snv')]
    
    # 合并所有文件
    all_files = count_files + snv_files
    
    print(f"在目录 {input_dir} 中找到:")
    print(f"  .counts 文件: {len(count_files)} 个")
    print(f"  .snv 文件: {len(snv_files)} 个")
    print(f"  总共: {len(all_files)} 个文件")
    print(f"输出目录: {output_dir}")
    print(f"操作模式: 移动文件（总是移动）")
    
    # 设置并行工作数
    if max_workers is None:
        max_workers = min(32, (multiprocessing.cpu_count() or 1) * 4)
    
    print(f"并行处理: 使用{method}模式，{max_workers}个工作线程")
    print("开始处理文件...")
    
    # 统计信息
    stats = {
        'total': len(all_files),
        'moved': 0,
        'unmatched': 0,
        'skipped': 0,
        'errors': 0
    }
    
    # 创建处理结果记录
    unmatched_files = []
    error_files = []
    
    # 创建固定参数的partial函数
    process_func = partial(process_single_file, 
                          input_dir=input_dir,
                          output_dir=output_dir,
                          barcode_to_sample=barcode_to_sample)
    
    # 根据选择的方法执行并行处理
    if method == 'thread':
        executor_class = ThreadPoolExecutor
    elif method == 'process':
        executor_class = ProcessPoolExecutor
    else:
        print(f"错误：不支持的并行方法 '{method}'，使用线程模式")
        executor_class = ThreadPoolExecutor
    
    with executor_class(max_workers=max_workers) as executor:
        # 提交所有任务
        future_to_file = {executor.submit(process_func, filename): filename 
                         for filename in all_files}
        
        # 处理完成的任务
        completed = 0
        for future in as_completed(future_to_file):
            completed += 1
            filename = future_to_file[future]
            
            try:
                status, fname, target_dir, error = future.result()
                
                # 更新统计信息
                if status == 'moved':
                    stats['moved'] += 1
                    if completed % 10 == 0 or completed == stats['total']:
                        print(f"进度: {completed}/{stats['total']} - 已移动 {filename} -> {target_dir}/")
                elif status == 'unmatched':
                    stats['unmatched'] += 1
                    unmatched_files.append(filename)
                elif status == 'skip':
                    stats['skipped'] += 1
                elif status == 'error':
                    stats['errors'] += 1
                    error_files.append((filename, error))
                    print(f"错误处理文件 {filename}: {error}")
                    
            except Exception as e:
                stats['errors'] += 1
                error_files.append((filename, str(e)))
                print(f"处理文件 {filename} 时发生异常: {e}")
    
    # 打印汇总信息
    print("\n" + "="*60)
    print("文件整理完成!")
    print(f"输入目录: {input_dir}")
    print(f"输出目录: {output_dir}")
    print(f"并行模式: {method} ({max_workers}个工作者)")
    print("\n统计信息:")
    print(f"  总文件数: {stats['total']}")
    print(f"  成功移动: {stats['moved']}")
    print(f"  未匹配: {stats['unmatched']}")
    print(f"  跳过: {stats['skipped']}")
    print(f"  错误: {stats['errors']}")
    
    if unmatched_files:
        print(f"\n未匹配的文件 ({len(unmatched_files)} 个):")
        # 只显示前20个未匹配文件
        for i, file in enumerate(unmatched_files[:20]):
            print(f"  {file}")
        if len(unmatched_files) > 20:
            print(f"  ... 还有 {len(unmatched_files) - 20} 个文件未显示")
    
    if error_files:
        print(f"\n处理错误的文件 ({len(error_files)} 个):")
        for i, (file, error) in enumerate(error_files[:10]):
            print(f"  {file}: {error}")
        if len(error_files) > 10:
            print(f"  ... 还有 {len(error_files) - 10} 个错误未显示")
    
    # 检查创建的文件夹
    print("\n创建的Donor文件夹:")
    donor_dirs = []
    for item in os.listdir(output_dir):
        item_path = os.path.join(output_dir, item)
        if os.path.isdir(item_path) and ('Donor' in item or item.startswith('Donor_')):
            donor_files = os.listdir(item_path)
            donor_dirs.append((item, len(donor_files)))
    
    if donor_dirs:
        for donor_name, file_count in donor_dirs:
            print(f"  {donor_name}: {file_count} 个文件")
    else:
        print("  未创建任何Donor文件夹")

def main():
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(
        description='根据CSV注解文件将文件整理到对应的Donor文件夹（并行版本，总是移动文件）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
使用示例:
  # 基本用法（在当前目录处理，结果也在当前目录）
  python organize_files_move.py -i . -a annotation.csv
  
  # 从输入目录移动到输出目录
  python organize_files_move.py -i /input -a annotation.csv -o /output
  
  # 指定并行处理
  python organize_files_move.py -i /input -a annotation.csv -o /output -w 8 -m process
  
  # 使用线程模式，自动确定工作线程数
  python organize_files_move.py -i /input -a annotation.csv -w 0 -m thread

重要说明:
  - 此脚本总是移动文件，而不是复制
  - 如果输入和输出目录相同，文件会在当前目录内移动
  - 如果输入和输出目录不同，文件会从输入目录移动到输出目录
  - 操作不可逆，请确保已备份重要数据

性能建议:
  - I/O密集型任务（大量小文件）：推荐使用线程模式 (-m thread)
  - CPU密集型任务：推荐使用进程模式 (-m process)
  - 默认线程数: CPU核心数 × 4
  - 对于大量文件，可以适当增加工作线程数（如16-32）
        '''
    )
    
    parser.add_argument('-i', '--input-dir', 
                       required=True,
                       help='包含输入文件的目录路径')
    
    parser.add_argument('-a', '--annotation-file', 
                       required=True,
                       help='CSV注解文件路径')
    
    parser.add_argument('-o', '--output-dir', 
                       default=None,
                       help='输出目录路径（可选，默认为输入目录，文件将在同一目录内移动）')
    
    parser.add_argument('-w', '--workers', 
                       type=int,
                       default=None,
                       help='并行工作线程/进程数（可选，默认自动计算）')
    
    parser.add_argument('-m', '--method', 
                       choices=['thread', 'process'],
                       default='thread',
                       help='并行处理方法：thread(线程)或process(进程，默认: thread)')
    
    parser.add_argument('--no-parallel', 
                       action='store_true',
                       help='禁用并行处理，使用单线程')
    
    # 添加确认选项，避免误操作
    parser.add_argument('-y', '--yes',
                       action='store_true',
                       help='跳过确认提示，直接执行移动操作')
    
    # 解析参数
    args = parser.parse_args()
    
    # 检查输入目录是否存在
    if not os.path.exists(args.input_dir):
        print(f"错误：输入目录不存在 - {args.input_dir}")
        sys.exit(1)
    
    # 检查注解文件是否存在
    if not os.path.exists(args.annotation_file):
        print(f"错误：注解文件不存在 - {args.annotation_file}")
        sys.exit(1)
    
    # 确认操作（如果不使用-y参数）
    if not args.yes:
        print(f"\n警告：此操作将移动文件，不可撤销！")
        print(f"输入目录: {args.input_dir}")
        print(f"输出目录: {args.output_dir if args.output_dir else args.input_dir}")
        print(f"注解文件: {args.annotation_file}")
        
        response = input("\n是否继续？(yes/no): ")
        if response.lower() not in ['yes', 'y']:
            print("操作已取消")
            sys.exit(0)
    
    # 如果禁用并行，使用单线程版本
    if args.no_parallel:
        print("警告：并行处理已禁用，使用单线程模式")
        args.workers = 1
    
    # 执行整理操作
    organize_files_by_donor_parallel(
        args.input_dir, 
        args.annotation_file, 
        args.output_dir,
        args.workers,
        args.method
    )

if __name__ == "__main__":
    main()