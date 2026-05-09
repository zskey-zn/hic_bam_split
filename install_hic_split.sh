#!/bin/bash
# Hi-C BAM文件C扩展一键安装脚本
module load python/3.7.10 samtools/
set -e

echo "=== Hi-C BAM文件C扩展极速拆分工具 ==="
echo "版本: 1.0.0"
echo "作者: zhengshang-zn@qq.com zhengshang@frasergen.com"
echo

# 检查依赖
echo "检查系统依赖..."
for cmd in gcc python3 samtools; do
    if ! command -v $cmd &> /dev/null; then
        echo "错误: $cmd 未安装"
        
        if [ "$cmd" == "samtools" ]; then
            echo "请安装samtools:"
            echo "  conda install -c bioconda samtools"
            echo "或"
            echo "  apt-get install samtools  # Ubuntu/Debian"
            echo "  yum install samtools      # CentOS/RHEL"
        fi
        
        if [ "$cmd" == "gcc" ]; then
            echo "请安装gcc:"
            echo "  apt-get install build-essential  # Ubuntu/Debian"
            echo "  yum groupinstall 'Development Tools'  # CentOS/RHEL"
            echo "  xcode-select --install  # macOS"
        fi
        
        exit 1
    fi
    echo "  ✓ $cmd"
done

echo
echo "所有依赖检查通过！"
echo

# 创建C源文件
echo "创建C源文件..."
cat > hic_split.c << 'EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <pthread.h>
#include <ctype.h>
#include <fcntl.h>
#include <errno.h>
#include <zlib.h>

// 编译命令：
// gcc -O3 -shared -fPIC -o libhic_split.so hic_split.c -lpthread -lz
// gcc -O3 -shared -fPIC -o libhic_split.dylib hic_split.c -lpthread -lz

// 定义常量
#define MAX_LINE_LEN 65536
#define MAX_READ_ID_LEN 1024
#define HASH_SPACE 1048576  // 哈希空间大小
#define BUFFER_SIZE (64 * 1024 * 1024)  // 64MB缓冲区
#define MAX_OUTPUT_FILES 10000  // 最大输出文件数

// 获取当前时间（毫秒）
long long get_time_ms() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (long long)tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

// 简单哈希函数（稳定，相同输入总是得到相同输出）
unsigned long hash_read_name(const char* str) {
    unsigned long hash = 5381;
    int c;
    
    while ((c = *str++)) {
        // 忽略read名称中的/1、/2后缀，确保配对reads有相同哈希
        if (c == '/' && (*str == '1' || *str == '2')) {
            break;
        }
        hash = ((hash << 5) + hash) + c; // hash * 33 + c
    }
    
    return hash % HASH_SPACE;
}

// 提取头部信息
char* extract_header_from_bam(const char* bam_file, size_t* header_len) {
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "samtools view -H '%s' 2>/dev/null", bam_file);
    
    FILE* fp = popen(cmd, "r");
    if (!fp) {
        fprintf(stderr, "无法执行命令: %s\n", cmd);
        return NULL;
    }
    
    // 分配内存存储头部
    char* header = malloc(10 * 1024 * 1024);  // 10MB缓冲区
    if (!header) {
        pclose(fp);
        return NULL;
    }
    
    size_t pos = 0;
    char buffer[4096];
    
    while (fgets(buffer, sizeof(buffer), fp)) {
        size_t len = strlen(buffer);
        if (pos + len >= 10 * 1024 * 1024) {
            break;  // 头部太大
        }
        memcpy(header + pos, buffer, len);
        pos += len;
    }
    
    pclose(fp);
    *header_len = pos;
    header[pos] = '\0';
    
    return header;
}

// 主拆分函数 - 哈希分片
int split_bam_by_read_hash(const char* input_file, int num_chunks, const char* output_prefix, int verbose) {
    if (verbose) {
        printf("开始处理Hi-C BAM文件: %s\n", input_file);
        printf("拆分成: %d 个文件\n", num_chunks);
        printf("输出前缀: %s\n", output_prefix);
        printf("时间: %lld ms\n", get_time_ms());
    }
    
    // 检查参数
    if (num_chunks <= 0 || num_chunks > MAX_OUTPUT_FILES) {
        fprintf(stderr, "错误: 无效的块数 %d (必须在1-%d之间)\n", num_chunks, MAX_OUTPUT_FILES);
        return -1;
    }
    
    // 提取头部
    if (verbose) printf("提取头部信息...\n");
    size_t header_len = 0;
    char* header = extract_header_from_bam(input_file, &header_len);
    if (!header) {
        fprintf(stderr, "错误: 无法提取头部信息\n");
        return -1;
    }
    
    if (verbose) printf("头部长度: %zu 字节\n", header_len);
    
    // 为每个输出文件写入头部
    if (verbose) printf("创建输出文件...\n");
    FILE* output_files[MAX_OUTPUT_FILES];
    
    for (int i = 0; i < num_chunks; i++) {
        char filename[256];
        snprintf(filename, sizeof(filename), "%s_%04d.sam", output_prefix, i + 1);
        
        output_files[i] = fopen(filename, "w");
        if (!output_files[i]) {
            fprintf(stderr, "错误: 无法创建输出文件 %s\n", filename);
            free(header);
            for (int j = 0; j < i; j++) fclose(output_files[j]);
            return -1;
        }
        
        // 写入头部
        if (header_len > 0) {
            fwrite(header, 1, header_len, output_files[i]);
        }
        
        if (verbose && (i < 10 || i % 100 == 0)) {
            printf("  创建文件: %s\n", filename);
        }
    }
    
    // 构建samtools命令
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "samtools view '%s' 2>/dev/null", input_file);
    
    FILE* fp = popen(cmd, "r");
    if (!fp) {
        fprintf(stderr, "错误: 无法执行samtools命令\n");
        free(header);
        for (int i = 0; i < num_chunks; i++) fclose(output_files[i]);
        return -1;
    }
    
    // 处理数据
    if (verbose) printf("开始处理数据...\n");
    
    long long start_time = get_time_ms();
    long total_records = 0;
    char line[MAX_LINE_LEN];
    
    while (fgets(line, sizeof(line), fp)) {
        // 跳过头部行
        if (line[0] == '@') continue;
        
        // 提取read名称
        char* tab_pos = strchr(line, '\t');
        if (!tab_pos) continue;
        
        // 计算哈希值
        int name_len = tab_pos - line;
        char read_name[MAX_READ_ID_LEN];
        if (name_len >= MAX_READ_ID_LEN) name_len = MAX_READ_ID_LEN - 1;
        
        strncpy(read_name, line, name_len);
        read_name[name_len] = '\0';
        
        unsigned long hash = hash_read_name(read_name);
        int file_idx = hash % num_chunks;
        
        // 写入对应文件
        if (output_files[file_idx]) {
            fputs(line, output_files[file_idx]);
        }
        
        total_records++;
        
        // 进度显示
        if (verbose && total_records % 1000000 == 0) {
            long long current_time = get_time_ms();
            long long elapsed = current_time - start_time;
            double speed = (elapsed > 0) ? (double)total_records * 1000 / elapsed : 0;
            printf("已处理: %ld 行, 速度: %.0f 行/秒\n", total_records, speed);
            
            // 定期刷新所有文件
            for (int i = 0; i < num_chunks; i++) {
                fflush(output_files[i]);
            }
        }
    }
    
    // 关闭管道
    int status = pclose(fp);
    if (status != 0) {
        fprintf(stderr, "警告: samtools命令返回非零状态: %d\n", status);
    }
    
    // 关闭所有输出文件
    for (int i = 0; i < num_chunks; i++) {
        fclose(output_files[i]);
    }
    
    // 清理
    free(header);
    
    long long end_time = get_time_ms();
    long long elapsed = end_time - start_time;
    double speed = (elapsed > 0) ? (double)total_records * 1000 / elapsed : 0;
    
    if (verbose) {
        printf("\n处理完成！\n");
        printf("总记录数: %ld\n", total_records);
        printf("创建文件数: %d\n", num_chunks);
        printf("总耗时: %.2f 秒\n", elapsed / 1000.0);
        printf("平均速度: %.0f 行/秒\n", speed);
    }
    
    return 0;
}

// Python调用的接口函数
int split_hic_bam(const char* input_file, int num_chunks, const char* output_prefix, int verbose) {
    return split_bam_by_read_hash(input_file, num_chunks, output_prefix, verbose);
}

// 主函数（用于独立测试）
int main(int argc, char** argv) {
    if (argc < 4) {
        printf("用法: %s <输入BAM文件> <块数> <输出前缀> [详细]\n", argv[0]);
        printf("  详细: 0=关闭, 1=开启\n");
        printf("示例: %s input.bam 10 output 1\n", argv[0]);
        return 1;
    }
    
    const char* input_file = argv[1];
    int num_chunks = atoi(argv[2]);
    const char* output_prefix = argv[3];
    int verbose = (argc > 4) ? atoi(argv[4]) : 1;
    
    return split_hic_bam(input_file, num_chunks, output_prefix, verbose);
}
EOF

echo "C源文件创建成功: hic_split.c"

# 编译C扩展
echo
echo "编译C扩展..."

if [[ "$OSTYPE" == "darwin"* ]]; then
    COMPILE_CMD="gcc -O3 -shared -fPIC -o libhic_split.dylib hic_split.c -lpthread -lz"
    LIB_NAME="libhic_split.dylib"
else
    COMPILE_CMD="gcc -O3 -shared -fPIC -o libhic_split.so hic_split.c -lpthread -lz"
    LIB_NAME="libhic_split.so"
fi

echo "执行: $COMPILE_CMD"
eval $COMPILE_CMD

if [ $? -eq 0 ]; then
    echo "编译成功: $LIB_NAME"
else
    echo "编译失败！"
    exit 1
fi

# 创建Python包装器
echo
echo "创建Python包装器..."

cat > hic_split.py << 'EOF'
#!/usr/bin/env python3
"""
Hi-C BAM文件C扩展极速拆分工具
"""

import os
import sys
import ctypes
import subprocess
import argparse
from datetime import datetime

def split_hic_bam_c(input_file, num_chunks, output_prefix, verbose=True):
    """使用C扩展拆分BAM文件"""
    
    # 加载库
    if sys.platform == 'darwin':
        lib_name = './libhic_split.dylib'
    else:
        lib_name = './libhic_split.so'
    
    if not os.path.exists(lib_name):
        print(f"错误: C扩展库不存在: {lib_name}")
        return False
    
    try:
        lib = ctypes.CDLL(lib_name)
    except Exception as e:
        print(f"加载C扩展失败: {e}")
        return False
    
    # 设置函数参数类型
    lib.split_hic_bam.argtypes = [
        ctypes.c_char_p,  # input_file
        ctypes.c_int,     # num_chunks
        ctypes.c_char_p,  # output_prefix
        ctypes.c_int      # verbose
    ]
    lib.split_hic_bam.restype = ctypes.c_int
    
    # 调用C函数
    verbose_code = 1 if verbose else 0
    
    print(f"开始处理: {input_file}")
    print(f"拆分成: {num_chunks} 个文件")
    print(f"输出前缀: {output_prefix}")
    print(f"时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    start_time = datetime.now()
    
    result = lib.split_hic_bam(
        input_file.encode('utf-8'),
        num_chunks,
        output_prefix.encode('utf-8'),
        verbose_code
    )
    
    elapsed = (datetime.now() - start_time).total_seconds()
    
    if result == 0:
        print(f"\nC扩展拆分成功！耗时: {elapsed:.1f}秒")
        return True
    else:
        print(f"\nC扩展拆分失败，错误码: {result}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Hi-C BAM文件C扩展极速拆分')
    parser.add_argument('input_file', help='输入BAM文件')
    parser.add_argument('--chunks', type=int, required=True, help='拆分数')
    parser.add_argument('--prefix', default='split', help='输出前缀')
    parser.add_argument('--no-verbose', action='store_true', help='关闭详细输出')
    
    args = parser.parse_args()
    
    # 检查文件
    if not os.path.exists(args.input_file):
        print(f"错误: 文件不存在: {args.input_file}")
        sys.exit(1)
    
    # 检查samtools
    if not shutil.which('samtools'):
        print("错误: samtools未安装")
        sys.exit(1)
    
    # 执行拆分
    success = split_hic_bam_c(
        args.input_file,
        args.chunks,
        args.prefix,
        verbose=not args.no_verbose
    )
    
    if success:
        print(f"\n处理完成！输出文件: {args.prefix}_*.sam")
    else:
        print("\n处理失败！")
        sys.exit(1)

if __name__ == "__main__":
    # 添加shutil到局部
    import shutil
    main()
EOF

# 创建bash包装器
echo
echo "创建bash包装器..."

cat > hic_split_cli << 'EOF'
#!/bin/bash
# Hi-C BAM文件C扩展命令行工具

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIB_PATH="$SCRIPT_DIR/libhic_split.so"

if [[ "$OSTYPE" == "darwin"* ]]; then
    LIB_PATH="$SCRIPT_DIR/libhic_split.dylib"
fi

if [[ ! -f "$LIB_PATH" ]]; then
    echo "错误: C扩展库不存在: $LIB_PATH"
    echo "请先运行 ./install_hic_split.sh"
    exit 1
fi

if [[ $# -lt 3 ]]; then
    echo "用法: $0 <输入BAM文件> <块数> <输出前缀> [详细]"
    echo "  详细: 0=关闭, 1=开启 (默认: 1)"
    echo "示例: $0 input.bam 10 output 1"
    exit 1
fi

INPUT="$1"
CHUNKS="$2"
PREFIX="$3"
VERBOSE="${4:1}"

# 检查文件
if [[ ! -f "$INPUT" ]]; then
    echo "错误: 文件不存在: $INPUT"
    exit 1
fi

# 检查samtools
if ! command -v samtools &> /dev/null; then
    echo "错误: samtools未安装"
    exit 1
fi

# 运行C程序
echo "开始处理: $INPUT"
echo "时间: $(date)"

"$SCRIPT_DIR/hic_split_exec" "$INPUT" "$CHUNKS" "$PREFIX" "$VERBOSE"

if [[ $? -eq 0 ]]; then
    echo "处理完成！"
    echo "输出文件: ${PREFIX}_*.sam"
else
    echo "处理失败！"
    exit 1
fi
EOF

# 编译独立的可执行文件
echo
echo "编译独立可执行文件..."

if [[ "$OSTYPE" == "darwin"* ]]; then
    gcc -O3 -o hic_split_exec hic_split.c -lpthread -lz
else
    gcc -O3 -o hic_split_exec hic_split.c -lpthread -lz
fi

chmod +x hic_split_exec
chmod +x hic_split.py
chmod +x hic_split_cli

echo
echo "=================================================="
echo "安装完成！"
echo
echo "可用工具："
echo "  1. C扩展库: $LIB_NAME"
echo "  2. Python包装器: ./hic_split.py"
echo "  3. 命令行工具: ./hic_split_cli"
echo "  4. 独立可执行文件: ./hic_split_exec"
echo
echo "使用方法："
echo "  # Python版本"
echo "  ./hic_split.py input.bam --chunks 100 --prefix split"
echo
echo "  # 命令行版本"
echo "  ./hic_split_cli input.bam 100 split"
echo
echo "  # 直接使用C程序"
echo "  ./hic_split_exec input.bam 100 split 1"
echo
echo "  # 对于14TB Hi-C文件建议"
echo "  ./hic_split_cli 14TB_hic.bam 1000 split_parts 1"
echo
echo "注意："
echo "  - 需要安装samtools来处理BAM文件"
echo "  - 保证同一对reads的所有记录在同一文件中"
echo "  - 使用哈希分片算法，速度快"
echo "=================================================="
