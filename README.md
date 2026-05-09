# hic_bam_split
Use hash split large Hi-C bam

# Install
bash install_hic_split.sh

# exec file
  1. C扩展库: libhic_split.so
  2. Python包装器: ./hic_split.py
  3. 命令行工具: ./hic_split_cli
  4. 独立可执行文件: ./hic_split_exec

# Usage
  ```bash
  # Python版本
  python3 hic_split.py input.bam --chunks 100 --prefix split

  # 命令行版本
  hic_split_cli input.bam 100 split

  # 对于14TB Hi-C文件建议直接使用C程序
  hic_split_exec input.bam 100 split 1
```
