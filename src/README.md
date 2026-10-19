# LJA-Bridge

LJA-Bridge 是基于 [La Jolla Assembler (LJA)](https://github.com/AntonBankevich/LJA) 的基因组组装工具。它保留 LJA 面向 PacBio HiFi 的 de Bruijn 图组装流程，并增加一轮保守的 **ONT gap bridging**：用 Oxford Nanopore 长读段把相邻 contig 两端锚定起来，再把中间间隔填回去，从而提高连续性。

当前版本：`0.1`

## 相对 LJA 的改动

LJA 用很大的 k 构建 de Bruijn 图，从而在哺乳动物基因组尺度上自动解开大部分重复。完整流程仍由三个核心模块组成：

1. **jumboDBG**：支持任意大 k 的 de Bruijn 图构建
2. **mowerDBG**：面向 HiFi 的纠错
3. **multiDBG**：在基因组不同区域使用可变 k，并利用 HiFi 全长做重复解析

随后 LJA 会对同聚物压缩后的结果做解压和抛光。LJA-Bridge 在此基础上增加：

- `--ont-reads`：输入 ONT 读段
- `--use-ont-in-polish`：在抛光阶段使用 ONT，并启用抛光后的 ONT gap bridging
- 保守桥接策略：要求足够的跨端 ONT 支持、端点锚定重叠和路径不冲突，避免把错误邻接拼进最终序列

当前 **ONT-guided MDBG** 已在代码中暂时关闭，不会改变 multiDBG 阶段的图解析行为。

## 依赖

- 64 位 Linux
- CMake 3.1+
- 支持 C++17 的编译器（GCC 7+）
- GNU make
- zlib
- OpenMP

推荐在 Linux 或 WSL 中编译。

## 编译

在源码目录下执行：

```bash
cmake -S . -B build
cmake --build build -j
```

也可以沿用原 LJA 的写法：

```bash
cmake .
make -j
```

编译完成后，主程序位于：

- out-of-source：`build/bin/LJA-Bridge`
- 若在源码树内构建：`bin/LJA-Bridge`

`jumboDBG` 仍会作为独立可执行文件一起生成。

安装到源码树的 `bin/`（`CMAKE_INSTALL_PREFIX` 默认指向源码目录）：

```bash
cmake --install build
```

请将包含 `LJA-Bridge` 的目录加入 `PATH`，或使用绝对路径调用。

## 输入

读段可以是 fasta / fastq，也可以是 gzip 压缩文件。扩展名必须和实际格式一致，例如：

- 未压缩 fasta：`.fa` / `.fasta`
- 压缩 fastq：`.fq.gz` / `.fastq.gz`

HiFi 读段通过 `--reads` 提供，可重复指定多个文件。ONT 读段通过 `--ont-reads` 提供，同样可以指定多个文件。

## 用法

仅用 HiFi（与 LJA 相同）：

```bash
LJA-Bridge -o <output_dir> --reads <hifi.fastq>
```

HiFi + ONT 抛光与 gap bridging：

```bash
LJA-Bridge -o <output_dir> \
    --reads <hifi.fastq> \
    --ont-reads <ont.fastq> \
    --use-ont-in-polish \
    -t 32
```

查看全部参数：

```bash
LJA-Bridge --help
```

### 常用参数

| 参数 | 说明 |
| --- | --- |
| `-o` / `--output-dir` | 输出目录 |
| `--reads` | HiFi 读段文件，可重复指定 |
| `--ont-reads` | ONT 读段文件，可重复指定 |
| `--use-ont-in-polish` | 在抛光中使用 ONT，并启用 ONT gap bridging |
| `-t` / `--threads` | 线程数，默认 16 |
| `-k` | 第一轮纠错使用的 k，应为奇数 |
| `-K` | 第二轮纠错和 multiDBG 初始化使用的 k，应为奇数 |
| `--diploid` | 二倍体基因组；默认按单倍体 / 近交系处理 |

只提供 `--ont-reads`、不打开 `--use-ont-in-polish` 时，组装仍按 HiFi 流程完成，ONT 不会进入抛光，也不会做 gap bridging。

## 输出

主要结果写在 `-o` 指定的目录中：

- `assembly.fasta`：最终组装序列（若启用 ONT gap bridging，已合并的 contig 会反映在该文件中）
- `mdbg.gfa`：最终 multiplex de Bruijn 图（实验性）
- `mdbg/assembly.hpc.fasta`：同聚物解压前的组装
- `mdbg/mdbg.hpc.gfa`：同聚物压缩边上的最终图
- `dbg.log`：运行日志

被 ONT 桥接合并的 contig 标识中会带有 `_ONTgap_`。

## 许可与引用

本工具沿用 LJA 的 BSD 许可，详见 [LICENSE](LICENSE)。

如果你在研究中使用 LJA-Bridge，请同时引用原始 LJA 工作：

- [LJA preprint](https://www.biorxiv.org/content/10.1101/2020.12.10.420448)

LJA 由 Anton Bankevich、Andrey Bzikadze（[Pavel Pevzner lab, UCSD](http://cseweb.ucsd.edu/~ppevzner/)）和 Dmitry Antipov（[Center for Algorithmic Biology, SPbSU](https://cab.spbu.ru/)）开发。LJA-Bridge 在此代码基础上增加了 ONT gap bridging。
