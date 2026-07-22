# vcf2phylip 2.9-mt6 并行优化版

该版本基于 `edgardomortiz/vcf2phylip` v2.9，在保留原有参数、过滤逻辑和输出格式的基础上，加入多进程解析、矩阵块式转置、CPU 自动检测、智能分块，以及针对不同压缩格式的自动加速后端。

## 自动输入路径

程序不再仅根据 `.gz` 扩展名判断格式，而是读取文件头：

| 实际输入格式 | 自动处理方式 |
|---|---|
| 普通 `.vcf` | 二进制流读取 + 多进程基因型解析 |
| 普通 gzip | 优先 `python-isal`，其次 `python-zlib-ng`，最后回退 Python `gzip` |
| BGZF、无 TBI/CSI | `bgzip -@` 多线程解压 + 多进程解析 |
| BGZF、有 TBI/CSI | `tabix` 按染色体/窗口并行读取、解压和解析 |

只有在 gzip extra field 中找到有效的 BGZF `BC` 子字段，并且首个 BGZF block 可以正确解压时，程序才会判定为 BGZF。文件名是 `.vcf.gz` 或旁边存在 `.tbi` 并不足以证明它是 BGZF。

## BGZF 索引并行模式

自动启用必须同时满足：

1. 输入文件经验证确实是 BGZF；
2. 存在 `<文件>.tbi` 或 `<文件>.csi`；
3. 系统可以找到 `tabix`；
4. `tabix -l` 可以正常读取索引中的染色体名称。

程序按照索引中的染色体顺序处理。若 VCF header 中含有 `##contig=<ID=...,length=...>`，大染色体会继续拆成有序窗口。由于 tabix 是区间重叠查询，程序会再按每条记录的 `POS` 过滤，使跨窗口的 `END`/结构变异记录不会被重复统计。

如果自动索引模式运行失败，例如索引过期或损坏，程序会清空部分临时结果，并自动回退到 BGZF 流式模式。手动指定 `--input-backend tabix` 时则直接报告错误，便于排查。

## 依赖

必需：

- Python 3.8 或更高版本

可选加速组件：

- HTSlib 的 `bgzip`、`tabix`：强烈建议用于 BGZF
- `python-isal`：加速普通 gzip
- `python-zlib-ng`：普通 gzip 的另一种加速后端

不安装任何可选组件时，程序仍可通过 Python 标准库正常运行。

```bash
mamba install -c conda-forge -c bioconda htslib python-isal python-zlib-ng
```

## CPU 含义

不指定 `-t` 时，程序使用当前任务允许的全部 CPU：

```text
workers = 调度器/CPU affinity/操作系统允许的最大 CPU 数
```

VCF 大小只用于调整分块，不会主动减少 CPU 数。

```bash
python3 vcf2phylip.py -i input.vcf.gz -m 380 -f
```

手动指定总 CPU 预算：

```bash
python3 vcf2phylip.py -i input.vcf.gz -m 380 -f -t 32
```

BGZF 无索引流式模式会把这 32 个 CPU 分配给两部分，例如：

```text
6 个 bgzip 解压线程 + 26 个基因型解析进程 = 32
```

BGZF 有索引时，32 个 worker 分别执行独立的 tabix 区域读取、解压和解析。

## 新增参数

```text
--input-backend {auto,plain,stdlib,isal,zlib-ng,bgzip,tabix}
    手动覆盖自动后端，主要用于测试或排错。

--decompression-threads N
    为 BGZF 流式 bgzip 解压保留 N 个线程。

--no-indexed-regions
    即使存在 TBI/CSI，也不使用区域并行。
```

原有 `-i/-m/-o/-p/-f/-n/-b/-r/-w/--output-folder/--output-prefix` 均保留。

## 使用示例

```bash
# 自动识别格式和后端
python3 vcf2phylip.py \
  -i population.vcf.gz \
  -m 380 \
  -f \
  -w \
  -t 32
```

```bash
# 强制使用 TBI/CSI 区域并行
python3 vcf2phylip.py \
  -i population.vcf.gz \
  --input-backend tabix \
  -t 32 \
  -f
```

```bash
# 存在索引，但本次只使用 bgzip 流式解压
python3 vcf2phylip.py \
  -i population.vcf.gz \
  --no-indexed-regions \
  --decompression-threads 6 \
  -t 32 \
  -f
```

## 准备 BGZF 和索引

```bash
bcftools sort input.vcf -Oz -o input.sorted.vcf.gz

# TBI
bcftools index -t input.sorted.vcf.gz

# 或 CSI，较长染色体建议使用 CSI
bcftools index -c input.sorted.vcf.gz
```

普通 gzip 不能通过改名或创建空 `.tbi` 变成 BGZF，必须使用 `bgzip` 或 `bcftools -Oz` 重新压缩。

## 功能兼容性

继续支持：

- PHYLIP、FASTA、NEXUS、binary NEXUS
- 任意倍性核苷酸矩阵
- 二倍体 biallelic SNP 的 SNAPP 编码
- IUPAC 杂合编码
- `--resolve-IUPAC`
- 缺失样本数过滤
- 外群优先输出
- used-sites 坐标文件
- 自定义目录和前缀

未启用 `--resolve-IUPAC` 时，不同 worker 数或不同输入后端生成的矩阵应逐字节一致。

## 测试

```bash
python3 tests/test_regression.py
python3 tests/test_backends.py
```

测试覆盖普通 VCF、普通 gzip、BGZF 无索引、BGZF+TBI/CSI、跨窗口重叠记录、防重复、索引失败自动回退，以及全部矩阵格式和原有过滤功能。
