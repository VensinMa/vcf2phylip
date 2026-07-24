# vcf2phylip 2.9-mt7 并行优化版

该版本基于 `edgardomortiz/vcf2phylip` v2.9，在保留原有参数、过滤逻辑和输出格式的基础上，加入多进程解析、矩阵块式转置、CPU 自动检测、智能分块、普通 VCF 直接字节区间并行读取，以及针对不同压缩格式的自动加速后端。

## 基准测试

i9-14900KF，1,080,920 SNPs × 412 样本，5 次平均。

| 场景 | 原版 | -t 8 | -t 16 | -t auto |
|------|------|------|-------|---------|
| 压缩 VCF + .tbi (tabix) | 241.11s (1×) | 12.91s (18.7×) | 8.54s (28.2×) | **7.11s (33.9×)** |
| 压缩 VCF + 无索引 (bgzip) | 241.11s (1×) | 14.80s (16.3×) | 14.81s (16.3×) | 14.88s (16.2×) |
| 未压缩 VCF (12 GB) | 218.51s (1×) | 9.69s (22.6×) | 6.65s (32.9×) | **5.52s (39.6×)** |

> **提示：** 压缩 VCF 建议提前构建 `.tbi` 或 `.csi` 索引，可解锁 tabix 区域并行模式（最高 30× 加速）：
> ```bash
> bcftools sort input.vcf -Oz -o input.sorted.vcf.gz
> bcftools index -t input.sorted.vcf.gz
> ```

## 自动输入路径

程序不再仅根据 `.gz` 扩展名判断格式，而是读取文件头：

| 实际输入格式 | 自动处理方式 |
|---|---|
| 普通 `.vcf`，单 worker | 顺序二进制读取 |
| 普通 `.vcf`，多 worker | 独立字节区间 + `pread`/seek 直接读取 + 并行解析 |
| 普通 gzip | 优先 `python-isal`，其次 `python-zlib-ng`，最后回退 Python `gzip` |
| BGZF、无 TBI/CSI | `bgzip -@` 多线程解压 + 多进程解析 |
| BGZF、有 TBI/CSI | `tabix` 按染色体/窗口并行读取、解压和解析 |

只有在 gzip extra field 中找到有效的 BGZF `BC` 子字段，并且首个 BGZF block 可以正确解压时，程序才会判定为 BGZF。文件名是 `.vcf.gz` 或旁边存在 `.tbi` 并不足以证明它是 BGZF。

## 普通 VCF 直接区间并行

当普通未压缩 VCF 使用两个或更多 worker 时，程序不再由主进程顺序读取整个文件并把大块文本发送给子进程，而是：

1. 根据文件大小、CPU 数和估算记录数生成有序字节区间；
2. 每个 worker 独立打开同一个 VCF；
3. 把理论起止位置对齐到完整换行记录；
4. 使用 `os.pread()`，不支持时回退 `seek/read`，直接读取自己的区间；
5. 在 worker 内完成解析和矩阵块转置；
6. 主进程按照原始字节区间顺序合并结果。

这样可以消除 mt6 普通 VCF 路径中的主进程顺序读取和大块 VCF 文本 IPC 传输。即使一条超长 VCF 记录跨越一个或多个理论边界，也不会重复或遗漏。

如需测试旧式顺序读取路径，可使用：

```bash
python3 vcf2phylip.py \
  -i input.vcf \
  -t 16 \
  --input-backend plain-stream
```

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
--input-backend {auto,plain,plain-stream,stdlib,isal,zlib-ng,bgzip,tabix}
    手动覆盖自动后端，主要用于测试或排错。

    plain
        多 worker 时使用普通 VCF 直接字节区间；-t 1 时顺序读取。

    plain-stream
        强制使用单个顺序读取器，再把数据块发送给解析进程。

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

```bash
# 普通 VCF 新旧读取路径对比
python3 vcf2phylip.py -i population.vcf -t 16 --input-backend plain
python3 vcf2phylip.py -i population.vcf -t 16 --input-backend plain-stream
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

测试覆盖普通 VCF 直接区间、强制顺序读取、边界落在数 MiB 超长记录内部、普通 gzip、BGZF 无索引、BGZF+TBI/CSI、跨窗口重叠记录、防重复、索引失败自动回退，以及全部矩阵格式和原有过滤功能。
