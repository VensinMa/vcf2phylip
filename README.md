# vcf2phylip (multithreaded, adaptive input backends)

Convert SNP genotypes in VCF format to relaxed PHYLIP, FASTA, NEXUS, or binary NEXUS matrices for phylogenetic analysis.

This repository is a performance-oriented fork of [`edgardomortiz/vcf2phylip`](https://github.com/edgardomortiz/vcf2phylip) v2.9. Version **2.9-mt7** preserves the original matrix formats and filtering behavior while adding multiprocessing, optimized matrix transposition, automatic CPU detection, adaptive chunking, direct plain-VCF range readers, and compressed-input backend selection.

## Automatic input paths

The input path is selected from the actual file bytes rather than from the `.gz` suffix:

| Input detected | Automatic processing path |
|---|---|
| Plain VCF, one worker | Sequential binary reader |
| Plain VCF, multiple workers | Independent aligned byte ranges → direct `pread`/seek → parallel parsing |
| Ordinary gzip | `python-isal` → `python-zlib-ng` → stdlib `gzip`, in that priority order |
| BGZF without TBI/CSI | HTSlib `bgzip -@` multithreaded decompression → multiprocessing parsing |
| BGZF with TBI/CSI | Ordered chromosome/window queries through `tabix`; decompression and parsing run in parallel |

A `.vcf.gz` file is treated as BGZF only when its gzip extra field contains a valid `BC` BGZF subfield and the first block can be decompressed successfully. Merely having a `.gz` suffix or a `.tbi` file is not enough.

### Direct plain-VCF mode

With more than one worker, an uncompressed VCF is divided into ordered raw
byte ranges. Each worker opens the file independently, aligns its nominal
start and end offsets to complete newline-delimited VCF records, reads the
aligned range directly, parses it, and writes a local transposed matrix shard.
The parent process merges shards in byte-range order.

This removes the mt6 plain-input bottleneck where one parent process read the
whole VCF and serialized large text chunks through multiprocessing queues.
Range size is derived from the same CPU/workload-aware chunk profile and is
bounded to avoid a few oversized tasks. The implementation also handles a
single VCF record spanning one or more nominal range boundaries without
duplication or omission.

Use `--input-backend plain-stream` to retain the mt6-style sequential input
reader for benchmarking or troubleshooting.

### Indexed BGZF mode

Indexed mode is enabled automatically only when all of the following are true:

1. The file is verified as BGZF.
2. A sidecar `<file>.tbi` or `<file>.csi` exists.
3. The `tabix` executable is available.
4. `tabix -l` can read sequence names from the index.

Contigs are processed in index order. When `##contig` lengths are available, large contigs are divided into ordered windows. Returned records are filtered by their VCF `POS`, preventing overlap-based tabix queries from duplicating records that span a window boundary. If indexed mode fails in automatic mode, the program discards partial temporary data and retries with sequential BGZF streaming.

HTSlib documents BGZF as concatenated gzip-compatible blocks smaller than 64 KiB and supports multithreaded `bgzip` operation with `-@`. Tabix requires position-sorted BGZF input and a TBI or CSI index for region retrieval:

- <https://www.htslib.org/doc/bgzip.html>
- <https://www.htslib.org/doc/tabix.html>

## Compatibility

The original options and outputs remain available:

- plain VCF and compressed VCF
- arbitrary-ploidy nucleotide matrices
- relaxed PHYLIP, FASTA, and NEXUS
- diploid biallelic binary NEXUS for SNAPP
- IUPAC heterozygous consensus
- random `--resolve-IUPAC`
- `--min-samples-locus`
- outgroup-first output
- `--write-used-sites`
- output folder and prefix selection

When `--resolve-IUPAC` is not used, worker count and input backend do not change matrix content or ordering.

## Requirements

Required:

- Python 3.8 or newer

Optional accelerators:

- HTSlib `bgzip` and `tabix`: recommended for BGZF input
- [`python-isal`](https://python-isal.readthedocs.io/): faster ordinary-gzip decompression
- [`python-zlib-ng`](https://python-zlib-ng.readthedocs.io/): ordinary-gzip fallback acceleration

The program remains functional with the Python standard library alone.

Example installations:

```bash
# Conda / mamba: full compressed-input acceleration
mamba install -c conda-forge -c bioconda htslib python-isal python-zlib-ng

# Pip: optional ordinary-gzip accelerators only
python3 -m pip install isal zlib-ng
```

## Basic usage

```bash
python3 vcf2phylip.py -i input.vcf.gz -m 380 -f -w
```

Without `-t`, the program uses 100% of the CPUs available to the current process or scheduler allocation. Input size changes chunk sizing but never reduces this CPU ceiling.

```bash
# Explicit total CPU budget
python3 vcf2phylip.py -i input.vcf.gz -m 380 -f -t 32
```

For a BGZF stream without an index, part of the `-t` budget is assigned to `bgzip` and the remainder to parser processes. For indexed BGZF, all workers are independent tabix region workers.

## New performance options

```text
-t, --threads N
    Total CPU budget. Default: all CPUs available to the job.

--chunk-size-mb N
    Override the automatically selected chunk byte limit.

--input-backend {auto,plain,plain-stream,stdlib,isal,zlib-ng,bgzip,tabix}
    Override backend selection for benchmarking or troubleshooting.

    plain        Direct byte-range mode when -t > 1; sequential when -t 1.
    plain-stream Force one sequential plain-VCF reader feeding parser workers.

--decompression-threads N
    Threads reserved for the streaming bgzip backend. The value must be
    smaller than the total -t budget when -t > 1.

--no-indexed-regions
    Ignore TBI/CSI region parallelism and use the best streaming backend.
```

Examples:

```bash
# Force indexed BGZF processing
python3 vcf2phylip.py -i input.vcf.gz -t 32 --input-backend tabix -f

# Use BGZF streaming even though an index exists
python3 vcf2phylip.py -i input.vcf.gz -t 32 --no-indexed-regions -f

# Reserve 6 of 32 CPUs for bgzip decompression
python3 vcf2phylip.py -i input.vcf.gz -t 32 \
  --no-indexed-regions --decompression-threads 6 -f

# Force the dependency-free gzip reader
python3 vcf2phylip.py -i input.vcf.gz --input-backend stdlib -f

# Compare the optimized plain-VCF reader with the old streaming path
python3 vcf2phylip.py -i input.vcf -t 16 --input-backend plain -f
python3 vcf2phylip.py -i input.vcf -t 16 --input-backend plain-stream -f
```

## Startup report

The selected path is printed before conversion:

```text
Parallel workers: 32 (100% auto-detected from SLURM_CPUS_PER_TASK; not reduced by VCF size)
Detected input format: BGZF
Input backend: indexed BGZF parallel regions through tabix
Indexed sidecar: population.vcf.gz.csi
Region workers: 32
VCF size: 4.8 GiB stored; about 20.1 GiB uncompressed (gzip/BGZF prefix estimate)
Estimated data records: 10,120,000; average sampled row: 2,134 bytes
VCF chunk limit: 53 MiB or 26,000 data records
```

For streaming BGZF, CPU allocation is reported separately:

```text
CPU allocation: 6 bgzip threads + 26 parser workers = 32 total
```

For an uncompressed VCF with multiple workers:

```text
Input backend: plain VCF direct parallel byte ranges
Direct range workers: 16
Direct byte ranges: 128
```

## Preparing indexed BGZF input

```bash
# Input VCF must be coordinate sorted before indexing
bcftools sort input.vcf -Oz -o input.sorted.vcf.gz

# TBI is suitable when all coordinates fit its range
bcftools index -t input.sorted.vcf.gz

# CSI is safer for very long chromosomes
bcftools index -c input.sorted.vcf.gz
```

An ordinary gzip file cannot become indexed merely by renaming it or creating an empty `.tbi`; it must first be recompressed as BGZF.

## SLURM example

```bash
#!/bin/bash
#SBATCH --job-name=vcf2phylip
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=vcf2phylip.%j.out
#SBATCH --error=vcf2phylip.%j.err

python3 vcf2phylip.py \
  -i population.filtered.vcf.gz \
  -m 380 \
  -f \
  -w \
  --output-folder results \
  --output-prefix population_snps
```

`SLURM_CPUS_PER_TASK` is detected automatically, so `-t` can be omitted.

## Tests

The included regression tests cover:

- single-worker and multi-worker output equality
- plain VCF direct ranges, forced sequential plain input, and ordinary gzip
- plain-range boundaries falling inside a multi-megabyte VCF record
- verified BGZF with and without an index
- indexed chromosome/window order
- overlap records spanning region boundaries
- automatic fallback from a failing/stale index
- PHYLIP, FASTA, NEXUS, binary NEXUS, and used-sites output
- outgroup order, missingness filtering, MNP exclusion, multiallelic SNPs, `<NON_REF>`, and malformed genotypes

Run:

```bash
python3 tests/test_regression.py
python3 tests/test_backends.py
```

## Benchmark

Tested on i9-14900KF, 1,080,920 SNPs × 412 samples (v2.9-mt7). System resources idle during test.

> **Tip:** For compressed VCF, build a `.tbi` or `.csi` index beforehand to unlock tabix region-parallel mode (up to 30× speedup):
> ```bash
> bcftools sort input.vcf -Oz -o input.sorted.vcf.gz
> bcftools index -t input.sorted.vcf.gz
> ```

### 1. Compressed VCF + .tbi index (tabix region-parallel)

| Threads | Time (s) | Speedup |
|---------|----------|---------|
| original v2.9 | 241.11 | 1.00× |
| mt7 -t 8 | 12.91 | **18.68×** |
| mt7 -t 16 | 8.54 | **28.23×** |
| mt7 -t auto | 7.11 | **33.91×** |

### 2. Compressed VCF + no index (BGZF streaming)

| Threads | Time (s) | Speedup |
|---------|----------|---------|
| original v2.9 | 241.11 | 1.00× |
| mt7 -t 8 | 14.80 | **16.29×** |
| mt7 -t 16 | 14.81 | **16.28×** |
| mt7 -t auto | 14.88 | **16.20×** |

### 3. Uncompressed VCF (12 GB)

| Threads | Time (s) | Speedup |
|---------|----------|---------|
| original v2.9 | 218.51 | 1.00× |
| mt7 -t 8 | 9.69 | **22.55×** |
| mt7 -t 16 | 6.65 | **32.86×** |
| mt7 -t auto | 5.52 | **39.59×** |

## Credits and citation

- Original code: Edgardo M. Ortiz
- Original data/testing: Juan D. Palacio-Mejía
- Multithreaded fork: Ma Wenxin

Please cite the original software:

> Ortiz, E.M. 2019. vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. DOI: 10.5281/zenodo.2540861
