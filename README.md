# vcf2phylip (multithreaded, adaptive compressed-input backends)

Convert SNP genotypes in VCF format to relaxed PHYLIP, FASTA, NEXUS, or binary NEXUS matrices for phylogenetic analysis.

This repository is a performance-oriented fork of [`edgardomortiz/vcf2phylip`](https://github.com/edgardomortiz/vcf2phylip) v2.9. Version **2.9-mt6** preserves the original matrix formats and filtering behavior while adding multiprocessing, optimized matrix transposition, automatic CPU detection, adaptive chunking, and compressed-input backend selection.

## What mt6 adds

The input path is selected from the actual file bytes rather than from the `.gz` suffix:

| Input detected | Automatic processing path |
|---|---|
| Plain VCF | Binary stream → multiprocessing genotype parsing |
| Ordinary gzip | `python-isal` → `python-zlib-ng` → stdlib `gzip`, in that priority order |
| BGZF without TBI/CSI | HTSlib `bgzip -@` multithreaded decompression → multiprocessing parsing |
| BGZF with TBI/CSI | Ordered chromosome/window queries through `tabix`; decompression and parsing run in parallel |

A `.vcf.gz` file is treated as BGZF only when its gzip extra field contains a valid `BC` BGZF subfield and the first block can be decompressed successfully. Merely having a `.gz` suffix or a `.tbi` file is not enough.

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

--input-backend {auto,plain,stdlib,isal,zlib-ng,bgzip,tabix}
    Override backend selection for benchmarking or troubleshooting.

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
- plain VCF and ordinary gzip
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

Tested on i9-14900KF (32 cores), 1,080,920 SNPs × 412 samples (v2.9-mt5):

### Gzipped VCF (1.9 GB)

| Threads | Time (s) | Speedup |
|---------|----------|---------|
| original v2.9 | 247.00 | 1.00× |
| mt5 -t 1 | 159.89 | 1.54× |
| mt5 -t 2 | 71.18 | 3.47× |
| mt5 -t 4 | 37.67 | **6.55×** |
| mt5 -t 8 | 37.25 | 6.63× |
| mt5 -t 16 | 36.64 | 6.74× |

> Gzip decompression is sequential — bottleneck caps at ~6.6× with 4 threads.
> mt6 with tabix index can bypass this limit via region-parallel decompression.

### Uncompressed VCF (12 GB, averaged over 2 runs)

| Threads | Time (s) | Speedup |
|---------|----------|---------|
| original v2.9 | 220.19 | 1.00× |
| mt5 -t 1 | 136.91 | 1.61× |
| mt5 -t 2 | 70.08 | 3.14× |
| mt5 -t 4 | 36.42 | 6.05× |
| mt5 -t 8 | 21.55 | **10.22×** |
| mt5 -t 16 | 17.27 | **12.75×** |

> Without gzip bottleneck, scales linearly to 12.8× at 16 threads.

## Credits and citation

- Original code: Edgardo M. Ortiz
- Original data/testing: Juan D. Palacio-Mejía
- Multithreaded fork: Ma Wenxin

Please cite the original software:

> Ortiz, E.M. 2019. vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. DOI: 10.5281/zenodo.2540861
