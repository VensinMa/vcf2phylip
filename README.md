# vcf2phylip (Multithreaded)

Convert SNPs in VCF format to PHYLIP, NEXUS, binary NEXUS, or FASTA alignments for phylogenetic analysis.

**Multithreaded fork** of [edgardomortiz/vcf2phylip](https://github.com/edgardomortiz/vcf2phylip) (v2.9) with parallel processing and adaptive chunking.

## What's new

Added `-t/--threads` and `--chunk-size-mb` parameters. All other parameters and output are **identical** to the original.

### Parallelized stages

| Stage | Description |
|-------|-------------|
| **Phase 1** | VCF parsing + genotype conversion — adaptive chunks profiled from VCF, streamed with backpressure via `ProcessPoolExecutor` |
| **Phase 2** | Matrix assembly — transposed blocks written during Phase 1, per-sample `seek()` via `ThreadPoolExecutor` |

### Adaptive chunking

The script profiles the input VCF (samples 32MB prefix) to estimate total size, record count, and average row width. Chunk limits are then calculated from CPU count and workload:

- Targets 8–12 tasks per worker depending on SNP count
- Compressed VCFs get slightly fewer tasks (gzip decompression is sequential)
- Dual limit: chunk ends at byte target OR record target, whichever comes first
- Override with `--chunk-size-mb` if needed

### Smart CPU detection

Automatically detects available CPUs respecting:
- HPC scheduler limits (SLURM, PBS, LSF)
- Linux CPU affinity / cgroup / cpuset
- Python 3.13+ `os.process_cpu_count()`

## Usage

```bash
# Auto-detect CPUs (default), adaptive chunking
python3 vcf2phylip.py -i myfile.vcf -f

# Specify worker count
python3 vcf2phylip.py -i myfile.vcf -f -t 8

# Single-threaded (original behavior)
python3 vcf2phylip.py -i myfile.vcf -f -t 1

# Override chunk size
python3 vcf2phylip.py -i myfile.vcf -f --chunk-size-mb 64
```

## Full usage

```
usage: vcf2phylip.py [-h] -i FILENAME [--output-folder FOLDER]
                     [--output-prefix PREFIX] [-m MIN_SAMPLES_LOCUS]
                     [-o OUTGROUP] [-p] [-f] [-n] [-b] [-r] [-w]
                     [-t THREADS] [--chunk-size-mb CHUNK_SIZE_MB] [-v]

optional arguments:
  -i FILENAME, --input FILENAME
  --output-folder FOLDER
  --output-prefix PREFIX
  -m MIN_SAMPLES_LOCUS, --min-samples-locus MIN_SAMPLES_LOCUS
  -o OUTGROUP, --outgroup OUTGROUP
  -p, --phylip-disable
  -f, --fasta
  -n, --nexus
  -b, --nexus-binary
  -r, --resolve-IUPAC
  -w, --write-used-sites
  -t THREADS, --threads THREADS         Parallel workers (default: 100% of available CPUs)
  --chunk-size-mb CHUNK_SIZE_MB         Override auto chunk size in MiB
  -v, --version
```

## Example output

```
Converting file '412samples.vcf.gz':

Number of samples in VCF: 412
Parallel workers: 32 (100% auto-detected from CPU affinity; not reduced by VCF size)
VCF size: 1.9 GiB stored; about 12.3 GiB uncompressed (gzip/BGZF prefix estimate)
Estimated data records: 1,234,567; average sampled row: 10,967 bytes
VCF chunk limit: 64 MiB or 15,000 data records; target about 320 tasks (auto from CPU count and estimated VCF workload)
```

## Benchmark

Tested on i9-14900KF (32 cores), 1,080,920 SNPs × 412 samples:

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

**Tip:** If your pipeline already decompresses VCF for other tools, run vcf2phylip on the uncompressed file for maximum speed.

## Credits

- Original code: [Edgardo M. Ortiz](https://github.com/edgardomortiz)
- Multithreaded version: [Ma Wenxin](https://github.com/VensinMa)

## Citation

Original: Ortiz, E.M. 2019. vcf2phylip v2.0. DOI:10.5281/zenodo.2540861
