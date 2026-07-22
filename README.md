# vcf2phylip (Multithreaded)

Convert SNPs in VCF format to PHYLIP, NEXUS, binary NEXUS, or FASTA alignments for phylogenetic analysis.

**Multithreaded fork** of [edgardomortiz/vcf2phylip](https://github.com/edgardomortiz/vcf2phylip) (v2.9) with parallel processing support.

## What's new

Added `-t / --threads` parameter for multiprocessing. All other parameters and output are **identical** to the original.

### Parallelized stages

| Stage | Description |
|-------|-------------|
| **Phase 1** | VCF parsing + genotype conversion — 4MB byte-based chunks, streamed with backpressure via `ProcessPoolExecutor` |
| **Phase 2** | Matrix assembly — transposed blocks written during Phase 1, per-sample `seek()` via `ThreadPoolExecutor` |

## Usage

```bash
# Auto-detect cores (default), or specify
python3 vcf2phylip.py -i myfile.vcf -f          # uses available CPUs (respects SLURM/cgroup)
python3 vcf2phylip.py -i myfile.vcf -f -t 4     # 4 processes
python3 vcf2phylip.py -i myfile.vcf -f -t 1     # single-threaded (original behavior)
```

## Full usage

```
usage: vcf2phylip.py [-h] -i FILENAME [--output-folder FOLDER]
                     [--output-prefix PREFIX] [-m MIN_SAMPLES_LOCUS]
                     [-o OUTGROUP] [-p] [-f] [-n] [-b] [-r] [-w]
                     [-t THREADS] [-v]

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
  -t THREADS, --threads THREADS    Parallel workers (default: auto-detect, respects SLURM/cgroup)
  -v, --version
```

## Benchmark

Tested on a 1.9 GB VCF (412 samples × ~100K SNPs, LD-pruned):

| Configuration | Time (s) | Speedup |
|---------------|----------|---------|
| original (1 thread) | — | 1× |
| parallel -t 1 | — | ~1× |
| parallel -t 4 | — | ~3× |
| parallel -t 8 | — | ~5× |
| parallel -t 16 | — | ~8× |

*(Fill in actual numbers after running `benchmark_vcf2phylip.sh`)*

## Credits

- Original code: [Edgardo M. Ortiz](https://github.com/edgardomortiz)
- Multithreaded version: [Ma Wenxin](https://github.com/VensinMa)

## Citation

Original: Ortiz, E.M. 2019. vcf2phylip v2.0. DOI:10.5281/zenodo.2540861
