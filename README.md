# vcf2phylip (Multithreaded)

Convert SNPs in VCF format to PHYLIP, NEXUS, binary NEXUS, or FASTA alignments for phylogenetic analysis.

**Multithreaded fork** of [edgardomortiz/vcf2phylip](https://github.com/edgardomortiz/vcf2phylip) (v2.9) with parallel processing support.

## What's new

Added `-t / --threads` parameter for multiprocessing. All other parameters and output are **identical** to the original.

### Parallelized stages

| Stage | Description |
|-------|-------------|
| **Phase 1** | VCF parsing + genotype conversion — batches of 10,000 lines processed in parallel |
| **Phase 2** | Matrix transposition + output — per-sample extraction parallelized, file read only once |

## Usage

```bash
# Original single-threaded (default, same as upstream)
python3 vcf2phylip.py -i myfile.vcf -f

# Parallel: 4 / 8 / 16 processes
python3 vcf2phylip.py -i myfile.vcf -f -t 4
python3 vcf2phylip.py -i myfile.vcf -f -t 8
python3 vcf2phylip.py -i myfile.vcf -f -t 16
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
  -t THREADS, --threads THREADS    Number of parallel processes (default=1)
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
- Multithreaded version: [Ma Wenxin](https://github.com/mawenxin)

## Citation

Original: Ortiz, E.M. 2019. vcf2phylip v2.0. DOI:10.5281/zenodo.2540861
