#!/usr/bin/env python3
import filecmp
import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "vcf2phylip.py"


def run(args, env=None):
    result = subprocess.run(
        [sys.executable, str(SCRIPT)] + args,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
        check=True,
    )
    return result.stdout


def make_vcf(path):
    text = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t1\t.\tA\tG\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.
chr1\t2\t.\tC\tT\t60\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/1
chr1\t3\t.\tG\tA,C\t60\tPASS\t.\tGT\t0/2\t1/2\t0/0\t2/2
chr1\t4\t.\tAT\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0
chr1\t5\t.\tT\tC\t60\tPASS\t.\tGT\t0/0\t0/1\t9/9\t1/1
chr2\t1\t.\tA\t<NON_REF>\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0
"""
    path.write_text(text, encoding="utf-8")


def compare_outputs(left, right):
    suffixes = [".phy", ".fasta", ".nexus", ".bin.nexus", ".used_sites.tsv"]
    for suffix in suffixes:
        a = next(left.glob("*" + suffix))
        b = next(right.glob("*" + suffix))
        if not filecmp.cmp(a, b, shallow=False):
            raise AssertionError("Output differs: {}".format(suffix))


def main():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        vcf = root / "small.vcf"
        make_vcf(vcf)
        gz = root / "small.vcf.gz"
        with vcf.open("rb") as source, gzip.open(gz, "wb") as target:
            target.write(source.read())

        out1 = root / "one"
        out4 = root / "four"
        common = ["-m", "3", "-f", "-n", "-b", "-w", "-o", "S4"]
        run(["-i", str(vcf), "-t", "1", "--output-folder", str(out1)] + common)
        run(["-i", str(gz), "-t", "4", "--chunk-size-mb", "8",
             "--output-folder", str(out4)] + common)
        compare_outputs(out1, out4)

        env = os.environ.copy()
        env["SLURM_CPUS_PER_TASK"] = "2"
        auto_out = run(["-i", str(vcf), "-p", "-b",
                        "--output-folder", str(root / "auto")], env=env)
        if "Parallel workers: 2 (100% auto-detected from SLURM_CPUS_PER_TASK; not reduced by VCF size)" not in auto_out:
            raise AssertionError("Scheduler CPU auto-detection failed")
        if "VCF size:" not in auto_out or "target about" not in auto_out:
            raise AssertionError("Adaptive VCF profiling output missing")

        manual_out = run(["-i", str(vcf), "-p", "-b", "-t", "2",
                          "--chunk-size-mb", "12",
                          "--output-folder", str(root / "manual")])
        if "user byte limit" not in manual_out:
            raise AssertionError("Manual chunk byte override was not reported")

        run(["-i", str(vcf), "-r", "-t", "2",
             "--output-folder", str(root / "random")])

    print("All tests passed")


if __name__ == "__main__":
    main()
