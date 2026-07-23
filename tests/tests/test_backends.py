#!/usr/bin/env python3
import filecmp
import gzip
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import zlib
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "vcf2phylip.py"


def run(args, env=None):
    result = subprocess.run(
        [sys.executable, str(SCRIPT)] + args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        check=True,
    )
    return result.stdout


def make_vcf(path):
    text = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=3000000>
##contig=<ID=chr2,length=2000000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t1\t.\tA\tG\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.
chr1\t999999\t.\tA\t<DEL>\t60\tPASS\tEND=1000002\tGT\t0/0\t0/1\t1/1\t0/0
chr1\t1000001\t.\tC\tT\t60\tPASS\t.\tGT\t0/1\t0/0\t1/1\t0/1
chr1\t2000001\t.\tG\tA,C\t60\tPASS\t.\tGT\t0/2\t1/2\t0/0\t2/2
chr1\t2500000\t.\tAT\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0
chr1\t2999999\t.\tT\tC\t60\tPASS\t.\tGT\t0/0\t0/1\t9/9\t1/1
chr2\t1\t.\tA\t<NON_REF>\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0
chr2\t1500000\t.\tC\tG\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0
"""
    path.write_text(text, encoding="utf-8")


def make_boundary_vcf(path):
    """Create a plain VCF whose byte-range boundary lands inside a long row."""
    long_info = "X=" + ("A" * (5 * 1024 * 1024))
    with path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##contig=<ID=chr1,length=1000>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n")
        handle.write("chr1\t1\t.\tA\tG\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\t./.\n")
        handle.write("chr1\t2\t.\tC\tT\t60\tPASS\t{}\tGT\t0/1\t0/0\t1/1\t0/1\n".format(long_info))
        handle.write("chr1\t3\t.\tG\tA\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0\n")


def make_bgzf(source, target):
    raw = source.read_bytes()

    def block(data):
        compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
        payload = compressor.compress(data) + compressor.flush()
        total = 18 + len(payload) + 8
        if total > 65536:
            raise AssertionError("test BGZF block too large")
        header = (
            b"\x1f\x8b\x08\x04" + b"\x00\x00\x00\x00" + b"\x00\xff" +
            b"\x06\x00" + b"BC" + b"\x02\x00" + struct.pack("<H", total - 1)
        )
        trailer = struct.pack("<II", zlib.crc32(data) & 0xffffffff, len(data) & 0xffffffff)
        return header + payload + trailer

    with target.open("wb") as out:
        for offset in range(0, len(raw), 60000):
            out.write(block(raw[offset:offset + 60000]))
        out.write(block(b""))


def make_fake_tools(folder):
    folder.mkdir()
    bgzip = folder / "bgzip"
    bgzip.write_text("""#!/usr/bin/env python3
import gzip, shutil, sys
path = sys.argv[-1]
with gzip.open(path, 'rb') as source:
    shutil.copyfileobj(source, sys.stdout.buffer)
""", encoding="utf-8")
    bgzip.chmod(0o755)

    tabix = folder / "tabix"
    tabix.write_text("""#!/usr/bin/env python3
import gzip, sys
args = sys.argv[1:]
if args[0] == '-l':
    path = args[1]
    seen = []
    with gzip.open(path, 'rb') as handle:
        for line in handle:
            if line.startswith(b'#'):
                continue
            chrom = line.split(b'\\t', 1)[0].decode()
            if chrom not in seen:
                seen.append(chrom)
    sys.stdout.write('\\n'.join(seen) + ('\\n' if seen else ''))
    raise SystemExit(0)
path, region = args[0], args[1]
if ':' in region:
    chrom, interval = region.rsplit(':', 1)
    start_text, end_text = interval.split('-', 1)
    start, end = int(start_text), int(end_text)
else:
    chrom, start, end = region, None, None
with gzip.open(path, 'rb') as handle:
    for line in handle:
        if line.startswith(b'#'):
            continue
        fields = line.split(b'\\t')
        if fields[0].decode() != chrom:
            continue
        pos = int(fields[1])
        record_end = pos + max(1, len(fields[3])) - 1
        for item in fields[7].split(b';'):
            if item.startswith(b'END='):
                record_end = int(item[4:])
        if start is not None and (record_end < start or pos > end):
            continue
        sys.stdout.buffer.write(line)
""", encoding="utf-8")
    tabix.chmod(0o755)


def compare_outputs(left, right):
    for suffix in (".phy", ".fasta", ".nexus", ".bin.nexus", ".used_sites.tsv"):
        a = next(left.glob("*" + suffix))
        b = next(right.glob("*" + suffix))
        if not filecmp.cmp(a, b, shallow=False):
            raise AssertionError("Output differs for {}\n{}\n{}".format(suffix, a, b))


def main():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        plain = root / "input.vcf"
        ordinary = root / "ordinary.vcf.gz"
        bgzf = root / "blocked.vcf.gz"
        make_vcf(plain)
        with plain.open("rb") as source, gzip.open(ordinary, "wb") as target:
            shutil.copyfileobj(source, target)
        make_bgzf(plain, bgzf)

        fakebin = root / "fakebin"
        make_fake_tools(fakebin)
        env = os.environ.copy()
        env["PATH"] = str(fakebin) + os.pathsep + env.get("PATH", "")

        common = ["-m", "3", "-f", "-n", "-b", "-w", "-o", "S4", "-t", "4"]
        out_plain = root / "plain"
        stdout = run(["-i", str(plain), "--output-folder", str(out_plain)] + common, env)
        assert "Detected input format: PLAIN" in stdout
        assert "plain VCF direct parallel byte ranges" in stdout

        out_plain_stream = root / "plain_stream"
        stdout = run(["-i", str(plain), "--input-backend", "plain-stream",
                      "--output-folder", str(out_plain_stream)] + common, env)
        assert "plain VCF sequential binary stream" in stdout
        compare_outputs(out_plain, out_plain_stream)

        # Confirm that raw byte boundaries inside a very long VCF row neither
        # duplicate nor omit the row when independent workers align to lines.
        boundary = root / "boundary.vcf"
        make_boundary_vcf(boundary)
        boundary_direct = root / "boundary_direct"
        stdout = run(["-i", str(boundary),
                      "--output-folder", str(boundary_direct)] + common, env)
        assert "Direct byte ranges: 2" in stdout
        boundary_stream = root / "boundary_stream"
        run(["-i", str(boundary), "--input-backend", "plain-stream",
             "--output-folder", str(boundary_stream)] + common, env)
        compare_outputs(boundary_direct, boundary_stream)

        out_gzip = root / "gzip"
        stdout = run(["-i", str(ordinary), "--output-folder", str(out_gzip)] + common, env)
        assert "Detected input format: GZIP" in stdout
        assert "Python stdlib gzip" in stdout
        compare_outputs(out_plain, out_gzip)

        out_bgzf = root / "bgzf"
        stdout = run(["-i", str(bgzf), "--output-folder", str(out_bgzf)] + common, env)
        assert "Detected input format: BGZF" in stdout
        assert "HTSlib bgzip multithreaded" in stdout
        compare_outputs(out_plain, out_bgzf)

        (Path(str(bgzf) + ".tbi")).touch()
        out_indexed = root / "indexed"
        stdout = run(["-i", str(bgzf), "--output-folder", str(out_indexed)] + common, env)
        assert "indexed BGZF parallel regions" in stdout
        assert "Indexed regions:" in stdout
        assert "Total of genotypes processed: 8" in stdout
        assert "excluded for being MNPs: 2" in stdout
        compare_outputs(out_plain, out_indexed)

        out_no_index = root / "no_index"
        stdout = run(["-i", str(bgzf), "--no-indexed-regions",
                      "--output-folder", str(out_no_index)] + common, env)
        assert "HTSlib bgzip multithreaded" in stdout
        compare_outputs(out_plain, out_no_index)

        # A stale/corrupt index must not break auto mode: tabix failure falls
        # back to the BGZF streaming backend and preserves the matrices.
        badbin = root / "badbin"
        shutil.copytree(fakebin, badbin)
        (badbin / "tabix").write_text("""#!/usr/bin/env python3
import gzip, sys
args = sys.argv[1:]
if args[0] == '-l':
    sys.stdout.write('chr1\nchr2\n')
    raise SystemExit(0)
sys.stderr.write('simulated stale index\n')
raise SystemExit(2)
""", encoding="utf-8")
        (badbin / "tabix").chmod(0o755)
        bad_env = env.copy()
        bad_env["PATH"] = str(badbin) + os.pathsep + os.environ.get("PATH", "")
        out_fallback = root / "fallback"
        stdout = run(["-i", str(bgzf), "--output-folder", str(out_fallback)] + common, bad_env)
        assert "indexed BGZF mode failed" in stdout
        assert "Fallback backend: HTSlib bgzip" in stdout
        compare_outputs(out_plain, out_fallback)

    print("All backend tests passed")


if __name__ == "__main__":
    main()
