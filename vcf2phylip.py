#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert SNPs in VCF format into PHYLIP, FASTA, NEXUS, or binary NEXUS
matrices for phylogenetic analysis.

This multithreaded/multiprocess edition preserves the command-line behavior
of vcf2phylip v2.9 and adds ``-t/--threads``. VCF chunks are parsed with
multiple worker processes (to bypass Python's GIL), while final sample
sequences are assembled concurrently with threads. The temporary matrix is
stored in transposed blocks, avoiding the original full temporary-file scan
for every sample.

Any ploidy is allowed, but binary NEXUS is produced only for diploid VCFs.
"""

__author__ = "Edgardo M. Ortiz"
__credits__ = "Juan D. Palacio-Mejia"
__modifier__ = "Ma Wenxin"
__version__ = "2.9-mt5"
__email__ = "e.ortiz.v@gmail.com"
__date__ = "2026-07-22"

import argparse
import gzip
import os
import random
import sys
import zlib
from collections import deque
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

# Automatic chunking targets enough tasks to keep all workers busy while
# avoiding thousands of tiny process submissions. Limits are intentionally
# broad because this edition is tuned for 300-400 samples and 1M-10M SNPs.
MIB = 1024 * 1024
MIN_CHUNK_TARGET_BYTES = 4 * MIB
MAX_CHUNK_TARGET_BYTES = 128 * MIB
MIN_CHUNK_RECORDS = 2000
MAX_CHUNK_RECORDS = 100000
PROFILE_SAMPLE_BYTES = 32 * MIB
PROFILE_SAMPLE_RECORDS = 20000

# Dictionary of IUPAC ambiguities for nucleotides.
# '*' is a deletion in GATK; deletions are ignored in consensus. Lowercase
# consensus is used when an 'N' or '*' is part of the genotype.
AMBIG = {
    "A": "A", "C": "C", "G": "G", "N": "N", "T": "T",
    "*A": "a", "*C": "c", "*G": "g", "*N": "n", "*T": "t",
    "AC": "M", "AG": "R", "AN": "a", "AT": "W", "CG": "S",
    "CN": "c", "CT": "Y", "GN": "g", "GT": "K", "NT": "t",
    "*AC": "m", "*AG": "r", "*AN": "a", "*AT": "w", "*CG": "s",
    "*CN": "c", "*CT": "y", "*GN": "g", "*GT": "k", "*NT": "t",
    "ACG": "V", "ACN": "m", "ACT": "H", "AGN": "r", "AGT": "D",
    "ANT": "w", "CGN": "s", "CGT": "B", "CNT": "y", "GNT": "k",
    "*ACG": "v", "*ACN": "m", "*ACT": "h", "*AGN": "r", "*AGT": "d",
    "*ANT": "w", "*CGN": "s", "*CGT": "b", "*CNT": "y", "*GNT": "k",
    "ACGN": "v", "ACGT": "N", "ACNT": "h", "AGNT": "d", "CGNT": "b",
    "*ACGN": "v", "*ACGT": "N", "*ACNT": "h", "*AGNT": "d", "*CGNT": "b",
    "*": "-", "*ACGNT": "N",
}

# Dictionary for translating biallelic SNPs into SNAPP states.
GEN_BIN = {
    "./.": "?", ".|.": "?",
    "0/0": "0", "0|0": "0",
    "0/1": "1", "0|1": "1", "1/0": "1", "1|0": "1",
    "1/1": "2", "1|1": "2",
}


def extract_sample_names(vcf_file):
    """Extract sample names from the #CHROM line of a VCF file."""
    opener = gzip.open if vcf_file.lower().endswith(".gz") else open
    sample_names = []
    with opener(vcf_file, "rt") as vcf:
        for line in vcf:
            line = line.strip("\n")
            if line.startswith("#CHROM"):
                record = line.split("\t")
                sample_names = [record[i].replace("./", "") for i in range(9, len(record))]
                break
    return sample_names


def is_anomalous(record, num_samples):
    """Return True when a VCF row has an unexpected number of columns."""
    return len(record) != num_samples + 9


def is_snp(record):
    """Return True for single-nucleotide REF/ALT alleles, including multiallelic SNPs."""
    alt = record[4].replace("<NON_REF>", record[3])
    return len(record[3]) == 1 and len(alt) - alt.count(",") == alt.count(",") + 1


def num_genotypes(record, num_samples):
    """Count samples whose genotype field does not start with a missing allele."""
    missing = sum(1 for field in record[9:num_samples + 9] if field.startswith("."))
    return num_samples - missing


def get_matrix_column(record, num_samples, resolve_iupac, rng=None):
    """Transform one VCF record into one nucleotide matrix column."""
    nt_dict = {"0": record[3].replace("-", "*").upper(), ".": "N"}
    alt = record[4].replace("-", "*").replace("<NON_REF>", nt_dict["0"])
    for n, allele in enumerate(alt.split(","), start=1):
        nt_dict[str(n)] = allele

    chars = []
    choose = random.choice if rng is None else rng.choice
    for i in range(9, num_samples + 9):
        geno_num = record[i].split(":", 1)[0].replace("/", "").replace("|", "")
        try:
            geno_nuc = "".join(sorted({nt_dict[j] for j in geno_num}))
        except KeyError:
            return "malformed"
        if resolve_iupac:
            chars.append(AMBIG[nt_dict[choose(geno_num)]])
        else:
            chars.append(AMBIG[geno_nuc])
    return "".join(chars)


def get_matrix_column_bin(record, num_samples):
    """Return one binary NEXUS column; unsupported genotypes become '?'."""
    chars = []
    for i in range(9, num_samples + 9):
        genotype = record[i].split(":", 1)[0]
        chars.append(GEN_BIN.get(genotype, "?"))
    return "".join(chars)


def iter_vcf_chunks(vcf, target_bytes, target_records):
    """Yield bounded VCF chunks using both text size and data-record limits."""
    chunk = []
    size = 0
    records = 0
    for line in vcf:
        chunk.append(line)
        size += len(line)
        if line and not line.startswith("#"):
            records += 1
        if size >= target_bytes or records >= target_records:
            yield chunk
            chunk = []
            size = 0
            records = 0
    if chunk:
        yield chunk


def process_vcf_chunk(lines, num_samples, min_samples_locus, need_nt,
                      need_bin, resolve_iupac, write_used, random_seed=None):
    """Parse one VCF chunk. This function is safe to run in a worker process."""
    rng = random.Random(random_seed) if random_seed is not None else None

    nt_rows = []
    bin_rows = []
    used_rows = []
    malformed_lines = []

    snp_num = 0
    snp_accepted = 0
    snp_shallow = 0
    mnp_num = 0
    snp_biallelic = 0

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        record = line.split("\t")
        snp_num += 1

        if is_anomalous(record, num_samples):
            malformed_lines.append(line)
            continue

        num_samples_locus = num_genotypes(record, num_samples)
        if num_samples_locus < min_samples_locus:
            snp_shallow += 1
            continue

        if not is_snp(record):
            mnp_num += 1
            continue

        # Preserve v2.9 behavior: if a nucleotide matrix is requested and its
        # genotype cannot be decoded, skip this site for all requested outputs.
        if need_nt:
            site_tmp = get_matrix_column(record, num_samples, resolve_iupac, rng)
            if site_tmp == "malformed":
                malformed_lines.append(line)
                continue
            snp_accepted += 1
            nt_rows.append(site_tmp)
            if write_used:
                used_rows.append((record[0], record[1], num_samples_locus))

        if need_bin and len(record[4]) == 1:
            snp_biallelic += 1
            bin_rows.append(get_matrix_column_bin(record, num_samples))

    counters = (snp_num, snp_accepted, snp_shallow, mnp_num, snp_biallelic)
    return nt_rows, bin_rows, used_rows, malformed_lines, counters


def ordered_parallel_results(chunks, workers, task_args, resolve_iupac):
    """Run chunk tasks in parallel while yielding results in input order."""
    # Keep workers busy without retaining two full waves of large chunks/results.
    max_pending = max(1, workers + 2)
    pending = deque()

    with ProcessPoolExecutor(max_workers=workers) as executor:
        for chunk in chunks:
            seed = random.getrandbits(64) if resolve_iupac else None
            future = executor.submit(process_vcf_chunk, chunk, *task_args, seed)
            pending.append(future)
            if len(pending) >= max_pending:
                yield pending.popleft().result()

        while pending:
            yield pending.popleft().result()


def write_transposed_block(handle, rows, num_samples):
    """
    Store one row-major matrix chunk as a sample-major block.

    Returns (block_start_offset, sites_in_block). The final sequence for sample
    s can later be read from offset + s * sites_in_block.
    """
    if not rows:
        return None

    width = len(rows)
    for row in rows:
        if len(row) != num_samples:
            raise ValueError("Internal error: matrix row width does not match sample count")

    start = handle.tell()
    for sample_chars in zip(*rows):
        handle.write("".join(sample_chars).encode("ascii"))
    return start, width


def assemble_sequence(temp_path, blocks, sample_index):
    """Assemble one sample sequence from transposed temporary blocks."""
    if temp_path is None or not blocks:
        return ""

    sequence = bytearray()
    with open(temp_path, "rb") as handle:
        for start, width in blocks:
            handle.seek(start + sample_index * width)
            data = handle.read(width)
            if len(data) != width:
                raise OSError("Temporary matrix ended unexpectedly")
            sequence.extend(data)
    return sequence.decode("ascii")


def assemble_sample(sample_index, nt_temp_path, nt_blocks, bin_temp_path, bin_blocks):
    """Assemble nucleotide and/or binary sequences for one sample."""
    nt_seq = assemble_sequence(nt_temp_path, nt_blocks, sample_index)
    bin_seq = assemble_sequence(bin_temp_path, bin_blocks, sample_index)
    return sample_index, nt_seq, bin_seq


def detect_available_workers():
    """Return 100% of the CPU count available to this job and its source.

    Input size and sample count never reduce this value. They only affect
    chunk sizing, so omitting ``-t`` always requests the maximum CPU
    parallelism allowed by the scheduler, CPU affinity, or operating system.
    """
    scheduler_vars = (
        "SLURM_CPUS_PER_TASK",
        "PBS_NP",
        "NSLOTS",
        "LSB_DJOB_NUMPROC",
    )
    for variable in scheduler_vars:
        value = os.environ.get(variable)
        if value:
            try:
                workers = int(value)
            except ValueError:
                continue
            if workers > 0:
                return workers, variable

    try:
        affinity_count = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        affinity_count = 0
    if affinity_count > 0:
        return affinity_count, "CPU affinity"

    return max(1, os.cpu_count() or 1), "os.cpu_count()"


def clamp(value, lower, upper):
    """Clamp a numeric value to an inclusive interval."""
    return max(lower, min(upper, value))


def format_binary_size(size_bytes):
    """Return a compact IEC size string."""
    size = float(max(0, size_bytes))
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024.0 or unit == "TiB":
            return "{:.1f} {}".format(size, unit)
        size /= 1024.0
    return "{:.1f} TiB".format(size)


def sample_gzip_prefix(vcf_path, output_limit=PROFILE_SAMPLE_BYTES):
    """
    Decompress a bounded gzip/BGZF prefix and count compressed bytes consumed.

    Using zlib directly avoids gzip.GzipFile read-ahead, which can otherwise
    make compression-ratio estimates inaccurate for highly compressible VCFs.
    Concatenated gzip members (including BGZF blocks) are handled explicitly.
    """
    output = bytearray()
    compressed_consumed = 0
    pending = b""
    decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
    raw_eof = False

    with open(vcf_path, "rb") as raw:
        while len(output) < output_limit:
            if not pending:
                pending = raw.read(64 * 1024)
                if not pending:
                    raw_eof = True
                    break

            remaining = output_limit - len(output)
            produced = decompressor.decompress(pending, remaining)
            output.extend(produced)

            unconsumed = decompressor.unconsumed_tail
            unused = decompressor.unused_data
            consumed_now = len(pending) - len(unconsumed) - len(unused)
            compressed_consumed += max(0, consumed_now)

            if len(output) >= output_limit:
                break

            if unconsumed:
                pending = unconsumed
                continue

            if decompressor.eof:
                pending = unused
                decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
                continue

            pending = b""

    reached_eof = raw_eof and not pending
    return bytes(output), max(1, compressed_consumed), reached_eof


def summarize_sample(sample_bytes, reached_eof):
    """Summarize complete VCF data rows in a sampled byte prefix."""
    lines = sample_bytes.splitlines(keepends=True)
    if lines and not reached_eof and not lines[-1].endswith((b"\n", b"\r")):
        lines.pop()

    sampled_record_bytes = 0
    sampled_records = 0
    for line in lines:
        if line and not line.startswith(b"#"):
            sampled_records += 1
            sampled_record_bytes += len(line)
            if sampled_records >= PROFILE_SAMPLE_RECORDS:
                break
    return sampled_record_bytes, sampled_records


def profile_vcf(vcf_path):
    """
    Estimate uncompressed VCF size, record count, and average record width.

    Plain VCF size is exact. For gzip/BGZF input, a bounded prefix is
    decompressed and its observed compression ratio is extrapolated to the
    whole file. The file is reopened for the real conversion, so profiling
    does not alter output or processing order.
    """
    path = Path(vcf_path)
    stored_bytes = path.stat().st_size
    is_compressed = str(path).lower().endswith(".gz")

    if is_compressed:
        sample, compressed_consumed, reached_eof = sample_gzip_prefix(path)
        sampled_uncompressed = len(sample)
        sampled_record_bytes, sampled_records = summarize_sample(
            sample, reached_eof)

        if reached_eof:
            estimated_uncompressed = sampled_uncompressed
            size_source = "exact gzip scan"
        else:
            observed_ratio = sampled_uncompressed / compressed_consumed
            observed_ratio = clamp(observed_ratio, 1.0, 1000.0)
            estimated_uncompressed = int(stored_bytes * observed_ratio)
            estimated_uncompressed = max(estimated_uncompressed,
                                         sampled_uncompressed)
            size_source = "gzip/BGZF prefix estimate"
    else:
        sampled = bytearray()
        with open(path, "rb") as raw:
            while len(sampled) < PROFILE_SAMPLE_BYTES:
                block = raw.read(min(1024 * 1024,
                                     PROFILE_SAMPLE_BYTES - len(sampled)))
                if not block:
                    break
                sampled.extend(block)
        reached_eof = len(sampled) >= stored_bytes
        sampled_uncompressed = len(sampled)
        sampled_record_bytes, sampled_records = summarize_sample(
            bytes(sampled), reached_eof)
        estimated_uncompressed = stored_bytes
        size_source = "exact file size"

    if sampled_records > 0:
        average_record_bytes = sampled_record_bytes / sampled_records
        estimated_records = int(round(estimated_uncompressed /
                                      average_record_bytes))
        estimated_records = max(sampled_records, estimated_records)
    else:
        average_record_bytes = 0.0
        estimated_records = 0

    return {
        "stored_bytes": stored_bytes,
        "estimated_uncompressed_bytes": estimated_uncompressed,
        "estimated_records": estimated_records,
        "average_record_bytes": average_record_bytes,
        "is_compressed": is_compressed,
        "size_source": size_source,
    }


def choose_chunk_settings(workers, vcf_path, chunk_size_mb=None):
    """
    Choose chunk limits from CPU count and estimated VCF workload.

    The target is several ordered tasks per worker: enough for load balancing,
    but not so many that process scheduling and pickling dominate runtime.
    """
    profile = profile_vcf(vcf_path)
    estimated_bytes = max(1, profile["estimated_uncompressed_bytes"])
    estimated_records = profile["estimated_records"]

    if estimated_records <= 1000000:
        tasks_per_worker = 8
    elif estimated_records <= 5000000:
        tasks_per_worker = 10
    else:
        tasks_per_worker = 12

    # gzip decompression is sequential in the reader process, so a slightly
    # smaller task count reduces IPC overhead without starving parser workers.
    if profile["is_compressed"]:
        tasks_per_worker = max(6, tasks_per_worker - 2)

    target_tasks = max(1, workers * tasks_per_worker)
    target_tasks = min(target_tasks, 4096)

    if chunk_size_mb is not None:
        target_bytes = chunk_size_mb * MIB
        source = "user byte limit; record limit estimated from VCF profile"
        if profile["average_record_bytes"] > 0:
            target_records = int(round(
                target_bytes / profile["average_record_bytes"]))
        elif estimated_records > 0:
            target_records = int(round(estimated_records / target_tasks))
        else:
            target_records = MAX_CHUNK_RECORDS
    else:
        target_bytes = int(round(estimated_bytes / target_tasks))
        if estimated_records > 0:
            target_records = int(round(estimated_records / target_tasks))
        elif profile["average_record_bytes"] > 0:
            target_records = int(round(
                target_bytes / profile["average_record_bytes"]))
        else:
            target_records = MAX_CHUNK_RECORDS
        source = "auto from CPU count and estimated VCF workload"

    if chunk_size_mb is None:
        target_bytes = int(clamp(target_bytes,
                                 MIN_CHUNK_TARGET_BYTES,
                                 MAX_CHUNK_TARGET_BYTES))
        # Round auto-selected values to whole MiB for stable, readable settings.
        target_bytes = max(MIN_CHUNK_TARGET_BYTES,
                           int(round(target_bytes / MIB)) * MIB)
    else:
        # A manual override is explicit and is therefore honored exactly.
        target_bytes = int(target_bytes)
    target_records = int(clamp(target_records,
                               MIN_CHUNK_RECORDS,
                               MAX_CHUNK_RECORDS))
    # Round record limits to reduce noisy differences from prefix sampling.
    record_quantum = 500 if target_records < 10000 else 1000
    target_records = max(MIN_CHUNK_RECORDS,
                         int(round(target_records / record_quantum)) * record_quantum)

    return target_bytes, target_records, source, target_tasks, profile


def positive_int(value):
    """argparse validator for positive integers."""
    try:
        integer = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if integer < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return integer


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", dest="filename", required=True,
                        help="Name of the input VCF file, can be gzipped")
    parser.add_argument("--output-folder", dest="folder", default="./",
                        help="Output folder name; created if it does not exist (default: ./)")
    parser.add_argument("--output-prefix", dest="prefix",
                        help="Prefix for output filenames (input VCF name by default)")
    parser.add_argument("-m", "--min-samples-locus", dest="min_samples_locus",
                        type=int, default=4,
                        help="Minimum samples required at a locus (default=4)")
    parser.add_argument("-o", "--outgroup", dest="outgroup", default="",
                        help="Outgroup name; its sequence is written first")
    parser.add_argument("-p", "--phylip-disable", action="store_true",
                        dest="phylipdisable",
                        help="Disable PHYLIP output, which is enabled by default")
    parser.add_argument("-f", "--fasta", action="store_true",
                        help="Write a FASTA matrix (disabled by default)")
    parser.add_argument("-n", "--nexus", action="store_true",
                        help="Write a NEXUS matrix (disabled by default)")
    parser.add_argument("-b", "--nexus-binary", action="store_true", dest="nexusbin",
                        help="Write binary NEXUS for biallelic diploid SNPs")
    parser.add_argument("-r", "--resolve-IUPAC", action="store_true",
                        dest="resolve_iupac",
                        help="Randomly resolve heterozygous genotypes")
    parser.add_argument("-w", "--write-used-sites", action="store_true", dest="write_used",
                        help="Save coordinates that passed filters and were used")
    parser.add_argument("-t", "--threads", type=positive_int, default=None,
                        help="Parallel workers; default: 100%% of CPUs available to the job; "
                             "VCF size never reduces this count")
    parser.add_argument("--chunk-size-mb", type=positive_int, default=None,
                        help="Override auto byte limit in MiB; record limit remains adaptive")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {version}".format(version=__version__))
    args = parser.parse_args()

    if args.threads is None:
        args.threads, worker_source = detect_available_workers()
        worker_description = "100% auto-detected from {}; not reduced by VCF size".format(
            worker_source)
    else:
        worker_description = "user specified"

    outgroup = args.outgroup.split(",")[0].split(";")[0]
    need_nt = args.fasta or args.nexus or not args.phylipdisable
    need_bin = args.nexusbin

    input_path = Path(args.filename)
    if not input_path.exists():
        print("\nInput VCF file not found, please verify the provided path")
        sys.exit(1)

    chunk_bytes, chunk_records, chunk_source, target_tasks, vcf_profile = \
        choose_chunk_settings(args.threads, args.filename, args.chunk_size_mb)

    sample_names = extract_sample_names(args.filename)
    num_samples = len(sample_names)
    if num_samples == 0:
        print("\nSample names not found in VCF, your file may be corrupt or missing the header.\n")
        sys.exit(1)

    print("\nConverting file '{}':\n".format(args.filename))
    print("Number of samples in VCF: {:d}".format(num_samples))
    print("Parallel workers: {:d} ({})".format(args.threads, worker_description))
    if vcf_profile["is_compressed"]:
        print("VCF size: {} stored; about {} uncompressed ({})".format(
            format_binary_size(vcf_profile["stored_bytes"]),
            format_binary_size(vcf_profile["estimated_uncompressed_bytes"]),
            vcf_profile["size_source"]))
    else:
        print("VCF size: {} ({})".format(
            format_binary_size(vcf_profile["estimated_uncompressed_bytes"]),
            vcf_profile["size_source"]))
    if vcf_profile["estimated_records"] > 0:
        print("Estimated data records: {:,}; average sampled row: {:.0f} bytes".format(
            vcf_profile["estimated_records"],
            vcf_profile["average_record_bytes"]))
    print("VCF chunk limit: {:.0f} MiB or {:,d} data records; target about {:,d} tasks ({})".format(
        chunk_bytes / MIB, chunk_records, target_tasks, chunk_source))

    args.min_samples_locus = min(num_samples, args.min_samples_locus)

    if not args.prefix:
        parts = input_path.name.split(".")
        prefix_parts = []
        for part in parts:
            if part.lower() == "vcf":
                break
            prefix_parts.append(part)
        args.prefix = ".".join(prefix_parts)
    args.prefix += ".min" + str(args.min_samples_locus)

    output_folder = Path(args.folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    outfile = str(output_folder / args.prefix)

    nt_temp_path = outfile + ".tmp" if need_nt else None
    bin_temp_path = outfile + ".bin.tmp" if need_bin else None
    nt_blocks = []
    bin_blocks = []

    used_sites = None
    nt_temp = None
    bin_temp = None

    try:
        if args.write_used:
            used_sites = open(outfile + ".used_sites.tsv", "w", encoding="utf-8")
            used_sites.write("#CHROM\tPOS\tNUM_SAMPLES\n")
        if need_nt:
            nt_temp = open(nt_temp_path, "wb")
        if need_bin:
            bin_temp = open(bin_temp_path, "wb")

        opener = gzip.open if args.filename.lower().endswith(".gz") else open

        snp_num = 0
        snp_accepted = 0
        snp_shallow = 0
        mnp_num = 0
        snp_biallelic = 0
        next_progress = 500000

        task_args = (
            num_samples,
            args.min_samples_locus,
            need_nt,
            need_bin,
            args.resolve_iupac,
            args.write_used,
        )

        with opener(args.filename, "rt") as vcf:
            chunks = iter_vcf_chunks(vcf, chunk_bytes, chunk_records)
            if args.threads == 1:
                results = (
                    process_vcf_chunk(chunk, *task_args, None)
                    for chunk in chunks
                )
            else:
                results = ordered_parallel_results(
                    chunks, args.threads, task_args, args.resolve_iupac)

            for nt_rows, bin_rows, used_rows, malformed_lines, counters in results:
                for malformed in malformed_lines:
                    print("Skipping malformed line:\n{}".format(malformed))

                c_total, c_accepted, c_shallow, c_mnp, c_biallelic = counters
                snp_num += c_total
                snp_accepted += c_accepted
                snp_shallow += c_shallow
                mnp_num += c_mnp
                snp_biallelic += c_biallelic

                while snp_num >= next_progress:
                    print("{:d} genotypes processed.".format(next_progress))
                    next_progress += 500000

                if nt_rows:
                    block = write_transposed_block(nt_temp, nt_rows, num_samples)
                    nt_blocks.append(block)
                if bin_rows:
                    block = write_transposed_block(bin_temp, bin_rows, num_samples)
                    bin_blocks.append(block)
                if used_sites is not None:
                    for chrom, pos, present in used_rows:
                        used_sites.write("{}\t{}\t{}\n".format(chrom, pos, present))

        if nt_temp is not None:
            nt_temp.close()
            nt_temp = None
        if bin_temp is not None:
            bin_temp.close()
            bin_temp = None
        if used_sites is not None:
            used_sites.close()
            used_sites = None

        print("Total of genotypes processed: {:d}".format(snp_num))
        print("Genotypes excluded because they exceeded the amount of missing data allowed: {:d}".format(snp_shallow))
        print("Genotypes that passed missing data filter but were excluded for being MNPs: {:d}".format(mnp_num))
        print("SNPs that passed the filters: {:d}".format(snp_accepted))
        if args.nexusbin:
            print("Biallelic SNPs selected for binary NEXUS: {:d}".format(snp_biallelic))
        if args.write_used:
            print("Used sites saved to: '" + outfile + ".used_sites.tsv'")
        print("")

        output_phy = open(outfile + ".phy", "w", encoding="utf-8") if not args.phylipdisable else None
        output_fas = open(outfile + ".fasta", "w", encoding="utf-8") if args.fasta else None
        output_nex = open(outfile + ".nexus", "w", encoding="utf-8") if args.nexus else None
        output_nexbin = open(outfile + ".bin.nexus", "w", encoding="utf-8") if args.nexusbin else None

        try:
            if output_phy is not None:
                output_phy.write("{:d} {:d}\n".format(num_samples, snp_accepted))
            if output_nex is not None:
                output_nex.write(
                    "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={:d} NCHAR={:d};\n"
                    "\tFORMAT DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n".format(
                        num_samples, snp_accepted))
            if output_nexbin is not None:
                output_nexbin.write(
                    "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={:d} NCHAR={:d};\n"
                    "\tFORMAT DATATYPE=SNP MISSING=? GAP=- ;\nMATRIX\n".format(
                        num_samples, snp_biallelic))

            len_longest_name = max(len(name) for name in sample_names)
            idx_outgroup = sample_names.index(outgroup) if outgroup in sample_names else None
            sample_order = []
            if idx_outgroup is not None:
                sample_order.append(idx_outgroup)
            sample_order.extend(i for i in range(num_samples) if i != idx_outgroup)

            # Keep the configured CPU ceiling unchanged during output assembly.
            # ThreadPoolExecutor creates threads lazily, so fewer samples/tasks may
            # naturally leave some CPUs idle, but the program does not downscale the
            # requested worker count based on sample number or VCF size.
            worker_count = args.threads
            output_executor = (
                ThreadPoolExecutor(max_workers=worker_count)
                if worker_count > 1 else None
            )
            try:
                for group_start in range(0, len(sample_order), worker_count):
                    group = sample_order[group_start:group_start + worker_count]
                    if output_executor is None:
                        assembled = [
                            assemble_sample(group[0], nt_temp_path, nt_blocks,
                                            bin_temp_path, bin_blocks)
                        ]
                    else:
                        assembled = list(output_executor.map(
                            lambda s: assemble_sample(s, nt_temp_path, nt_blocks,
                                                      bin_temp_path, bin_blocks),
                            group))

                    for sample_index, nt_seq, bin_seq in assembled:
                        name = sample_names[sample_index]
                        padding = (len_longest_name + 3 - len(name)) * " "

                        if output_fas is not None:
                            output_fas.write(">{}\n{}\n".format(name, nt_seq))
                        if output_phy is not None:
                            output_phy.write(name + padding + nt_seq + "\n")
                        if output_nex is not None:
                            output_nex.write(name + padding + nt_seq + "\n")
                        if output_nexbin is not None:
                            output_nexbin.write(name + padding + bin_seq + "\n")

                        if sample_index == idx_outgroup:
                            if need_nt:
                                print("Outgroup, '{}', added to the matrix(ces).".format(outgroup))
                            if need_bin:
                                print("Outgroup, '{}', added to the binary matrix.".format(outgroup))
                        else:
                            if need_nt:
                                print("Sample {:d} of {:d}, '{}', added to the nucleotide matrix(ces).".format(
                                    sample_index + 1, num_samples, name))
                            if need_bin:
                                print("Sample {:d} of {:d}, '{}', added to the binary matrix.".format(
                                    sample_index + 1, num_samples, name))

            finally:
                if output_executor is not None:
                    output_executor.shutdown(wait=True)

            print()
            if output_nex is not None:
                output_nex.write(";\nEND;\n")
            if output_nexbin is not None:
                output_nexbin.write(";\nEND;\n")
        finally:
            if output_phy is not None:
                output_phy.close()
            if output_fas is not None:
                output_fas.close()
            if output_nex is not None:
                output_nex.close()
            if output_nexbin is not None:
                output_nexbin.close()

        if not args.phylipdisable:
            print("PHYLIP matrix saved to: " + outfile + ".phy")
        if args.fasta:
            print("FASTA matrix saved to: " + outfile + ".fasta")
        if args.nexus:
            print("NEXUS matrix saved to: " + outfile + ".nexus")
        if args.nexusbin:
            print("BINARY NEXUS matrix saved to: " + outfile + ".bin.nexus")

    finally:
        if nt_temp is not None:
            nt_temp.close()
        if bin_temp is not None:
            bin_temp.close()
        if used_sites is not None:
            used_sites.close()
        if nt_temp_path is not None:
            Path(nt_temp_path).unlink(missing_ok=True)
        if bin_temp_path is not None:
            Path(bin_temp_path).unlink(missing_ok=True)

    print("\nDone!\n")


if __name__ == "__main__":
    main()
