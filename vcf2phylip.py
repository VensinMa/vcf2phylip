#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert SNPs in VCF format into PHYLIP, FASTA, NEXUS, or binary NEXUS
matrices for phylogenetic analysis.

This multithreaded/multiprocess edition preserves the command-line behavior
of vcf2phylip v2.9. It auto-detects plain VCF, ordinary gzip, and real BGZF
from file bytes; selects direct byte-range, ISA-L, zlib-ng, bgzip, or tabix
backends; and parses genotypes with multiple worker processes. Plain VCFs can
be read directly by independent ordered byte ranges, while indexed BGZF VCFs
can be decompressed and parsed concurrently by ordered genomic regions. Final
sample sequences are assembled concurrently from transposed temporary blocks.

Any ploidy is allowed, but binary NEXUS is produced only for diploid VCFs.
"""

__author__ = "Edgardo M. Ortiz"
__credits__ = "Juan D. Palacio-Mejia"
__modifier__ = "Ma Wenxin"
__version__ = "2.9-mt7"
__email__ = "e.ortiz.v@gmail.com"
__date__ = "2026-07-23"

import argparse
import gzip
import importlib.util
import math
import os
import random
import shutil
import struct
import subprocess
import sys
import tempfile
import zlib
from concurrent.futures import (FIRST_COMPLETED, ProcessPoolExecutor,
                                ThreadPoolExecutor, wait)
from contextlib import contextmanager
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
MIN_PLAIN_RANGE_BYTES = 4 * MIB

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

AMBIG_BYTES = {
    key.encode("ascii"): value.encode("ascii") for key, value in AMBIG.items()
}
GEN_BIN_BYTES = {
    key.encode("ascii"): value.encode("ascii") for key, value in GEN_BIN.items()
}


def inspect_input_format(vcf_file):
    """Return ``plain``, ``gzip``, or ``bgzf`` from file bytes, not suffix."""
    path = Path(vcf_file)
    with path.open("rb") as handle:
        fixed = handle.read(12)
        if len(fixed) < 3 or fixed[:3] != b"\x1f\x8b\x08":
            return "plain"
        if len(fixed) < 12 or not (fixed[3] & 0x04):
            return "gzip"

        xlen = struct.unpack("<H", fixed[10:12])[0]
        extra = handle.read(xlen)
        offset = 0
        block_size = None
        while offset + 4 <= len(extra):
            si1, si2 = extra[offset:offset + 2]
            slen = struct.unpack("<H", extra[offset + 2:offset + 4])[0]
            value_start = offset + 4
            value_end = value_start + slen
            if value_end > len(extra):
                break
            if si1 == ord("B") and si2 == ord("C") and slen == 2:
                block_size = struct.unpack("<H", extra[value_start:value_end])[0] + 1
                break
            offset = value_end

        if block_size is None or not (18 <= block_size <= 65536):
            return "gzip"

        handle.seek(0)
        first_block = handle.read(block_size)
        if len(first_block) != block_size:
            return "gzip"
        try:
            zlib.decompress(first_block, 16 + zlib.MAX_WBITS)
        except zlib.error:
            return "gzip"
        return "bgzf"


def find_vcf_index(vcf_file):
    """Return a standard TBI/CSI sidecar path when present."""
    path = Path(vcf_file)
    for suffix in (".csi", ".tbi"):
        candidate = Path(str(path) + suffix)
        if candidate.is_file():
            return candidate
    return None


def read_vcf_header(vcf_file, input_format=None):
    """Read VCF header, sample names, and declared contig lengths."""
    input_format = input_format or inspect_input_format(vcf_file)
    opener = open if input_format == "plain" else gzip.open
    sample_names = []
    contig_lengths = {}
    header_lines = []

    with opener(vcf_file, "rb") as vcf:
        for raw_line in vcf:
            if not raw_line.startswith(b"#"):
                break
            header_lines.append(raw_line)
            line = raw_line.rstrip(b"\r\n")
            if line.startswith(b"##contig=<") and line.endswith(b">"):
                inner = line[len(b"##contig=<"):-1]
                fields = {}
                for item in inner.split(b","):
                    if b"=" in item:
                        key, value = item.split(b"=", 1)
                        fields[key.strip().lower()] = value.strip()
                contig_id = fields.get(b"id")
                length = fields.get(b"length")
                if contig_id and length:
                    try:
                        contig_lengths[contig_id.decode("utf-8")] = int(length)
                    except (UnicodeDecodeError, ValueError):
                        pass
            if line.startswith(b"#CHROM"):
                record = line.split(b"\t")
                sample_names = [
                    field.decode("utf-8").replace("./", "")
                    for field in record[9:]
                ]
                break

    return sample_names, contig_lengths, header_lines


def find_plain_data_start(vcf_file):
    """Return the byte offset immediately after the plain-VCF header."""
    with open(vcf_file, "rb") as handle:
        while True:
            line = handle.readline()
            if not line:
                return handle.tell()
            if line.startswith(b"#CHROM"):
                return handle.tell()
            if not line.startswith(b"#"):
                raise ValueError(
                    "VCF data started before a #CHROM header line was found")


def build_plain_ranges(vcf_file, data_start, target_tasks, max_range_bytes):
    """Build ordered raw byte ranges for independent plain-VCF readers."""
    file_size = Path(vcf_file).stat().st_size
    data_bytes = max(0, file_size - data_start)
    if data_bytes == 0:
        return []

    tasks_for_size = max(1, int(math.ceil(
        data_bytes / float(max(1, max_range_bytes)))))
    max_useful_tasks = max(1, int(math.ceil(
        data_bytes / float(MIN_PLAIN_RANGE_BYTES))))
    desired_tasks = max(int(target_tasks), tasks_for_size)
    task_count = max(1, min(desired_tasks, max_useful_tasks, 4096))
    ranges = []
    for order in range(task_count):
        start = data_start + (data_bytes * order) // task_count
        end = data_start + (data_bytes * (order + 1)) // task_count
        ranges.append({
            "order": order,
            "start": start,
            "end": end,
            "data_start": data_start,
            "is_last": order == task_count - 1,
        })
    return ranges


def align_plain_offset(handle, offset, data_start, file_size):
    """Move a nominal byte offset to the next complete VCF line boundary."""
    if offset <= data_start:
        return data_start
    if offset >= file_size:
        return file_size
    handle.seek(offset - 1)
    if handle.read(1) == b"\n":
        return offset
    handle.seek(offset)
    handle.readline()
    return handle.tell()


def read_plain_slice(vcf_file, start, end):
    """Read an aligned regular-file slice, using pread when available."""
    length = max(0, end - start)
    if length == 0:
        return b""
    with open(vcf_file, "rb", buffering=0) as handle:
        if hasattr(os, "pread"):
            parts = []
            remaining = length
            offset = start
            while remaining:
                block = os.pread(handle.fileno(), remaining, offset)
                if not block:
                    break
                parts.append(block)
                offset += len(block)
                remaining -= len(block)
            data = b"".join(parts)
        else:
            handle.seek(start)
            data = handle.read(length)
    if len(data) != length:
        raise OSError(
            "Plain VCF ended unexpectedly while reading byte range {}-{}"
            .format(start, end))
    return data


def extract_sample_names(vcf_file):
    """Compatibility wrapper returning sample names from the VCF header."""
    return read_vcf_header(vcf_file)[0]


def is_anomalous(record, num_samples):
    """Return True when a VCF row has an unexpected number of columns."""
    return len(record) != num_samples + 9


def is_snp(record):
    """Return True for single-nucleotide REF/ALT alleles, including multiallelic SNPs."""
    alt = record[4].replace(b"<NON_REF>", record[3])
    return len(record[3]) == 1 and len(alt) - alt.count(b",") == alt.count(b",") + 1


def num_genotypes(record, num_samples):
    """Count samples whose genotype field does not start with a missing allele."""
    missing = sum(1 for field in record[9:num_samples + 9] if field.startswith(b"."))
    return num_samples - missing


def get_matrix_column(record, num_samples, resolve_iupac, rng=None):
    """Transform one VCF record into one nucleotide matrix column as bytes."""
    nt_dict = {b"0": record[3].replace(b"-", b"*").upper(), b".": b"N"}
    alt = record[4].replace(b"-", b"*").replace(b"<NON_REF>", nt_dict[b"0"])
    for n, allele in enumerate(alt.split(b","), start=1):
        nt_dict[str(n).encode("ascii")] = allele

    chars = bytearray()
    choose = random.choice if rng is None else rng.choice
    genotype_cache = {}
    cleaned_cache = {}
    for field in record[9:num_samples + 9]:
        genotype = field.split(b":", 1)[0]
        if not resolve_iupac and genotype in genotype_cache:
            chars.extend(genotype_cache[genotype])
            continue
        geno_num = cleaned_cache.get(genotype)
        if geno_num is None:
            geno_num = genotype.replace(b"/", b"").replace(b"|", b"")
            cleaned_cache[genotype] = geno_num
        try:
            if resolve_iupac:
                selected = nt_dict[bytes((choose(geno_num),))]
                chars.extend(AMBIG_BYTES[selected])
            else:
                alleles = {nt_dict[bytes((code,))] for code in geno_num}
                geno_nuc = b"".join(sorted(alleles))
                consensus = AMBIG_BYTES[geno_nuc]
                genotype_cache[genotype] = consensus
                chars.extend(consensus)
        except KeyError:
            return None
    return bytes(chars)


def get_matrix_column_bin(record, num_samples):
    """Return one binary NEXUS column; unsupported genotypes become '?'."""
    chars = bytearray()
    for field in record[9:num_samples + 9]:
        genotype = field.split(b":", 1)[0]
        chars.extend(GEN_BIN_BYTES.get(genotype, b"?"))
    return bytes(chars)


def iter_vcf_chunks(vcf, target_bytes, target_records):
    """Yield contiguous binary VCF chunks using byte and record limits."""
    chunk = bytearray()
    records = 0
    for line in vcf:
        chunk.extend(line)
        if line and not line.startswith(b"#"):
            records += 1
        if len(chunk) >= target_bytes or records >= target_records:
            yield bytes(chunk)
            chunk.clear()
            records = 0
    if chunk:
        yield bytes(chunk)


def process_vcf_chunk(chunk, num_samples, min_samples_locus, need_nt,
                      need_bin, resolve_iupac, write_used, random_seed=None):
    """Parse one binary VCF chunk in a worker process."""
    rng = random.Random(random_seed) if random_seed is not None else None
    nt_data = bytearray()
    bin_data = bytearray()
    used_data = bytearray()
    malformed_lines = []

    snp_num = 0
    snp_accepted = 0
    snp_shallow = 0
    mnp_num = 0
    snp_biallelic = 0

    for raw_line in chunk.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(b"#"):
            continue

        record = line.split(b"\t")
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

        if need_nt:
            site = get_matrix_column(record, num_samples, resolve_iupac, rng)
            if site is None:
                malformed_lines.append(line)
                continue
            snp_accepted += 1
            nt_data.extend(site)
            if write_used:
                used_data.extend(record[0])
                used_data.extend(b"\t")
                used_data.extend(record[1])
                used_data.extend(b"\t")
                used_data.extend(str(num_samples_locus).encode("ascii"))
                used_data.extend(b"\n")

        if need_bin and len(record[4]) == 1:
            snp_biallelic += 1
            bin_data.extend(get_matrix_column_bin(record, num_samples))

    counters = (snp_num, snp_accepted, snp_shallow, mnp_num, snp_biallelic)
    return (bytes(nt_data), bytes(bin_data), bytes(used_data),
            malformed_lines, counters)


def ordered_parallel_results(chunks, workers, task_args, resolve_iupac):
    """Process chunks concurrently while yielding results in input order."""
    max_pending = max(1, workers + 2)
    chunk_iter = iter(enumerate(chunks))
    pending = {}
    completed = {}
    next_yield = 0

    with ProcessPoolExecutor(max_workers=workers) as executor:
        def submit_one():
            try:
                chunk_id, chunk = next(chunk_iter)
            except StopIteration:
                return False
            seed = random.getrandbits(64) if resolve_iupac else None
            future = executor.submit(process_vcf_chunk, chunk, *task_args, seed)
            pending[future] = chunk_id
            return True

        while len(pending) < max_pending and submit_one():
            pass

        while pending:
            done, _ = wait(tuple(pending), return_when=FIRST_COMPLETED)
            for future in done:
                chunk_id = pending.pop(future)
                completed[chunk_id] = future.result()
                submit_one()

            while next_yield in completed:
                yield completed.pop(next_yield)
                next_yield += 1


def write_transposed_block(handle, row_major, num_samples):
    """Write flattened row-major site data as a sample-major temporary block."""
    if not row_major:
        return None
    if len(row_major) % num_samples:
        raise ValueError("Internal error: matrix row width does not match sample count")

    width = len(row_major) // num_samples
    start = handle.tell()
    for sample_index in range(num_samples):
        handle.write(row_major[sample_index::num_samples])
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


def module_available(name):
    """Return whether an optional acceleration module can be imported."""
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, AttributeError, ValueError):
        return False


def automatic_bgzip_threads(total_workers):
    """Allocate part of the total CPU budget to HTSlib BGZF decompression."""
    if total_workers <= 2:
        return 1
    if total_workers <= 8:
        return 2
    if total_workers <= 16:
        return 4
    return min(8, max(4, total_workers // 4))


def select_input_plan(vcf_file, total_workers, requested_backend="auto",
                      decompression_threads=None, no_indexed_regions=False):
    """Inspect the file and select the fastest available compatible backend."""
    input_format = inspect_input_format(vcf_file)
    index_path = find_vcf_index(vcf_file) if input_format == "bgzf" else None
    bgzip_path = shutil.which("bgzip")
    tabix_path = shutil.which("tabix")
    isal_available = module_available("isal")
    zlib_ng_available = module_available("zlib_ng")

    if requested_backend == "plain":
        if input_format != "plain":
            raise ValueError("--input-backend plain requires an uncompressed VCF")
        backend = "plain-range" if total_workers > 1 else "plain"
    elif requested_backend == "plain-stream":
        if input_format != "plain":
            raise ValueError(
                "--input-backend plain-stream requires an uncompressed VCF")
        backend = "plain"
    elif requested_backend == "tabix":
        if input_format != "bgzf":
            raise ValueError("--input-backend tabix requires a real BGZF file")
        if index_path is None:
            raise ValueError("--input-backend tabix requires a .tbi or .csi index")
        if tabix_path is None:
            raise ValueError("--input-backend tabix requires the tabix executable")
        backend = "tabix"
    elif requested_backend == "bgzip":
        if input_format != "bgzf":
            raise ValueError("--input-backend bgzip requires a real BGZF file")
        if bgzip_path is None:
            raise ValueError("--input-backend bgzip requires the bgzip executable")
        backend = "bgzip"
    elif requested_backend == "isal":
        if input_format == "plain":
            raise ValueError("--input-backend isal requires gzip-compressed input")
        if not isal_available:
            raise ValueError("--input-backend isal requested, but python-isal is unavailable")
        backend = "isal"
    elif requested_backend == "zlib-ng":
        if input_format == "plain":
            raise ValueError("--input-backend zlib-ng requires gzip-compressed input")
        if not zlib_ng_available:
            raise ValueError("--input-backend zlib-ng requested, but python-zlib-ng is unavailable")
        backend = "zlib-ng"
    elif requested_backend == "stdlib":
        backend = "plain" if input_format == "plain" else "stdlib-gzip"
    elif requested_backend != "auto":
        raise ValueError("Unsupported input backend: {}".format(requested_backend))
    elif input_format == "plain":
        backend = "plain-range" if total_workers > 1 else "plain"
    elif input_format == "gzip":
        if isal_available:
            backend = "isal"
        elif zlib_ng_available:
            backend = "zlib-ng"
        else:
            backend = "stdlib-gzip"
    elif index_path is not None and tabix_path is not None and not no_indexed_regions:
        backend = "tabix"
    elif bgzip_path is not None and total_workers > 1:
        backend = "bgzip"
    else:
        # gzip.open correctly reads concatenated BGZF members and is the safest
        # dependency-free fallback when HTSlib is unavailable.
        backend = "stdlib-gzip"

    parser_workers = total_workers
    bgzf_threads = 0
    if backend == "bgzip":
        if decompression_threads is None:
            bgzf_threads = automatic_bgzip_threads(total_workers)
        else:
            bgzf_threads = decompression_threads
        if total_workers > 1 and bgzf_threads >= total_workers:
            raise ValueError(
                "BGZF decompression threads must be smaller than the total -t CPU budget")
        parser_workers = max(1, total_workers - bgzf_threads)
    elif decompression_threads is not None:
        raise ValueError(
            "--decompression-threads is only used by the bgzip streaming backend")

    descriptions = {
        "plain": "plain VCF sequential binary stream",
        "plain-range": "plain VCF direct parallel byte ranges",
        "stdlib-gzip": "Python stdlib gzip sequential decompression",
        "isal": "python-isal IGzip sequential decompression",
        "zlib-ng": "python-zlib-ng gzip sequential decompression",
        "bgzip": "HTSlib bgzip multithreaded BGZF stream",
        "tabix": "indexed BGZF parallel regions through tabix",
    }
    return {
        "input_format": input_format,
        "backend": backend,
        "description": descriptions[backend],
        "total_workers": total_workers,
        "parser_workers": parser_workers,
        "decompression_threads": bgzf_threads,
        "index_path": str(index_path) if index_path else None,
        "bgzip_path": bgzip_path,
        "tabix_path": tabix_path,
    }


@contextmanager
def open_vcf_stream(vcf_file, plan):
    """Yield a binary sequential stream for non-indexed input backends."""
    backend = plan["backend"]
    if backend == "plain":
        with open(vcf_file, "rb") as handle:
            yield handle
        return
    if backend == "stdlib-gzip":
        with gzip.open(vcf_file, "rb") as handle:
            yield handle
        return
    if backend == "isal":
        from isal import igzip
        with igzip.open(vcf_file, "rb") as handle:
            yield handle
        return
    if backend == "zlib-ng":
        from zlib_ng import gzip_ng
        with gzip_ng.open(vcf_file, "rb") as handle:
            yield handle
        return
    if backend != "bgzip":
        raise ValueError("Backend {} is not a sequential stream".format(backend))

    command = [
        plan["bgzip_path"], "-d", "-c", "-@",
        str(plan["decompression_threads"]), str(vcf_file),
    ]
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1024 * 1024)
    active_exception = None
    try:
        if process.stdout is None:
            raise OSError("Failed to open bgzip stdout")
        yield process.stdout
    except BaseException as exc:
        active_exception = exc
        raise
    finally:
        if process.stdout is not None:
            process.stdout.close()
        if active_exception is not None and process.poll() is None:
            process.terminate()
        stderr = process.stderr.read() if process.stderr is not None else b""
        return_code = process.wait()
        if process.stderr is not None:
            process.stderr.close()
        if active_exception is None and return_code != 0:
            raise RuntimeError(
                "bgzip decompression failed: {}".format(
                    stderr.decode("utf-8", errors="replace").strip()))


def list_tabix_contigs(tabix_path, vcf_file):
    """Return indexed contig names in the order reported by tabix."""
    result = subprocess.run(
        [tabix_path, "-l", str(vcf_file)], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            "tabix could not read the index: {}".format(
                result.stderr.decode("utf-8", errors="replace").strip()))
    contigs = [
        line.decode("utf-8") for line in result.stdout.splitlines() if line
    ]
    if not contigs:
        raise RuntimeError("tabix index contains no sequence names")
    return contigs


def build_tabix_regions(contigs, contig_lengths, workers):
    """Split indexed contigs into ordered coordinate windows when lengths exist."""
    target_regions = min(4096, max(len(contigs), workers * 4))
    known_total = sum(
        contig_lengths.get(contig, 0) for contig in contigs
        if contig_lengths.get(contig, 0) > 0 and ":" not in contig)
    if known_total > 0:
        window_size = max(1000000, int(math.ceil(known_total / target_regions)))
    else:
        window_size = None

    regions = []
    for contig in contigs:
        length = contig_lengths.get(contig)
        if not window_size or not length or length <= 0 or ":" in contig:
            regions.append({
                "query": contig, "contig": contig,
                "start": None, "end": None,
            })
            continue
        for start in range(1, length + 1, window_size):
            end = min(length, start + window_size - 1)
            regions.append({
                "query": "{}:{}-{}".format(contig, start, end),
                "contig": contig, "start": start, "end": end,
            })
    return regions


def _accumulate_counters(total, addition):
    return tuple(left + right for left, right in zip(total, addition))


def process_tabix_region(task):
    """Query, decompress, parse, and transpose one indexed BGZF region."""
    order = task["order"]
    prefix = os.path.join(task["temp_dir"], "region_{:06d}".format(order))
    nt_path = prefix + ".nt" if task["need_nt"] else None
    bin_path = prefix + ".bin" if task["need_bin"] else None
    used_path = prefix + ".used" if task["write_used"] else None
    malformed_path = prefix + ".malformed"
    nt_widths = []
    bin_widths = []
    counters = (0, 0, 0, 0, 0)

    nt_handle = open(nt_path, "wb") if nt_path else None
    bin_handle = open(bin_path, "wb") if bin_path else None
    used_handle = open(used_path, "wb") if used_path else None
    malformed_handle = open(malformed_path, "wb")
    process = None

    try:
        process = subprocess.Popen(
            [task["tabix_path"], task["vcf_file"], task["query"]],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1024 * 1024)
        if process.stdout is None:
            raise OSError("Failed to open tabix stdout")

        chunk = bytearray()
        record_count = 0
        seed_rng = random.Random(task["random_seed"])

        def flush_chunk():
            nonlocal chunk, record_count, counters
            if not chunk:
                return
            seed = seed_rng.getrandbits(64) if task["resolve_iupac"] else None
            result = process_vcf_chunk(
                bytes(chunk), task["num_samples"], task["min_samples_locus"],
                task["need_nt"], task["need_bin"], task["resolve_iupac"],
                task["write_used"], seed)
            nt_data, bin_data, used_data, malformed_lines, chunk_counters = result
            counters = _accumulate_counters(counters, chunk_counters)
            if nt_data:
                block = write_transposed_block(nt_handle, nt_data, task["num_samples"])
                nt_widths.append(block[1])
            if bin_data:
                block = write_transposed_block(bin_handle, bin_data, task["num_samples"])
                bin_widths.append(block[1])
            if used_data:
                used_handle.write(used_data)
            for malformed in malformed_lines:
                malformed_handle.write(malformed + b"\n")
            chunk = bytearray()
            record_count = 0

        for line in process.stdout:
            if task["start"] is not None:
                columns = line.split(b"\t", 2)
                if len(columns) < 2:
                    continue
                try:
                    position = int(columns[1])
                except ValueError:
                    continue
                if position < task["start"] or position > task["end"]:
                    continue
            chunk.extend(line)
            record_count += 1
            if len(chunk) >= task["chunk_bytes"] or record_count >= task["chunk_records"]:
                flush_chunk()
        flush_chunk()

        process.stdout.close()
        stderr = process.stderr.read() if process.stderr is not None else b""
        return_code = process.wait()
        if process.stderr is not None:
            process.stderr.close()
        if return_code != 0:
            raise RuntimeError(
                "tabix query '{}' failed: {}".format(
                    task["query"], stderr.decode("utf-8", errors="replace").strip()))
    except BaseException:
        if process is not None and process.poll() is None:
            process.terminate()
            process.wait()
        for path in (nt_path, bin_path, used_path, malformed_path):
            if path:
                Path(path).unlink(missing_ok=True)
        raise
    finally:
        for handle in (nt_handle, bin_handle, used_handle, malformed_handle):
            if handle is not None and not handle.closed:
                handle.close()

    return {
        "order": order,
        "query": task["query"],
        "nt_path": nt_path,
        "nt_widths": nt_widths,
        "bin_path": bin_path,
        "bin_widths": bin_widths,
        "used_path": used_path,
        "malformed_path": malformed_path,
        "counters": counters,
    }


def process_plain_range(task):
    """Read, parse, and transpose one ordered byte range of a plain VCF."""
    order = task["order"]
    prefix = os.path.join(task["temp_dir"], "range_{:06d}".format(order))
    nt_path = prefix + ".nt" if task["need_nt"] else None
    bin_path = prefix + ".bin" if task["need_bin"] else None
    used_path = prefix + ".used" if task["write_used"] else None
    malformed_path = prefix + ".malformed"
    nt_widths = []
    bin_widths = []
    counters = (0, 0, 0, 0, 0)

    nt_handle = open(nt_path, "wb") if nt_path else None
    bin_handle = open(bin_path, "wb") if bin_path else None
    used_handle = open(used_path, "wb") if used_path else None
    malformed_handle = open(malformed_path, "wb")

    try:
        file_size = Path(task["vcf_file"]).stat().st_size
        with open(task["vcf_file"], "rb", buffering=8 * MIB) as handle:
            aligned_start = align_plain_offset(
                handle, task["start"], task["data_start"], file_size)
            aligned_end = align_plain_offset(
                handle, task["end"], task["data_start"], file_size)
            if hasattr(os, "posix_fadvise") and hasattr(os, "POSIX_FADV_SEQUENTIAL"):
                try:
                    os.posix_fadvise(
                        handle.fileno(), aligned_start,
                        max(0, aligned_end - aligned_start),
                        os.POSIX_FADV_SEQUENTIAL)
                except OSError:
                    pass

        chunk = read_plain_slice(task["vcf_file"], aligned_start, aligned_end)
        seed = task["random_seed"] if task["resolve_iupac"] else None
        result = process_vcf_chunk(
            chunk, task["num_samples"], task["min_samples_locus"],
            task["need_nt"], task["need_bin"], task["resolve_iupac"],
            task["write_used"], seed)
        nt_data, bin_data, used_data, malformed_lines, counters = result
        if nt_data:
            block = write_transposed_block(
                nt_handle, nt_data, task["num_samples"])
            nt_widths.append(block[1])
        if bin_data:
            block = write_transposed_block(
                bin_handle, bin_data, task["num_samples"])
            bin_widths.append(block[1])
        if used_data:
            used_handle.write(used_data)
        for malformed in malformed_lines:
            malformed_handle.write(malformed + b"\n")
    except BaseException:
        for path in (nt_path, bin_path, used_path, malformed_path):
            if path:
                Path(path).unlink(missing_ok=True)
        raise
    finally:
        for handle in (nt_handle, bin_handle, used_handle, malformed_handle):
            if handle is not None and not handle.closed:
                handle.close()

    return {
        "order": order,
        "start": task["start"],
        "end": task["end"],
        "nt_path": nt_path,
        "nt_widths": nt_widths,
        "bin_path": bin_path,
        "bin_widths": bin_widths,
        "used_path": used_path,
        "malformed_path": malformed_path,
        "counters": counters,
    }


def merge_region_blocks(global_handle, local_path, widths, global_blocks,
                        num_samples):
    """Append worker-generated transposed blocks and register global offsets."""
    if local_path is None:
        return
    local = Path(local_path)
    if not local.exists():
        return
    offset = global_handle.tell()
    with local.open("rb") as source:
        shutil.copyfileobj(source, global_handle, length=8 * MIB)
    for width in widths:
        global_blocks.append((offset, width))
        offset += num_samples * width
    local.unlink(missing_ok=True)


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


def profile_vcf(vcf_path, input_format=None):
    """Estimate uncompressed size, record count, and average record width."""
    path = Path(vcf_path)
    stored_bytes = path.stat().st_size
    input_format = input_format or inspect_input_format(path)
    is_compressed = input_format in ("gzip", "bgzf")

    if is_compressed:
        sample, compressed_consumed, reached_eof = sample_gzip_prefix(path)
        sampled_uncompressed = len(sample)
        sampled_record_bytes, sampled_records = summarize_sample(
            sample, reached_eof)
        if reached_eof:
            estimated_uncompressed = sampled_uncompressed
            size_source = "exact compressed scan"
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
        "input_format": input_format,
        "size_source": size_source,
    }


def choose_chunk_settings(workers, vcf_path, chunk_size_mb=None, profile=None):
    """Choose chunk limits from parser workers and estimated VCF workload."""
    profile = profile or profile_vcf(vcf_path)
    estimated_bytes = max(1, profile["estimated_uncompressed_bytes"])
    estimated_records = profile["estimated_records"]

    if estimated_records <= 1000000:
        tasks_per_worker = 8
    elif estimated_records <= 5000000:
        tasks_per_worker = 10
    else:
        tasks_per_worker = 12

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
        source = "auto from parser CPU count and estimated VCF workload"

    if chunk_size_mb is None:
        target_bytes = int(clamp(target_bytes,
                                 MIN_CHUNK_TARGET_BYTES,
                                 MAX_CHUNK_TARGET_BYTES))
        target_bytes = max(MIN_CHUNK_TARGET_BYTES,
                           int(round(target_bytes / MIB)) * MIB)
    else:
        target_bytes = int(target_bytes)

    target_records = int(clamp(target_records,
                               MIN_CHUNK_RECORDS,
                               MAX_CHUNK_RECORDS))
    record_quantum = 500 if target_records < 10000 else 1000
    target_records = max(
        MIN_CHUNK_RECORDS,
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


def print_malformed_lines(lines):
    for malformed in lines:
        print("Skipping malformed line:\n{}".format(
            malformed.decode("utf-8", errors="replace")))


def print_malformed_file(path):
    if not path:
        return
    malformed_path = Path(path)
    if malformed_path.exists():
        with malformed_path.open("rb") as handle:
            for malformed in handle:
                print("Skipping malformed line:\n{}".format(
                    malformed.rstrip(b"\r\n").decode("utf-8", errors="replace")))
        malformed_path.unlink(missing_ok=True)


def run_stream_conversion(vcf_file, plan, chunk_bytes, chunk_records,
                          task_args, nt_temp, nt_blocks, bin_temp, bin_blocks,
                          used_sites):
    """Run plain, ordinary-gzip, or streaming-BGZF conversion."""
    counters = (0, 0, 0, 0, 0)
    next_progress = 500000
    with open_vcf_stream(vcf_file, plan) as vcf:
        chunks = iter_vcf_chunks(vcf, chunk_bytes, chunk_records)
        if plan["parser_workers"] == 1:
            results = (
                process_vcf_chunk(chunk, *task_args, None)
                for chunk in chunks
            )
        else:
            results = ordered_parallel_results(
                chunks, plan["parser_workers"], task_args, task_args[4])

        for nt_data, bin_data, used_data, malformed_lines, addition in results:
            print_malformed_lines(malformed_lines)
            counters = _accumulate_counters(counters, addition)
            while counters[0] >= next_progress:
                print("{:d} genotypes processed.".format(next_progress))
                next_progress += 500000
            if nt_data:
                nt_blocks.append(write_transposed_block(
                    nt_temp, nt_data, task_args[0]))
            if bin_data:
                bin_blocks.append(write_transposed_block(
                    bin_temp, bin_data, task_args[0]))
            if used_data and used_sites is not None:
                used_sites.write(used_data)
    return counters


def run_plain_range_conversion(vcf_file, plan, chunk_bytes, chunk_records,
                               target_tasks, task_args, nt_temp, nt_blocks,
                               bin_temp, bin_blocks, used_sites,
                               scratch_parent):
    """Run direct file-range parallelism for an uncompressed VCF."""
    data_start = find_plain_data_start(vcf_file)
    ranges = build_plain_ranges(
        vcf_file, data_start, target_tasks, chunk_bytes)
    print("Direct byte ranges: {:,d}".format(len(ranges)))
    if not ranges:
        return (0, 0, 0, 0, 0)

    range_temp_dir = tempfile.mkdtemp(
        prefix=".vcf2phylip_ranges_", dir=str(scratch_parent))
    counters = (0, 0, 0, 0, 0)
    next_progress = 500000

    tasks = []
    for range_spec in ranges:
        task = dict(range_spec)
        task.update({
            "vcf_file": str(vcf_file),
            "temp_dir": range_temp_dir,
            "chunk_bytes": chunk_bytes,
            "chunk_records": chunk_records,
            "num_samples": task_args[0],
            "min_samples_locus": task_args[1],
            "need_nt": task_args[2],
            "need_bin": task_args[3],
            "resolve_iupac": task_args[4],
            "write_used": task_args[5],
            "random_seed": random.getrandbits(64),
        })
        tasks.append(task)

    try:
        if plan["parser_workers"] == 1:
            results = map(process_plain_range, tasks)
            executor = None
        else:
            executor = ProcessPoolExecutor(max_workers=plan["parser_workers"])
            results = executor.map(process_plain_range, tasks, chunksize=1)

        try:
            for result in results:
                print_malformed_file(result["malformed_path"])
                merge_region_blocks(
                    nt_temp, result["nt_path"], result["nt_widths"],
                    nt_blocks, task_args[0])
                merge_region_blocks(
                    bin_temp, result["bin_path"], result["bin_widths"],
                    bin_blocks, task_args[0])
                if result["used_path"]:
                    used_path = Path(result["used_path"])
                    if used_path.exists() and used_sites is not None:
                        with used_path.open("rb") as source:
                            shutil.copyfileobj(source, used_sites, length=8 * MIB)
                    used_path.unlink(missing_ok=True)

                counters = _accumulate_counters(counters, result["counters"])
                while counters[0] >= next_progress:
                    print("{:d} genotypes processed.".format(next_progress))
                    next_progress += 500000
        finally:
            if executor is not None:
                executor.shutdown(wait=True)
    finally:
        shutil.rmtree(range_temp_dir, ignore_errors=True)
    return counters


def run_indexed_conversion(vcf_file, plan, contig_lengths, chunk_bytes,
                           chunk_records, task_args, nt_temp, nt_blocks,
                           bin_temp, bin_blocks, used_sites, scratch_parent):
    """Run true region-parallel BGZF decompression and genotype parsing."""
    contigs = list_tabix_contigs(plan["tabix_path"], vcf_file)
    regions = build_tabix_regions(contigs, contig_lengths,
                                  plan["parser_workers"])
    print("Indexed regions: {:,d} across {:,d} contigs".format(
        len(regions), len(contigs)))

    region_temp_dir = tempfile.mkdtemp(
        prefix=".vcf2phylip_regions_", dir=str(scratch_parent))
    counters = (0, 0, 0, 0, 0)
    next_progress = 500000

    tasks = []
    for order, region in enumerate(regions):
        tasks.append({
            "order": order,
            "tabix_path": plan["tabix_path"],
            "vcf_file": str(vcf_file),
            "query": region["query"],
            "start": region["start"],
            "end": region["end"],
            "temp_dir": region_temp_dir,
            "chunk_bytes": chunk_bytes,
            "chunk_records": chunk_records,
            "num_samples": task_args[0],
            "min_samples_locus": task_args[1],
            "need_nt": task_args[2],
            "need_bin": task_args[3],
            "resolve_iupac": task_args[4],
            "write_used": task_args[5],
            "random_seed": random.getrandbits(64),
        })

    try:
        if plan["parser_workers"] == 1:
            results = map(process_tabix_region, tasks)
            executor = None
        else:
            executor = ProcessPoolExecutor(max_workers=plan["parser_workers"])
            results = executor.map(process_tabix_region, tasks, chunksize=1)

        try:
            for result in results:
                print_malformed_file(result["malformed_path"])
                merge_region_blocks(
                    nt_temp, result["nt_path"], result["nt_widths"],
                    nt_blocks, task_args[0])
                merge_region_blocks(
                    bin_temp, result["bin_path"], result["bin_widths"],
                    bin_blocks, task_args[0])
                if result["used_path"]:
                    used_path = Path(result["used_path"])
                    if used_path.exists() and used_sites is not None:
                        with used_path.open("rb") as source:
                            shutil.copyfileobj(source, used_sites, length=8 * MIB)
                    used_path.unlink(missing_ok=True)

                counters = _accumulate_counters(counters, result["counters"])
                while counters[0] >= next_progress:
                    print("{:d} genotypes processed.".format(next_progress))
                    next_progress += 500000
        finally:
            if executor is not None:
                executor.shutdown(wait=True)
    finally:
        shutil.rmtree(region_temp_dir, ignore_errors=True)
    return counters


def write_final_matrices(args, outfile, sample_names, outgroup,
                         snp_accepted, snp_biallelic,
                         nt_temp_path, nt_blocks, bin_temp_path, bin_blocks):
    """Write all requested matrix formats using the optimized temp blocks."""
    num_samples = len(sample_names)
    output_phy = open(outfile + ".phy", "w", encoding="utf-8") \
        if not args.phylipdisable else None
    output_fas = open(outfile + ".fasta", "w", encoding="utf-8") \
        if args.fasta else None
    output_nex = open(outfile + ".nexus", "w", encoding="utf-8") \
        if args.nexus else None
    output_nexbin = open(outfile + ".bin.nexus", "w", encoding="utf-8") \
        if args.nexusbin else None

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
                        lambda sample_index: assemble_sample(
                            sample_index, nt_temp_path, nt_blocks,
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
                        if nt_temp_path is not None:
                            print("Outgroup, '{}', added to the matrix(ces).".format(outgroup))
                        if bin_temp_path is not None:
                            print("Outgroup, '{}', added to the binary matrix.".format(outgroup))
                    else:
                        if nt_temp_path is not None:
                            print("Sample {:d} of {:d}, '{}', added to the nucleotide matrix(ces).".format(
                                sample_index + 1, num_samples, name))
                        if bin_temp_path is not None:
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
        for output in (output_phy, output_fas, output_nex, output_nexbin):
            if output is not None:
                output.close()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", dest="filename", required=True,
                        help="Name of the input VCF file; plain gzip and BGZF are auto-detected")
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
                        help="Total CPU budget; default: 100%% of CPUs available to the job")
    parser.add_argument("--chunk-size-mb", type=positive_int, default=None,
                        help="Override auto byte limit in MiB; record limit remains adaptive")
    parser.add_argument(
        "--input-backend", default="auto",
        choices=("auto", "plain", "plain-stream", "stdlib", "isal",
                 "zlib-ng", "bgzip", "tabix"),
        help="Override automatic input backend selection (default: auto)")
    parser.add_argument(
        "--decompression-threads", type=positive_int, default=None,
        help="Threads reserved for bgzip streaming decompression")
    parser.add_argument(
        "--no-indexed-regions", action="store_true",
        help="Do not use TBI/CSI region parallelism even when available")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {version}".format(version=__version__))
    args = parser.parse_args()

    if args.threads is None:
        args.threads, worker_source = detect_available_workers()
        worker_description = "100% auto-detected from {}; not reduced by VCF size".format(
            worker_source)
    else:
        worker_description = "user specified total CPU budget"

    outgroup = args.outgroup.split(",")[0].split(";")[0]
    need_nt = args.fasta or args.nexus or not args.phylipdisable
    need_bin = args.nexusbin
    input_path = Path(args.filename)
    if not input_path.exists():
        print("\nInput VCF file not found, please verify the provided path")
        sys.exit(1)

    try:
        plan = select_input_plan(
            args.filename, args.threads, args.input_backend,
            args.decompression_threads, args.no_indexed_regions)
    except ValueError as error:
        parser.error(str(error))

    profile = profile_vcf(args.filename, plan["input_format"])
    chunk_bytes, chunk_records, chunk_source, target_tasks, profile = \
        choose_chunk_settings(
            plan["parser_workers"], args.filename, args.chunk_size_mb, profile)
    sample_names, contig_lengths, _ = read_vcf_header(
        args.filename, plan["input_format"])
    num_samples = len(sample_names)
    if num_samples == 0:
        print("\nSample names not found in VCF, your file may be corrupt or missing the header.\n")
        sys.exit(1)

    print("\nConverting file '{}':\n".format(args.filename))
    print("Number of samples in VCF: {:d}".format(num_samples))
    print("Parallel workers: {:d} ({})".format(args.threads, worker_description))
    print("Detected input format: {}".format(plan["input_format"].upper()))
    print("Input backend: {}".format(plan["description"]))
    if plan["backend"] == "bgzip":
        print("CPU allocation: {:d} bgzip threads + {:d} parser workers = {:d} total".format(
            plan["decompression_threads"], plan["parser_workers"], args.threads))
    elif plan["backend"] == "tabix":
        print("Indexed sidecar: {}".format(plan["index_path"]))
        print("Region workers: {:d}".format(plan["parser_workers"]))
    elif plan["backend"] == "plain-range":
        print("Direct range workers: {:d}".format(plan["parser_workers"]))
    else:
        print("Parser workers: {:d}".format(plan["parser_workers"]))

    if profile["is_compressed"]:
        print("VCF size: {} stored; about {} uncompressed ({})".format(
            format_binary_size(profile["stored_bytes"]),
            format_binary_size(profile["estimated_uncompressed_bytes"]),
            profile["size_source"]))
    else:
        print("VCF size: {} ({})".format(
            format_binary_size(profile["estimated_uncompressed_bytes"]),
            profile["size_source"]))
    if profile["estimated_records"] > 0:
        print("Estimated data records: {:,}; average sampled row: {:.0f} bytes".format(
            profile["estimated_records"], profile["average_record_bytes"]))
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
            used_sites = open(outfile + ".used_sites.tsv", "wb")
            used_sites.write(b"#CHROM\tPOS\tNUM_SAMPLES\n")
        if need_nt:
            nt_temp = open(nt_temp_path, "wb")
        if need_bin:
            bin_temp = open(bin_temp_path, "wb")

        task_args = (
            num_samples, args.min_samples_locus, need_nt, need_bin,
            args.resolve_iupac, args.write_used,
        )
        if plan["backend"] == "plain-range":
            try:
                counters = run_plain_range_conversion(
                    args.filename, plan, chunk_bytes, chunk_records,
                    target_tasks, task_args, nt_temp, nt_blocks, bin_temp,
                    bin_blocks, used_sites, output_folder)
            except (OSError, RuntimeError, ValueError) as error:
                if args.input_backend not in ("auto", "plain"):
                    raise
                print("WARNING: direct plain-VCF range mode failed: {}".format(error))
                print("Falling back to the sequential plain-VCF reader without changing output semantics.")
                if nt_temp is not None:
                    nt_temp.seek(0)
                    nt_temp.truncate()
                if bin_temp is not None:
                    bin_temp.seek(0)
                    bin_temp.truncate()
                nt_blocks.clear()
                bin_blocks.clear()
                if used_sites is not None:
                    used_sites.seek(0)
                    used_sites.truncate()
                    used_sites.write(b"#CHROM\tPOS\tNUM_SAMPLES\n")

                fallback_plan = select_input_plan(
                    args.filename, args.threads, "plain-stream")
                print("Fallback backend: {}".format(
                    fallback_plan["description"]))
                counters = run_stream_conversion(
                    args.filename, fallback_plan, chunk_bytes, chunk_records,
                    task_args, nt_temp, nt_blocks, bin_temp, bin_blocks,
                    used_sites)
        elif plan["backend"] == "tabix":
            try:
                counters = run_indexed_conversion(
                    args.filename, plan, contig_lengths, chunk_bytes, chunk_records,
                    task_args, nt_temp, nt_blocks, bin_temp, bin_blocks,
                    used_sites, output_folder)
            except (OSError, RuntimeError, subprocess.SubprocessError) as error:
                if args.input_backend != "auto":
                    raise
                print("WARNING: indexed BGZF mode failed: {}".format(error))
                print("Falling back to sequential BGZF streaming without changing output semantics.")
                if nt_temp is not None:
                    nt_temp.seek(0)
                    nt_temp.truncate()
                if bin_temp is not None:
                    bin_temp.seek(0)
                    bin_temp.truncate()
                nt_blocks.clear()
                bin_blocks.clear()
                if used_sites is not None:
                    used_sites.seek(0)
                    used_sites.truncate()
                    used_sites.write(b"#CHROM\tPOS\tNUM_SAMPLES\n")

                fallback_plan = select_input_plan(
                    args.filename, args.threads, "auto", None, True)
                fallback_chunk_bytes, fallback_chunk_records, _, _, _ = \
                    choose_chunk_settings(
                        fallback_plan["parser_workers"], args.filename,
                        args.chunk_size_mb, profile)
                print("Fallback backend: {}".format(fallback_plan["description"]))
                if fallback_plan["backend"] == "bgzip":
                    print("Fallback CPU allocation: {:d} bgzip threads + {:d} parser workers".format(
                        fallback_plan["decompression_threads"],
                        fallback_plan["parser_workers"]))
                counters = run_stream_conversion(
                    args.filename, fallback_plan, fallback_chunk_bytes,
                    fallback_chunk_records, task_args, nt_temp, nt_blocks,
                    bin_temp, bin_blocks, used_sites)
        else:
            counters = run_stream_conversion(
                args.filename, plan, chunk_bytes, chunk_records, task_args,
                nt_temp, nt_blocks, bin_temp, bin_blocks, used_sites)

        if nt_temp is not None:
            nt_temp.close()
            nt_temp = None
        if bin_temp is not None:
            bin_temp.close()
            bin_temp = None
        if used_sites is not None:
            used_sites.close()
            used_sites = None

        snp_num, snp_accepted, snp_shallow, mnp_num, snp_biallelic = counters
        print("Total of genotypes processed: {:d}".format(snp_num))
        print("Genotypes excluded because they exceeded the amount of missing data allowed: {:d}".format(snp_shallow))
        print("Genotypes that passed missing data filter but were excluded for being MNPs: {:d}".format(mnp_num))
        print("SNPs that passed the filters: {:d}".format(snp_accepted))
        if args.nexusbin:
            print("Biallelic SNPs selected for binary NEXUS: {:d}".format(snp_biallelic))
        if args.write_used:
            print("Used sites saved to: '" + outfile + ".used_sites.tsv'")
        print("")

        write_final_matrices(
            args, outfile, sample_names, outgroup, snp_accepted,
            snp_biallelic, nt_temp_path, nt_blocks,
            bin_temp_path, bin_blocks)

        if not args.phylipdisable:
            print("PHYLIP matrix saved to: " + outfile + ".phy")
        if args.fasta:
            print("FASTA matrix saved to: " + outfile + ".fasta")
        if args.nexus:
            print("NEXUS matrix saved to: " + outfile + ".nexus")
        if args.nexusbin:
            print("BINARY NEXUS matrix saved to: " + outfile + ".bin.nexus")
    finally:
        for handle in (nt_temp, bin_temp, used_sites):
            if handle is not None:
                handle.close()
        if nt_temp_path is not None:
            Path(nt_temp_path).unlink(missing_ok=True)
        if bin_temp_path is not None:
            Path(bin_temp_path).unlink(missing_ok=True)

    print("\nDone!\n")


if __name__ == "__main__":
    main()
