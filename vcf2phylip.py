#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Convert SNPs in VCF format into PHYLIP, FASTA, NEXUS, or binary NEXUS
matrices for phylogenetic analysis.

Multithreaded edition preserving the command-line behavior of vcf2phylip v2.9
with ``-t/--threads`` (default=0, auto-detect CPU cores). VCF chunks are
parsed with multiple worker processes; final sample sequences are assembled
from transposed blocks stored during Phase 1, avoiding a full-file scan per
sample.

Any ploidy is allowed, but binary NEXUS is produced only for diploid VCFs.
"""

__author__      = "Edgardo M. Ortiz"
__credits__     = "Juan D. Palacio-Mejía"
__modifier__    = "Ma Wenxin"
__version__     = "2.9-mt2"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2026-07-22"

import argparse
import gzip
import os
import random
import sys
from collections import deque
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

# Target amount of VCF text submitted to each parsing task.
CHUNK_TARGET_BYTES = 4 * 1024 * 1024

# Dictionary of IUPAC ambiguities for nucleotides
AMBIG = {
    "A"    :"A", "C"    :"C", "G"    :"G", "N"    :"N", "T"    :"T",
    "*A"   :"a", "*C"   :"c", "*G"   :"g", "*N"   :"n", "*T"   :"t",
    "AC"   :"M", "AG"   :"R", "AN"   :"a", "AT"   :"W", "CG"   :"S",
    "CN"   :"c", "CT"   :"Y", "GN"   :"g", "GT"   :"K", "NT"   :"t",
    "*AC"  :"m", "*AG"  :"r", "*AN"  :"a", "*AT"  :"w", "*CG"  :"s",
    "*CN"  :"c", "*CT"  :"y", "*GN"  :"g", "*GT"  :"k", "*NT"  :"t",
    "ACG"  :"V", "ACN"  :"m", "ACT"  :"H", "AGN"  :"r", "AGT"  :"D",
    "ANT"  :"w", "CGN"  :"s", "CGT"  :"B", "CNT"  :"y", "GNT"  :"k",
    "*ACG" :"v", "*ACN" :"m", "*ACT" :"h", "*AGN" :"r", "*AGT" :"d",
    "*ANT" :"w", "*CGN" :"s", "*CGT" :"b", "*CNT" :"y", "*GNT" :"k",
    "ACGN" :"v", "ACGT" :"N", "ACNT" :"h", "AGNT" :"d", "CGNT" :"b",
    "*ACGN":"v", "*ACGT":"N", "*ACNT":"h", "*AGNT":"d", "*CGNT":"b",
    "*"    :"-", "*ACGNT":"N",
}

# Dictionary for translating biallelic SNPs into SNAPP, only for diploid VCF
GEN_BIN = {
    "./.":"?",
    ".|.":"?",
    "0/0":"0",
    "0|0":"0",
    "0/1":"1",
    "0|1":"1",
    "1/0":"1",
    "1|0":"1",
    "1/1":"2",
    "1|1":"2",
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
    """Return True for single-nucleotide REF/ALT alleles."""
    alt = record[4].replace("<NON_REF>", record[3])
    return len(record[3]) == 1 and len(alt) - alt.count(",") == alt.count(",") + 1


def num_genotypes(record, num_samples):
    """Count samples whose genotype field does not start with a missing allele."""
    missing = sum(1 for field in record[9:num_samples + 9] if field.startswith("."))
    return num_samples - missing


def get_matrix_column(record, num_samples, resolve_IUPAC, rng=None):
    """Transform one VCF record into one nucleotide matrix column."""
    nt_dict = {str(0): record[3].replace("-", "*").upper(), ".": "N"}
    alt = record[4].replace("-", "*").replace("<NON_REF>", nt_dict["0"])
    alt = alt.split(",")
    for n in range(len(alt)):
        nt_dict[str(n + 1)] = alt[n]
    column = ""
    choose = random.choice if rng is None else rng.choice
    for i in range(9, num_samples + 9):
        geno_num = record[i].split(":")[0].replace("/", "").replace("|", "")
        try:
            geno_nuc = "".join(sorted(set([nt_dict[j] for j in geno_num])))
        except KeyError:
            return "malformed"
        if resolve_IUPAC is False:
            column += AMBIG[geno_nuc]
        else:
            column += AMBIG[nt_dict[choose(geno_num)]]
    return column


def get_matrix_column_bin(record, num_samples):
    """Return one binary NEXUS column; unsupported genotypes become '?'."""
    column = ""
    for i in range(9, num_samples + 9):
        genotype = record[i].split(":")[0]
        if genotype in GEN_BIN:
            column += GEN_BIN[genotype]
        else:
            column += "?"
    return column


# ---------------------------------------------------------------------------
# Phase 1 helpers — VCF parsing
# ---------------------------------------------------------------------------

def iter_vcf_chunks(vcf, target_bytes=CHUNK_TARGET_BYTES):
    """Yield bounded lists of VCF lines without loading the full file."""
    chunk = []
    size = 0
    for line in vcf:
        chunk.append(line)
        size += len(line)
        if size >= target_bytes:
            yield chunk
            chunk = []
            size = 0
    if chunk:
        yield chunk


def process_vcf_chunk(lines, num_samples, min_samples_locus, need_nt,
                      need_bin, resolve_iupac, write_used, random_seed=None):
    """Parse one VCF chunk. Safe to run in a worker process."""
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

        # If nucleotide matrices are requested and genotype cannot be decoded,
        # skip this site for all outputs (preserves v2.9 behavior).
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
    max_pending = max(1, workers * 2)
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


# ---------------------------------------------------------------------------
# Phase 2 helpers — transposed-block assembly
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def positive_int(value):
    """argparse validator for non-negative integers."""
    try:
        integer = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if integer < 0:
        raise argparse.ArgumentTypeError("must be >= 0 (0 = auto-detect CPU cores)")
    return integer


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input",
        action = "store",
        dest = "filename",
        required = True,
        help = "Name of the input VCF file, can be gzipped")
    parser.add_argument("--output-folder",
        action = "store",
        dest = "folder",
        default = "./",
        help = "Output folder name, it will be created if it does not exist (same folder as input by "
               "default)")
    parser.add_argument("--output-prefix",
        action = "store",
        dest = "prefix",
        help = "Prefix for output filenames (same as the input VCF filename without the extension by "
               "default)")
    parser.add_argument("-m", "--min-samples-locus",
        action = "store",
        dest = "min_samples_locus",
        type = int,
        default = 4,
        help = "Minimum of samples required to be present at a locus (default=4)")
    parser.add_argument("-o", "--outgroup",
        action = "store",
        dest = "outgroup",
        default = "",
        help = "Name of the outgroup in the matrix. Sequence will be written as first taxon in the "
               "alignment.")
    parser.add_argument("-p", "--phylip-disable",
        action = "store_true",
        dest = "phylipdisable",
        help = "A PHYLIP matrix is written by default unless you enable this flag")
    parser.add_argument("-f", "--fasta",
        action = "store_true",
        dest = "fasta",
        help = "Write a FASTA matrix (disabled by default)")
    parser.add_argument("-n", "--nexus",
        action = "store_true",
        dest = "nexus",
        help = "Write a NEXUS matrix (disabled by default)")
    parser.add_argument("-b", "--nexus-binary",
        action = "store_true",
        dest = "nexusbin",
        help = "Write a binary NEXUS matrix for analysis of biallelic SNPs in SNAPP, only diploid "
               "genotypes will be processed (disabled by default)")
    parser.add_argument("-r", "--resolve-IUPAC",
        action = "store_true",
        dest = "resolve_IUPAC",
        help = "Randomly resolve heterozygous genotypes to avoid IUPAC ambiguities in the matrices "
               "(disabled by default)")
    parser.add_argument("-w", "--write-used-sites",
        action = "store_true",
        dest = "write_used",
        help = "Save the list of coordinates that passed the filters and were used in the alignments "
               "(disabled by default)")
    parser.add_argument("-t", "--threads",
        type = positive_int,
        default = 0,
        help = "Number of parallel processes (default=0, auto-detect CPU cores)")
    parser.add_argument("-v", "--version",
        action = "version",
        version = "%(prog)s {version}".format(version=__version__))
    args = parser.parse_args()

    # Auto-detect CPU cores if threads=0
    if args.threads <= 0:
        args.threads = os.cpu_count() or 1

    outgroup = args.outgroup.split(",")[0].split(";")[0]
    need_nt = args.fasta or args.nexus or not args.phylipdisable
    need_bin = args.nexusbin

    # Get samples names and number of samples in VCF
    input_path = Path(args.filename)
    if not input_path.exists():
        print("\nInput VCF file not found, please verify the provided path")
        sys.exit(1)

    sample_names = extract_sample_names(args.filename)
    num_samples = len(sample_names)
    if num_samples == 0:
        print("\nSample names not found in VCF, your file may be corrupt or missing the header.\n")
        sys.exit(1)
    print("\nConverting file '{}':\n".format(args.filename))
    print("Number of samples in VCF: {:d}".format(num_samples))
    print("Parallel workers: {:d}".format(args.threads))

    # If the 'min_samples_locus' is larger than the actual number of samples in VCF readjust it
    args.min_samples_locus = min(num_samples, args.min_samples_locus)

    # Output filename will be the same as input file, indicating the minimum of samples specified
    if not args.prefix:
        parts = input_path.name.split(".")
        args.prefix = []
        for p in parts:
            if p.lower() == "vcf":
                break
            else:
                args.prefix.append(p)
        args.prefix = ".".join(args.prefix)
    args.prefix += ".min" + str(args.min_samples_locus)

    # Check if outfolder exists, create it if it doesn't
    output_folder = Path(args.folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    outfile = str(output_folder / args.prefix)

    # Temp file paths and block index lists for Phase 2 assembly
    nt_temp_path = outfile + ".tmp" if need_nt else None
    bin_temp_path = outfile + ".bin.tmp" if need_bin else None
    nt_blocks = []
    bin_blocks = []

    # Resource handles — all cleaned up in finally
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

        ##########################
        # PROCESS GENOTYPES IN VCF

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
            args.resolve_IUPAC,
            args.write_used,
        )

        with opener(args.filename, "rt") as vcf:
            chunks = iter_vcf_chunks(vcf)
            if args.threads == 1:
                results = (
                    process_vcf_chunk(chunk, *task_args, None)
                    for chunk in chunks
                )
            else:
                results = ordered_parallel_results(
                    chunks, args.threads, task_args, args.resolve_IUPAC)

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

                # Write transposed blocks to temp files
                if nt_rows:
                    block = write_transposed_block(nt_temp, nt_rows, num_samples)
                    nt_blocks.append(block)
                if bin_rows:
                    block = write_transposed_block(bin_temp, bin_rows, num_samples)
                    bin_blocks.append(block)
                if used_sites is not None:
                    for chrom, pos, present in used_rows:
                        used_sites.write("{}\t{}\t{}\n".format(chrom, pos, present))

        # Close Phase 1 temp/output files before Phase 2
        if nt_temp is not None:
            nt_temp.close()
            nt_temp = None
        if bin_temp is not None:
            bin_temp.close()
            bin_temp = None
        if used_sites is not None:
            used_sites.close()
            used_sites = None

        # Print useful information about filtering of SNPs
        print("Total of genotypes processed: {:d}".format(snp_num))
        print("Genotypes excluded because they exceeded the amount "
              "of missing data allowed: {:d}".format(snp_shallow))
        print("Genotypes that passed missing data filter but were "
              "excluded for being MNPs: {:d}".format(mnp_num))
        print("SNPs that passed the filters: {:d}".format(snp_accepted))
        if need_bin:
            print("Biallelic SNPs selected for binary NEXUS: {:d}".format(snp_biallelic))
        if args.write_used:
            print("Used sites saved to: '" + outfile + ".used_sites.tsv'")
        print("")

        #######################
        # WRITE OUTPUT MATRICES

        output_phy = open(outfile + ".phy", "w", encoding="utf-8") if not args.phylipdisable else None
        output_fas = open(outfile + ".fasta", "w", encoding="utf-8") if args.fasta else None
        output_nex = open(outfile + ".nexus", "w", encoding="utf-8") if args.nexus else None
        output_nexbin = open(outfile + ".bin.nexus", "w", encoding="utf-8") if args.nexusbin else None

        try:
            if output_phy is not None:
                output_phy.write("{:d} {:d}\n".format(num_samples, snp_accepted))
            if output_nex is not None:
                output_nex.write("#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={:d} NCHAR={:d};\n\tFORMAT "
                                 "DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n".format(num_samples, snp_accepted))
            if output_nexbin is not None:
                output_nexbin.write("#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={:d} NCHAR={:d};\n\tFORMAT "
                                    "DATATYPE=SNP MISSING=? GAP=- ;\nMATRIX\n".format(num_samples, snp_biallelic))

            # Get length of longest sequence name
            len_longest_name = max(len(name) for name in sample_names)

            # Write outgroup as first sequence in alignment if the name is specified
            idx_outgroup = sample_names.index(outgroup) if outgroup in sample_names else None

            # Build ordered sample list: outgroup first, then ingroup
            sample_order = []
            if idx_outgroup is not None:
                sample_order.append(idx_outgroup)
            sample_order.extend(i for i in range(num_samples) if i != idx_outgroup)

            # Assemble sequences from transposed blocks using threads
            worker_count = min(args.threads, max(1, len(sample_order)))
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
        # Clean up any handles still open (e.g. if an exception occurred)
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
