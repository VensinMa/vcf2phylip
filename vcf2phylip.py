#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The script converts a collection of SNPs in VCF format into a PHYLIP, FASTA,
NEXUS, or binary NEXUS file for phylogenetic analysis. The code is optimized
to process VCF files with sizes >1GB. For small VCF files the algorithm slows
down as the number of taxa increases (but is still fast).

Any ploidy is allowed, but binary NEXUS is produced only for diploid VCFs.

Multithreaded version based on vcf2phylip v2.9 by Edgardo M. Ortiz.
Parallelized:
  - Phase 1: VCF parsing and genotype conversion (multiprocessing batches)
  - Phase 2: Matrix transposition and output writing (multiprocessing per sample)
"""

__author__      = "Edgardo M. Ortiz"
__credits__     = "Juan D. Palacio-Mejía"
__version__     = "2.9-parallel"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2023-07-07"

import argparse
import gzip
import os
import random
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Dictionary of IUPAC ambiguities for nucleotides
# '*' is a deletion in GATK, deletions are ignored in consensus, lowercase consensus is used when an
# 'N' or '*' is part of the genotype. Capitalization is used by some software but ignored by Geneious
# for example
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
# 0 is homozygous reference
# 1 is heterozygous
# 2 is homozygous alternative
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
    """
    Extract sample names from VCF file
    """
    if vcf_file.lower().endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
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
    """
    Determine if the number of samples in current record corresponds to number of samples described
    in the line '#CHROM'
    """
    return bool(len(record) != num_samples + 9)


def is_snp(record):
    """
    Determine if current VCF record is a SNP (single nucleotide polymorphism) as opposed to MNP
    (multinucleotide polymorphism)
    """
    # <NON_REF> must be replaced by the REF in the ALT field for GVCFs from GATK
    alt = record[4].replace("<NON_REF>", record[3])
    return bool(len(record[3]) == 1 and len(alt) - alt.count(",") == alt.count(",") + 1)


def num_genotypes(record, num_samples):
    """
    Get number of genotypes in VCF record, total number of samples - missing genotypes
    """
    missing = 0
    for i in range(9, num_samples + 9):
        if record[i].startswith("."):
            missing += 1
    return num_samples - missing


def get_matrix_column(record, num_samples, resolve_IUPAC):
    """
    Transform a VCF record into a phylogenetic matrix column with nucleotides instead of numbers
    """
    nt_dict = {str(0): record[3].replace("-","*").upper(), ".": "N"}
    # <NON_REF> must be replaced by the REF in the ALT field for GVCFs from GATK
    alt = record[4].replace("-", "*").replace("<NON_REF>", nt_dict["0"])
    alt = alt.split(",")
    for n in range(len(alt)):
        nt_dict[str(n+1)] = alt[n]
    column = ""
    for i in range(9, num_samples + 9):
        geno_num = record[i].split(":")[0].replace("/", "").replace("|", "")
        try:
            geno_nuc = "".join(sorted(set([nt_dict[j] for j in geno_num])))
        except KeyError:
            return "malformed"
        if resolve_IUPAC is False:
            column += AMBIG[geno_nuc]
        else:
            column += AMBIG[nt_dict[random.choice(geno_num)]]
    return column


def get_matrix_column_bin(record, num_samples):
    """
    Return an alignment column in NEXUS binary from a VCF record, if genotype is not diploid with at
    most two alleles it will return '?' as state
    """
    column = ""
    for i in range(9, num_samples + 9):
        genotype = record[i].split(":")[0]
        if genotype in GEN_BIN:
            column += GEN_BIN[genotype]
        else:
            column += "?"
    return column


# ──────────────────────────────────────────────────────────────────────────────
# Worker functions for multiprocessing
# ──────────────────────────────────────────────────────────────────────────────

def _process_lines_chunk(args):
    """
    Process a batch of VCF data lines (already stripped/split).
    Returns list of tuples: (site_tmp_or_None, binsite_tmp_or_None, is_accepted,
                              is_shallow, is_mnp, is_biallelic, chrom, pos,
                              num_samples_locus, line_text_for_error)
    """
    (lines, num_samples, min_samples_locus, resolve_IUPAC, write_nucleotide,
     write_binary) = args

    results = []
    for line in lines:
        record = line.split("\t")

        if is_anomalous(record, num_samples):
            results.append((None, None, False, False, False, False, None, None, 0, line))
            continue

        num_samples_locus = num_genotypes(record, num_samples)
        if num_samples_locus < min_samples_locus:
            results.append((None, None, False, True, False, False, None, None, num_samples_locus, None))
            continue

        if not is_snp(record):
            results.append((None, None, False, False, True, False, None, None, num_samples_locus, None))
            continue

        # Passed filters — it's an accepted SNP
        site_tmp = None
        if write_nucleotide:
            site_tmp = get_matrix_column(record, num_samples, resolve_IUPAC)
            if site_tmp == "malformed":
                results.append((None, None, False, False, False, False, None, None, 0, line))
                continue

        binsite_tmp = None
        is_biallelic = False
        if write_binary and len(record[4]) == 1:
            binsite_tmp = get_matrix_column_bin(record, num_samples)
            is_biallelic = True

        chrom = record[0]
        pos = record[1]
        results.append((site_tmp, binsite_tmp, True, False, False, is_biallelic,
                         chrom, pos, num_samples_locus, None))

    return results


def _extract_sample_sequences(args):
    """
    Worker for Phase 2: read the temp file and extract characters at given sample indices.
    Returns dict: {sample_index: sequence_string}
    """
    tmp_path, sample_indices = args
    seqs = {idx: [] for idx in sample_indices}
    with open(tmp_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            for idx in sample_indices:
                if idx < len(line):
                    seqs[idx].append(line[idx])
    return {idx: "".join(chars) for idx, chars in seqs.items()}


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
        action = "store",
        dest = "threads",
        type = int,
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

    # Get samples names and number of samples in VCF
    if Path(args.filename).exists():
        sample_names = extract_sample_names(args.filename)
    else:
        print("\nInput VCF file not found, please verify the provided path")
        sys.exit()
    num_samples = len(sample_names)
    if num_samples == 0:
        print("\nSample names not found in VCF, your file may be corrupt or missing the header.\n")
        sys.exit()
    print("\nConverting file '{}':\n".format(args.filename))
    print("Number of samples in VCF: {:d}".format(num_samples))
    if args.threads > 1:
        print("Using {:d} parallel processes".format(args.threads))

    # If the 'min_samples_locus' is larger than the actual number of samples in VCF readjust it
    args.min_samples_locus = min(num_samples, args.min_samples_locus)

    # Output filename will be the same as input file, indicating the minimum of samples specified
    if not args.prefix:
        parts = Path(args.filename).name.split(".")
        args.prefix = []
        for p in parts:
            if p.lower() == "vcf":
                break
            else:
                args.prefix.append(p)
        args.prefix = ".".join(args.prefix)
    args.prefix += ".min" + str(args.min_samples_locus)

    # Check if outfolder exists, create it if it doesn't
    if not Path(args.folder).exists():
        Path(args.folder).mkdir(parents=True)

    outfile = str(Path(args.folder, args.prefix))

    # Flags for which outputs are needed
    write_nucleotide = args.fasta or args.nexus or not args.phylipdisable
    write_binary = args.nexusbin


    ##########################
    # PROCESS GENOTYPES IN VCF (Phase 1)

    # We need to create an intermediate file to hold the sequence data vertically and then transpose
    # it to create the matrices
    if write_nucleotide:
        temporal = open(outfile+".tmp", "w")

    # If binary NEXUS is selected also create a separate temporal
    if write_binary:
        temporalbin = open(outfile+".bin.tmp", "w")

    if args.write_used:
        used_sites = open(outfile+".used_sites.tsv", "w")
        used_sites.write("#CHROM\tPOS\tNUM_SAMPLES\n")

    if args.filename.lower().endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    # Counters
    snp_num = 0
    snp_accepted = 0
    snp_shallow = 0
    mnp_num = 0
    snp_biallelic = 0

    BATCH_SIZE = 10000  # Lines per batch for parallel processing

    if args.threads <= 1:
        # ── Original single-threaded path ──
        with opener(args.filename, "rt") as vcf:
            while 1:
                vcf_chunk = vcf.readlines(50000)
                if not vcf_chunk:
                    break

                for line in vcf_chunk:
                    line = line.strip()

                    if line and not line.startswith("#"):
                        record = line.split("\t")
                        snp_num += 1
                        if snp_num % 500000 == 0:
                            print("{:d} genotypes processed.".format(snp_num))
                        if is_anomalous(record, num_samples):
                            print("Skipping malformed line:\n{}".format(line))
                            continue
                        else:
                            num_samples_locus = num_genotypes(record, num_samples)
                            if num_samples_locus < args.min_samples_locus:
                                snp_shallow += 1
                                continue
                            else:
                                if is_snp(record):
                                    if write_nucleotide:
                                        site_tmp = get_matrix_column(record, num_samples,
                                                                     args.resolve_IUPAC)
                                        if site_tmp == "malformed":
                                            print("Skipping malformed line:\n{}".format(line))
                                            continue
                                        else:
                                            snp_accepted += 1
                                            temporal.write(site_tmp+"\n")
                                            if args.write_used:
                                                used_sites.write(record[0] + "\t"
                                                                 + record[1] + "\t"
                                                                 + str(num_samples_locus) + "\n")
                                    if write_binary:
                                        if len(record[4]) == 1:
                                            snp_biallelic += 1
                                            binsite_tmp = get_matrix_column_bin(record, num_samples)
                                            temporalbin.write(binsite_tmp+"\n")
                                else:
                                    mnp_num += 1

    else:
        # ── Parallel path ──
        with opener(args.filename, "rt") as vcf:
            batch = []
            futures = []
            write_nucleotide_flag = write_nucleotide
            write_binary_flag = write_binary

            with ProcessPoolExecutor(max_workers=args.threads) as executor:
                # Read and submit batches
                while True:
                    vcf_chunk = vcf.readlines(50000)
                    if not vcf_chunk:
                        break

                    for line in vcf_chunk:
                        line = line.strip()
                        if line and not line.startswith("#"):
                            batch.append(line)

                            if len(batch) >= BATCH_SIZE:
                                futures.append(executor.submit(
                                    _process_lines_chunk,
                                    (batch, num_samples, args.min_samples_locus,
                                     args.resolve_IUPAC, write_nucleotide_flag,
                                     write_binary_flag)
                                ))
                                batch = []

                # Submit remaining lines
                if batch:
                    futures.append(executor.submit(
                        _process_lines_chunk,
                        (batch, num_samples, args.min_samples_locus,
                         args.resolve_IUPAC, write_nucleotide_flag,
                         write_binary_flag)
                    ))

                # Collect results IN ORDER (futures are in submission order)
                for future in futures:
                    results = future.result()
                    for (site_tmp, binsite_tmp, accepted, shallow, mnp, biallelic,
                         chrom, pos, num_samples_locus, error_line) in results:

                        snp_num += 1
                        if snp_num % 500000 == 0:
                            print("{:d} genotypes processed.".format(snp_num))

                        if error_line is not None:
                            print("Skipping malformed line:\n{}".format(error_line))
                            continue

                        if shallow:
                            snp_shallow += 1
                            continue

                        if mnp:
                            mnp_num += 1
                            continue

                        if accepted:
                            if site_tmp is not None and write_nucleotide:
                                snp_accepted += 1
                                temporal.write(site_tmp + "\n")
                                if args.write_used:
                                    used_sites.write(chrom + "\t" + pos + "\t"
                                                     + str(num_samples_locus) + "\n")
                            if binsite_tmp is not None and write_binary:
                                snp_biallelic += 1
                                temporalbin.write(binsite_tmp + "\n")

    # Print useful information about filtering of SNPs
    print("Total of genotypes processed: {:d}".format(snp_num))
    print("Genotypes excluded because they exceeded the amount "
          "of missing data allowed: {:d}".format(snp_shallow))
    print("Genotypes that passed missing data filter but were "
          "excluded for being MNPs: {:d}".format(mnp_num))
    print("SNPs that passed the filters: {:d}".format(snp_accepted))
    if write_binary:
        print("Biallelic SNPs selected for binary NEXUS: {:d}".format(snp_biallelic))

    if args.write_used:
        print("Used sites saved to: '" + outfile + ".used_sites.tsv'")
        used_sites.close()
    print("")

    if write_nucleotide:
        temporal.close()
    if write_binary:
        temporalbin.close()


    #######################
    # WRITE OUTPUT MATRICES (Phase 2)

    if not args.phylipdisable:
        output_phy = open(outfile+".phy", "w")
        output_phy.write("{:d} {:d}\n".format(len(sample_names), snp_accepted))

    if args.fasta:
        output_fas = open(outfile+".fasta", "w")

    if args.nexus:
        output_nex = open(outfile+".nexus", "w")
        output_nex.write("#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={:d} NCHAR={:d};\n\tFORMAT "
                         "DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n".format(len(sample_names),
                                                                                      snp_accepted))

    if write_binary:
        output_nexbin = open(outfile+".bin.nexus", "w")
        output_nexbin.write("#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={:d} NCHAR={:d};\n\tFORMAT "
                            "DATATYPE=SNP MISSING=? GAP=- ;\nMATRIX\n".format(len(sample_names),
                                                                                     snp_biallelic))

    # Get length of longest sequence name
    len_longest_name = 0
    for name in sample_names:
        if len(name) > len_longest_name:
            len_longest_name = len(name)

    # Write outgroup as first sequence in alignment if the name is specified
    idx_outgroup = None
    if outgroup in sample_names:
        idx_outgroup = sample_names.index(outgroup)

    if args.threads <= 1 or (not write_nucleotide and not write_binary):
        # ── Original single-threaded transpose path ──
        if idx_outgroup is not None:
            if write_nucleotide:
                with open(outfile+".tmp") as tmp_seq:
                    seqout = ""
                    for line in tmp_seq:
                        seqout += line[idx_outgroup]
                if args.fasta:
                    output_fas.write(">"+sample_names[idx_outgroup]+"\n"+seqout+"\n")
                padding = (len_longest_name + 3 - len(sample_names[idx_outgroup])) * " "
                if not args.phylipdisable:
                    output_phy.write(sample_names[idx_outgroup]+padding+seqout+"\n")
                if args.nexus:
                    output_nex.write(sample_names[idx_outgroup]+padding+seqout+"\n")
                print("Outgroup, '{}', added to the matrix(ces).".format(outgroup))

            if write_binary:
                with open(outfile+".bin.tmp") as bin_tmp_seq:
                    seqout = ""
                    for line in bin_tmp_seq:
                        seqout += line[idx_outgroup]
                padding = (len_longest_name + 3 - len(sample_names[idx_outgroup])) * " "
                output_nexbin.write(sample_names[idx_outgroup]+padding+seqout+"\n")
                print("Outgroup, '{}', added to the binary matrix.".format(outgroup))

        for s in range(0, len(sample_names)):
            if s != idx_outgroup:
                if write_nucleotide:
                    with open(outfile+".tmp") as tmp_seq:
                        seqout = ""
                        for line in tmp_seq:
                            seqout += line[s]
                    if args.fasta:
                        output_fas.write(">"+sample_names[s]+"\n"+seqout+"\n")
                    padding = (len_longest_name + 3 - len(sample_names[s])) * " "
                    if not args.phylipdisable:
                        output_phy.write(sample_names[s]+padding+seqout+"\n")
                    if args.nexus:
                        output_nex.write(sample_names[s]+padding+seqout+"\n")
                    print("Sample {:d} of {:d}, '{}', added to the nucleotide matrix(ces).".format(
                                                           s+1, len(sample_names), sample_names[s]))

                if write_binary:
                    with open(outfile+".bin.tmp") as bin_tmp_seq:
                        seqout = ""
                        for line in bin_tmp_seq:
                            seqout += line[s]
                    padding = (len_longest_name + 3 - len(sample_names[s])) * " "
                    output_nexbin.write(sample_names[s]+padding+seqout+"\n")
                    print("Sample {:d} of {:d}, '{}', added to the binary matrix.".format(
                                                           s+1, len(sample_names), sample_names[s]))

    else:
        # ── Parallel transpose path ──
        # Build ordered list of sample indices to process
        sample_order = []
        if idx_outgroup is not None:
            sample_order.append(idx_outgroup)
        for s in range(len(sample_names)):
            if s != idx_outgroup:
                sample_order.append(s)

        def write_sample_seq(sample_idx, seqout, is_outgroup=False):
            """Write a single sample's sequence to all enabled output files."""
            tag = "Outgroup" if is_outgroup else "Sample"
            if write_nucleotide:
                if args.fasta:
                    output_fas.write(">"+sample_names[sample_idx]+"\n"+seqout+"\n")
                padding = (len_longest_name + 3 - len(sample_names[sample_idx])) * " "
                if not args.phylipdisable:
                    output_phy.write(sample_names[sample_idx]+padding+seqout+"\n")
                if args.nexus:
                    output_nex.write(sample_names[sample_idx]+padding+seqout+"\n")
                if is_outgroup:
                    print("Outgroup, '{}', added to the matrix(ces).".format(outgroup))
                else:
                    print("Sample, '{}', added to the nucleotide matrix(ces).".format(
                              sample_names[sample_idx]))
            if write_binary:
                padding = (len_longest_name + 3 - len(sample_names[sample_idx])) * " "
                output_nexbin.write(sample_names[sample_idx]+padding+seqout+"\n")
                if is_outgroup:
                    print("Outgroup, '{}', added to the binary matrix.".format(outgroup))
                else:
                    print("Sample, '{}', added to the binary matrix.".format(
                              sample_names[sample_idx]))

        # Read temp files once and extract all sequences in parallel
        num_workers = min(args.threads, len(sample_order))

        # ── Nucleotide matrix ──
        if write_nucleotide and Path(outfile+".tmp").exists():
            print("Transposing nucleotide matrix with {:d} workers...".format(num_workers))

            # Split sample_order into chunks for workers
            chunks = [[] for _ in range(num_workers)]
            for i, idx in enumerate(sample_order):
                chunks[i % num_workers].append(idx)

            # Filter out empty chunks
            chunks = [c for c in chunks if c]

            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                futures_map = {}
                for chunk in chunks:
                    future = executor.submit(_extract_sample_sequences,
                                             (outfile+".tmp", chunk))
                    futures_map[future] = chunk

                # Collect all sequences
                all_seqs = {}
                for future in as_completed(futures_map):
                    result = future.result()
                    all_seqs.update(result)

            # Write in original order
            for idx in sample_order:
                is_outgroup = (idx == idx_outgroup)
                seqout = all_seqs.get(idx, "")
                if write_nucleotide:
                    if args.fasta:
                        output_fas.write(">"+sample_names[idx]+"\n"+seqout+"\n")
                    padding = (len_longest_name + 3 - len(sample_names[idx])) * " "
                    if not args.phylipdisable:
                        output_phy.write(sample_names[idx]+padding+seqout+"\n")
                    if args.nexus:
                        output_nex.write(sample_names[idx]+padding+seqout+"\n")
                    if is_outgroup:
                        print("Outgroup, '{}', added to the matrix(ces).".format(outgroup))
                    else:
                        print("Sample, '{}', added to the nucleotide matrix(ces).".format(
                                  sample_names[idx]))

        # ── Binary matrix ──
        if write_binary and Path(outfile+".bin.tmp").exists():
            print("Transposing binary matrix with {:d} workers...".format(num_workers))

            chunks_bin = [[] for _ in range(num_workers)]
            for i, idx in enumerate(sample_order):
                chunks_bin[i % num_workers].append(idx)
            chunks_bin = [c for c in chunks_bin if c]

            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                futures_map = {}
                for chunk in chunks_bin:
                    future = executor.submit(_extract_sample_sequences,
                                             (outfile+".bin.tmp", chunk))
                    futures_map[future] = chunk

                all_seqs_bin = {}
                for future in as_completed(futures_map):
                    result = future.result()
                    all_seqs_bin.update(result)

            for idx in sample_order:
                is_outgroup = (idx == idx_outgroup)
                seqout = all_seqs_bin.get(idx, "")
                padding = (len_longest_name + 3 - len(sample_names[idx])) * " "
                output_nexbin.write(sample_names[idx]+padding+seqout+"\n")
                if is_outgroup:
                    print("Outgroup, '{}', added to the binary matrix.".format(outgroup))
                else:
                    print("Sample, '{}', added to the binary matrix.".format(
                              sample_names[idx]))

    print()
    if not args.phylipdisable:
        print("PHYLIP matrix saved to: " + outfile+".phy")
        output_phy.close()
    if args.fasta:
        print("FASTA matrix saved to: " + outfile+".fasta")
        output_fas.close()
    if args.nexus:
        output_nex.write(";\nEND;\n")
        print("NEXUS matrix saved to: " + outfile+".nex")
        output_nex.close()
    if write_binary:
        output_nexbin.write(";\nEND;\n")
        print("BINARY NEXUS matrix saved to: " + outfile+".bin.nex")
        output_nexbin.close()

    if write_nucleotide:
        Path(outfile+".tmp").unlink()
    if write_binary:
        Path(outfile+".bin.tmp").unlink()

    print( "\nDone!\n")


if __name__ == "__main__":
    main()
