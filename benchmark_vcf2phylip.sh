#!/bin/bash
# ============================================================
# Benchmark: vcf2phylip original (1 thread) vs parallel (4/8/16)
# ============================================================

set -e

VCF="/home/vensin/workspace/Sweet_Osmanthus/05.variant_filter/03.LD_prune_SNP_plink2/412samples.SNP.biallelic.minGQ10.minQ30.meanDP6.maxmiss0.8.maf0.05.LDpruned.vcf.gz"
ORIGINAL="/tmp/vcf2phylip_original.py"
PARALLEL="/home/vensin/software/vcf2phylip/vcf2phylip.py"
BASEDIR="/tmp/vcf2phylip_benchmark"

# Clean up previous runs
rm -rf "$BASEDIR"
mkdir -p "$BASEDIR"

echo "=============================================="
echo "  vcf2phylip Benchmark"
echo "=============================================="
echo "Input VCF: $(basename $VCF)"
echo "VCF size : $(du -h $VCF | cut -f1)"
echo "Date     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "CPU      : $(lscpu | grep 'Model name' | sed 's/.*:\s*//')"
echo "=============================================="
echo ""

run_benchmark() {
    local label="$1"
    local script="$2"
    local threads="$3"
    local outdir="$BASEDIR/$label"
    mkdir -p "$outdir"

    echo -n "Running: $label ... "

    local start_time=$(date +%s%N)

    python3 "$script" \
        -i "$VCF" \
        -f \
        -m 4 \
        --output-folder "$outdir" \
        --output-prefix "bench" \
        ${threads:+-t $threads} \
        > "$outdir/stdout.log" 2>&1

    local end_time=$(date +%s%N)
    local elapsed_ms=$(( (end_time - start_time) / 1000000 ))
    local elapsed_sec=$(echo "scale=2; $elapsed_ms / 1000" | bc)

    echo "${elapsed_sec}s"

    # Store result
    echo "$label|$elapsed_sec" >> "$BASEDIR/results.tsv"

    # Extract key stats
    local snp_passed=$(grep "SNPs that passed" "$outdir/stdout.log" | grep -oP '\d+$')
    echo "  SNPs passed: $snp_passed | Time: ${elapsed_sec}s"
    echo ""
}

echo "Starting benchmarks..."
echo ""

# 1. Original script, 1 thread
run_benchmark "original_t1" "$ORIGINAL" ""

# 2. Parallel script, 1 thread (baseline)
run_benchmark "parallel_t1" "$PARALLEL" "1"

# 3. Parallel script, 4 threads
run_benchmark "parallel_t4" "$PARALLEL" "4"

# 4. Parallel script, 8 threads
run_benchmark "parallel_t8" "$PARALLEL" "8"

# 5. Parallel script, 16 threads
run_benchmark "parallel_t16" "$PARALLEL" "16"

# ── Verify correctness ──
echo "=============================================="
echo "  Correctness Check"
echo "=============================================="

phy_orig="$BASEDIR/original_t1/bench.min4.phy"
phy_t1="$BASEDIR/parallel_t1/bench.min4.phy"
phy_t4="$BASEDIR/parallel_t4/bench.min4.phy"
phy_t8="$BASEDIR/parallel_t8/bench.min4.phy"
phy_t16="$BASEDIR/parallel_t16/bench.min4.phy"

all_match=true
for f in "$phy_t1" "$phy_t4" "$phy_t8" "$phy_t16"; do
    if ! diff -q "$phy_orig" "$f" > /dev/null 2>&1; then
        echo "MISMATCH: $(basename $f) differs from original!"
        all_match=false
    fi
done
if $all_match; then
    echo "All PHYLIP outputs are IDENTICAL ✓"
fi

fast_orig="$BASEDIR/original_t1/bench.min4.fasta"
fast_t1="$BASEDIR/parallel_t1/bench.min4.fasta"
fast_t4="$BASEDIR/parallel_t4/bench.min4.fasta"
fast_t8="$BASEDIR/parallel_t8/bench.min4.fasta"
fast_t16="$BASEDIR/parallel_t16/bench.min4.fasta"

all_match_fasta=true
for f in "$fast_t1" "$fast_t4" "$fast_t8" "$fast_t16"; do
    if ! diff -q "$fast_orig" "$f" > /dev/null 2>&1; then
        echo "MISMATCH: $(basename $f) FASTA differs from original!"
        all_match_fasta=false
    fi
done
if $all_match_fasta; then
    echo "All FASTA outputs are IDENTICAL ✓"
fi
echo ""

# ── Summary table ──
echo "=============================================="
echo "  Results Summary"
echo "=============================================="
printf "%-20s %10s %10s\n" "Configuration" "Time(s)" "Speedup"
printf "%-20s %10s %10s\n" "─────────────────" "───────" "───────"

baseline=$(grep "original_t1" "$BASEDIR/results.tsv" | cut -d'|' -f2)

while IFS='|' read -r label time_sec; do
    speedup=$(echo "scale=2; $baseline / $time_sec" | bc)
    printf "%-20s %10s %10sx\n" "$label" "$time_sec" "$speedup"
done < "$BASEDIR/results.tsv"

echo ""
echo "Temp dir: $BASEDIR"
echo "Done!"
