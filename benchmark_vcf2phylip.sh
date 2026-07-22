#!/bin/bash
# ============================================================
# Benchmark: vcf2phylip original v2.9 vs merged mt2
# ============================================================

set -e

VCF="/home/vensin/workspace/Sweet_Osmanthus/05.variant_filter/03.LD_prune_SNP_plink2/412samples.SNP.biallelic.minGQ10.minQ30.meanDP6.maxmiss0.8.maf0.05.LDpruned.vcf.gz"
ORIGINAL="/tmp/vcf2phylip_orig_v29.py"
MERGED="/home/vensin/software/vcf2phylip/vcf2phylip.py"
BASEDIR="/tmp/vcf2phylip_bench"

# Download original v2.9 for comparison
curl -sL https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py -o "$ORIGINAL"

rm -rf "$BASEDIR"
mkdir -p "$BASEDIR"

echo "=============================================="
echo "  vcf2phylip Benchmark"
echo "=============================================="
echo "Input VCF : $(basename $VCF)"
echo "VCF size  : $(du -h $VCF | cut -f1)"
echo "Date      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "CPU       : $(lscpu | grep 'Model name' | sed 's/.*:\s*//')"
echo "CPU cores : $(nproc)"
echo "=============================================="
echo ""

run_benchmark() {
    local label="$1"
    local script="$2"
    local threads="$3"
    local outdir="$BASEDIR/$label"
    mkdir -p "$outdir"

    printf "%-22s" "$label"

    local start_ns=$(date +%s%N)

    python3 "$script" \
        -i "$VCF" \
        -f \
        -m 4 \
        --output-folder "$outdir" \
        --output-prefix "bench" \
        ${threads:+-t $threads} \
        > "$outdir/stdout.log" 2>&1

    local end_ns=$(date +%s%N)
    local elapsed_ms=$(( (end_ns - start_ns) / 1000000 ))
    local elapsed_sec=$(echo "scale=2; $elapsed_ms / 1000" | bc)

    echo "${elapsed_sec}s"
    echo "$label|$elapsed_sec" >> "$BASEDIR/results.tsv"
}

echo "Starting benchmarks..."
echo ""
printf "%-22s %10s\n" "Configuration" "Time"
printf "%-22s %10s\n" "─────────────────────" "──────────"

run_benchmark "original_v2.9"     "$ORIGINAL" ""
run_benchmark "merged_t1"         "$MERGED"   "1"
run_benchmark "merged_auto"       "$MERGED"   "0"
run_benchmark "merged_t4"         "$MERGED"   "4"
run_benchmark "merged_t8"         "$MERGED"   "8"
run_benchmark "merged_t16"        "$MERGED"   "16"

echo ""

# ── Verify correctness ──
echo "=============================================="
echo "  Correctness Check (PHYLIP)"
echo "=============================================="

phy_orig="$BASEDIR/original_v2.9/bench.min4.phy"
all_match=true
for label in merged_t1 merged_auto merged_t4 merged_t8 merged_t16; do
    f="$BASEDIR/$label/bench.min4.phy"
    if ! diff -q "$phy_orig" "$f" > /dev/null 2>&1; then
        echo "MISMATCH: $label differs from original!"
        all_match=false
    fi
done
if $all_match; then
    echo "All PHYLIP outputs are IDENTICAL ✓"
fi

all_match_fasta=true
fast_orig="$BASEDIR/original_v2.9/bench.min4.fasta"
for label in merged_t1 merged_auto merged_t4 merged_t8 merged_t16; do
    f="$BASEDIR/$label/bench.min4.fasta"
    if ! diff -q "$fast_orig" "$f" > /dev/null 2>&1; then
        echo "MISMATCH: $label FASTA differs from original!"
        all_match_fasta=false
    fi
done
if $all_match_fasta; then
    echo "All FASTA outputs are IDENTICAL ✓"
fi
echo ""

# ── Summary table with speedup ──
echo "=============================================="
echo "  Results Summary"
echo "=============================================="
printf "%-22s %10s %10s\n" "Configuration" "Time(s)" "Speedup"
printf "%-22s %10s %10s\n" "─────────────────────" "───────" "──────"

baseline=$(grep "original_v2.9" "$BASEDIR/results.tsv" | cut -d'|' -f2)

while IFS='|' read -r label time_sec; do
    speedup=$(echo "scale=2; $baseline / $time_sec" | bc)
    printf "%-22s %10s %10sx\n" "$label" "$time_sec" "$speedup"
done < "$BASEDIR/results.tsv"

echo ""
echo "Temp dir: $BASEDIR"
echo "Done!"
