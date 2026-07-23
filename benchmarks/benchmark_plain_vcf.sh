#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 /path/to/vcf2phylip.py /path/to/input.vcf [output_dir]" >&2
    exit 2
fi

SCRIPT=$(realpath "$1")
VCF=$(realpath "$2")
OUT_ROOT=${3:-plain_vcf_benchmark}
REPEATS=${REPEATS:-5}
THREAD_LIST=${THREAD_LIST:-"1 2 4 8 16"}

mkdir -p "$OUT_ROOT"
RESULTS="$OUT_ROOT/times.tsv"
printf 'backend\tthreads\trepeat\telapsed_s\tuser_s\tsystem_s\tmax_rss_kb\n' > "$RESULTS"

cleanup_dir=""
cleanup() {
    if [[ -n "$cleanup_dir" && -d "$cleanup_dir" ]]; then
        rm -rf "$cleanup_dir"
    fi
}
trap cleanup EXIT INT TERM

for backend in plain plain-stream; do
    for threads in $THREAD_LIST; do
        for repeat in $(seq 1 "$REPEATS"); do
            cleanup_dir="$OUT_ROOT/.tmp_${backend}_t${threads}_r${repeat}"
            rm -rf "$cleanup_dir"
            mkdir -p "$cleanup_dir"
            /usr/bin/time \
                -f "$backend\t$threads\t$repeat\t%e\t%U\t%S\t%M" \
                -o "$RESULTS" -a \
                python3 "$SCRIPT" \
                    -i "$VCF" \
                    -t "$threads" \
                    --input-backend "$backend" \
                    --output-folder "$cleanup_dir" \
                    >/dev/null
            rm -rf "$cleanup_dir"
            cleanup_dir=""
        done
    done
done

printf 'Results: %s\n' "$RESULTS"
