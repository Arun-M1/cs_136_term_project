#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR/.."
BINARY="$SCRIPT_DIR/fom"
PRED_DIR="$SCRIPT_DIR/output_bsds"
GT_DIR="$ROOT/Dataset_pgm/BSDS500/ground_truth"

# Build binary if not present
if [ ! -f "$BINARY" ]; then
    echo "Compiling fom..."
    gcc -Wall -O2 -I"$ROOT" -o "$BINARY" "$ROOT/Figure_of_Merit/fom.c" "$ROOT/netpbm.c" -lm
fi

IMAGES="2018 3063 5096 6046 8068"
total=0
count=0

for id in $IMAGES; do
    pred="$PRED_DIR/${id}_edges.pgm"
    gt="$GT_DIR/${id}_edges.pgm"
    result=$("$BINARY" "$pred" "$gt")
    score=$(echo "$result" | awk '{print $2}')
    printf "Image %s: %s\n" "$id" "$result"
    total=$(awk -v t="$total" -v s="$score" 'BEGIN {printf "%.6f", t+s}')
    count=$((count + 1))
done

mean=$(awk -v t="$total" -v n="$count" 'BEGIN {printf "%.6f", t/n}')
echo "---"
echo "Mean FOM: $mean"
