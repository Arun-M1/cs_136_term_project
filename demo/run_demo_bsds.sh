#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR/.."
BINARY="$SCRIPT_DIR/gaussian_canny_combo"
INPUT_DIR="$ROOT/Dataset_pgm/BSDS500/images"
OUTPUT_DIR="$SCRIPT_DIR/output_bsds"

# Build binary if not present
if [ ! -f "$BINARY" ]; then
    echo "Compiling gaussian_canny_combo..."
    gcc -Wall -O2 -I"$ROOT" -o "$BINARY" "$ROOT/Combined_algorithms/gaussian_plus_canny.c" "$ROOT/netpbm.c" -lm
fi

mkdir -p "$OUTPUT_DIR"

run() {
    local input="$1"
    local output="$2"
    echo "Processing: $input -> $output"
    "$BINARY" "$input" "$output"
}

run "$INPUT_DIR/2018.pgm"  "$OUTPUT_DIR/2018_edges.pgm"
run "$INPUT_DIR/3063.pgm"  "$OUTPUT_DIR/3063_edges.pgm"
run "$INPUT_DIR/5096.pgm"  "$OUTPUT_DIR/5096_edges.pgm"
run "$INPUT_DIR/6046.pgm"  "$OUTPUT_DIR/6046_edges.pgm"
run "$INPUT_DIR/8068.pgm"  "$OUTPUT_DIR/8068_edges.pgm"

echo "Done. Output files in $OUTPUT_DIR"
