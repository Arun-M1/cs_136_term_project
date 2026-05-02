#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$SCRIPT_DIR/.."
BINARY="$SCRIPT_DIR/canny_edge"
INPUT_DIR="$ROOT/Dataset_pgm"
OUTPUT_DIR="$SCRIPT_DIR/output"

# Build binary if not present
if [ ! -f "$BINARY" ]; then
    echo "Compiling canny_edge..."
    gcc -Wall -O2 -I"$ROOT" -o "$BINARY" "$ROOT/Canny_Edge_Detection/canny_edge.c" "$ROOT/netpbm.c" -lm
fi

mkdir -p "$OUTPUT_DIR/GeologyPGM"
mkdir -p "$OUTPUT_DIR/MarineBiologyPGM"
mkdir -p "$OUTPUT_DIR/New_DatasetPGM"

run() {
    local input="$1"
    local output="$2"
    echo "Processing: $input -> $output"
    "$BINARY" "$input" "$output"
}

# GeologyPGM
run "$INPUT_DIR/GeologyPGM/coastline-1.pgm"   "$OUTPUT_DIR/GeologyPGM/coastline-1_canny.pgm"
run "$INPUT_DIR/GeologyPGM/tahoe.pgm"          "$OUTPUT_DIR/GeologyPGM/tahoe_canny.pgm"

# MarineBiologyPGM
run "$INPUT_DIR/MarineBiologyPGM/32_T2_6_timepoint0.pgm"       "$OUTPUT_DIR/MarineBiologyPGM/32_T2_6_timepoint0_canny.pgm"
run "$INPUT_DIR/MarineBiologyPGM/32_T2_6_timepoint1.pgm"       "$OUTPUT_DIR/MarineBiologyPGM/32_T2_6_timepoint1_canny.pgm"
run "$INPUT_DIR/MarineBiologyPGM/33.5_T1_4_timepoint0.pgm"     "$OUTPUT_DIR/MarineBiologyPGM/33.5_T1_4_timepoint0_canny.pgm"
run "$INPUT_DIR/MarineBiologyPGM/33.5_T1_4_timepoint1.pgm"     "$OUTPUT_DIR/MarineBiologyPGM/33.5_T1_4_timepoint1_canny.pgm"
run "$INPUT_DIR/MarineBiologyPGM/35_T2_5_timepoint0.pgm"       "$OUTPUT_DIR/MarineBiologyPGM/35_T2_5_timepoint0_canny.pgm"
run "$INPUT_DIR/MarineBiologyPGM/35_T2_5_timepoint1.pgm"       "$OUTPUT_DIR/MarineBiologyPGM/35_T2_5_timepoint1_canny.pgm"

# New_DatasetPGM
run "$INPUT_DIR/New_DatasetPGM/DW.pgm"                                                "$OUTPUT_DIR/New_DatasetPGM/DW_canny.pgm"
run "$INPUT_DIR/New_DatasetPGM/T-101_DHEL-11__DHEL-12_20240802_BROOD.pgm"            "$OUTPUT_DIR/New_DatasetPGM/T-101_DHEL-11__DHEL-12_20240802_BROOD_canny.pgm"

echo "Done. Output files in $OUTPUT_DIR"
