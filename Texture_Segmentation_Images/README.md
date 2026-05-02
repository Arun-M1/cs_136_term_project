Compile:
clang texture_segmentation.c ../netpbm.c -I.. -o texture_segmentation -lm

Run:
./texture_segmentation input.pgm output/texture_segments.ppm 4
