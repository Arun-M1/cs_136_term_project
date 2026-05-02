Compile:
clang sobel_filter.c ../netpbm.c -I.. -o sobel_filter -lm

Run:
./sobel_filter input.pgm output/output_sobel.pgm
