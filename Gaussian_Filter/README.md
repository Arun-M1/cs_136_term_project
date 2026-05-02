Compile:
clang gaussian_filter.c ../netpbm.c -I.. -o gaussian_filter -lm

Run:
./gaussian_filter input.pgm output/gaussian.pgm
