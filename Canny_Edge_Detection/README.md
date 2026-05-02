Compile:
clang canny_edge.c ../netpbm.c -I.. -o canny_edge -lm

Run:
./canny_edge input.pgm output/canny.pgm
