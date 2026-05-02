Compile:
clang hough_transform.c ../netpbm.c -I.. -o hough_transform -lm

Run:
./hough_transform input.pgm output/edges.pgm output/hough.pgm output/circles.ppm
