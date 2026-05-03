CC = gcc
CFLAGS = -Wall -O2
LIBS = -lm

canny_edge: Canny_Edge_Detection/canny_edge.c netpbm.c
	$(CC) $(CFLAGS) -I. -o demo/canny_edge Canny_Edge_Detection/canny_edge.c netpbm.c $(LIBS)

gaussian_canny_combo: Combined_algorithms/gaussian_plus_canny.c netpbm.c
	$(CC) $(CFLAGS) -I. -o demo/gaussian_canny_combo Combined_algorithms/gaussian_plus_canny.c netpbm.c $(LIBS)

demo:
	bash demo/run_demo.sh

demo_bsds:
	bash demo/run_demo_bsds.sh

fom: Figure_of_Merit/fom.c netpbm.c
	$(CC) $(CFLAGS) -I. -o demo/fom Figure_of_Merit/fom.c netpbm.c $(LIBS)

demo_fom:
	bash demo/run_fom.sh

clean:
	rm -f demo/canny_edge demo/gaussian_canny_combo demo/fom
