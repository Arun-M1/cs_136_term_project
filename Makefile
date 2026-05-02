CC = gcc
CFLAGS = -Wall -O2
LIBS = -lm

canny_edge: Canny_Edge_Detection/canny_edge.c netpbm.c
	$(CC) $(CFLAGS) -I. -o demo/canny_edge Canny_Edge_Detection/canny_edge.c netpbm.c $(LIBS)

demo:
	bash demo/run_demo.sh

clean:
	rm -f demo/canny_edge
