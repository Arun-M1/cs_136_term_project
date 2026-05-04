#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../netpbm.h"

void make_low_contrast(Image img, double factor) {
    int i, j;

    for (i = 0; i < img.height; i++) {
        for (j = 0; j < img.width; j++) {
            
            int r = (int)((img.map[i][j].r - 128) * factor + 128);
            int g = (int)((img.map[i][j].g - 128) * factor + 128);
            int b = (int)((img.map[i][j].b - 128) * factor + 128);
            int intensity = (int)((img.map[i][j].i - 128) * factor + 128);

            setPixel(img, i, j, r, g, b, intensity);
        }
    }
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s input.pgm output.pgm\n", argv[0]);
        return 1;
    }

    Image img = readImage(argv[1]);
    make_low_contrast(img, 0.25);

    writeImage(img, argv[2]);

    deleteImage(img);
    return 0;
}
