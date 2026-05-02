#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "netpbm.h"

Matrix smoothing_filter(Matrix m1, Matrix m2);
Matrix createGaussianFilter();

Matrix smoothing_filter(Matrix m1, Matrix m2) {
    Matrix resultMatrix = createMatrix(m1.height, m1.width);
    int centerR = (m2.height - 1) / 2;
    int centerC = (m2.width - 1) / 2;

    for (int r = 0; r < m1.height; r++) {
        for (int c = 0; c < m1.width; c++) {
            int top = r - centerR;
            int bottom = r + (m2.height - centerR - 1);
            int left = c - centerC;
            int right = c + (m2.width - centerC - 1);

            //check bounds
            if (top < 0 || bottom >= m1.height || left < 0 || right >= m1.width) {
                resultMatrix.map[r][c] = 0.0;
                continue;
            }

            double sum = 0.0;
            //convolution
            for (int fr = 0; fr < m2.height; fr++) {
                for (int fc = 0; fc < m2.width; fc++) {
                    int image_row = r + (fr - centerR);
                    int image_col = c + (fc - centerC);
                    sum += m1.map[image_row][image_col] * m2.map[fr][fc];
                }
            }
            resultMatrix.map[r][c] = sum;
        }
    }

    return resultMatrix;
}

Matrix createGaussianFilter() {
    Matrix filter = createMatrix(3, 3);
    double values[3][3] = {
        {1, 2, 1},
        {2, 4, 2},
        {1, 2, 1}
    };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            filter.map[i][j] = values[i][j] / 16.0;  // Normalize so sum = 1
        }
    }
    return filter;
}


int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s input.pgm output.pgm\n", argv[0]);
        return 1;
    }

    Image img = readImage(argv[1]);
    Matrix imageMatrix = image2Matrix(img);
    Matrix gaussianFilter = createGaussianFilter();
    Matrix smoothedImageMatrix = smoothing_filter(imageMatrix, gaussianFilter);
    Image smoothedImage = matrix2Image(smoothedImageMatrix, 1, 1.0);

    writeImage(smoothedImage, argv[2]);

    deleteImage(img);
    deleteMatrix(imageMatrix);
    deleteMatrix(gaussianFilter);
    deleteMatrix(smoothedImageMatrix);
    deleteImage(smoothedImage);
    return 0;
}
