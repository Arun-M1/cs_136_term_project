#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../netpbm.h"

Matrix convolve(Matrix m1, Matrix m2);
Image sobel(Image img);

Matrix convolve(Matrix m1, Matrix m2) {
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

Image sobel(Image img) {
    Matrix horizontalFilter = createMatrix(3, 3);
    double horizontalValues[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            horizontalFilter.map[i][j] = horizontalValues[i][j];
        }
    }

    Matrix verticalFilter = createMatrix(3, 3);
    double verticalValues[3][3] = {
        {-1, -2, -1},
        {0, 0, 0},
        {1, 2, 1}
    };

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            verticalFilter.map[i][j] = verticalValues[i][j];
        }
    }

    Matrix imgM = image2Matrix(img);
    Matrix horizontal = convolve(imgM, horizontalFilter);
    Matrix vertical = convolve(imgM, verticalFilter);

    Matrix strength = createMatrix(img.height, img.width);
    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            strength.map[r][c] = sqrt(pow(horizontal.map[r][c], 2) + pow(vertical.map[r][c], 2));
        }
    }

    Image res = matrix2Image(strength, 1, 1.0);

    deleteMatrix(horizontalFilter);
    deleteMatrix(verticalFilter);
    deleteMatrix(imgM);
    deleteMatrix(horizontal);
    deleteMatrix(vertical);
    deleteMatrix(strength);
    
    return res;
}


int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s input.pgm output.pgm\n", argv[0]);
        return 1;
    }

    Image img = readImage(argv[1]);
    Image result = sobel(img);
    writeImage(result, argv[2]);
    deleteImage(img);
    deleteImage(result);
    return 0;
}
