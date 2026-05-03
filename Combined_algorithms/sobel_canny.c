#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "netpbm.h"

#define SOBEL_THRESHOLD 130.0
#define CANNY_LOW_RATIO 0.10
#define CANNY_HIGH_RATIO 0.22

Matrix convolve(Matrix m1, Matrix m2);
Matrix createGaussianFilter();
Matrix sobel(Image img);
Matrix cannyEdges(Image img);
Matrix combineSobelCanny(Matrix sobelMatrix, Matrix cannyMatrix);

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

            if (top < 0 || bottom >= m1.height || left < 0 || right >= m1.width) {
                resultMatrix.map[r][c] = 0.0;
                continue;
            }

            double sum = 0.0;
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
            filter.map[i][j] = values[i][j] / 16.0;
        }
    }
    return filter;
}

Matrix sobel(Image img) {
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

    deleteMatrix(horizontalFilter);
    deleteMatrix(verticalFilter);
    deleteMatrix(imgM);
    deleteMatrix(horizontal);
    deleteMatrix(vertical);

    return strength;
}

Matrix canny(Image img) {
    Matrix matrixImg = image2Matrix(img);
    Matrix gaussianFilter = createGaussianFilter();
    Matrix smoothedMatrix = convolve(matrixImg, gaussianFilter);

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

    Matrix gX = convolve(smoothedMatrix, horizontalFilter);
    Matrix gY = convolve(smoothedMatrix, verticalFilter);

    Matrix magnitude = createMatrix(img.height, img.width);
    Matrix direction = createMatrix(img.height, img.width);

    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            double gx = gX.map[r][c];
            double gy = gY.map[r][c];

            magnitude.map[r][c] = sqrt(gx * gx + gy * gy);
            direction.map[r][c] = atan2(gy, gx);
        }
    }

    Matrix sector = createMatrix(img.height, img.width);
    Matrix suppressedMatrix = createMatrix(img.height, img.width);

    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            double angle = direction.map[r][c] * 180.0 / PI;

            if (angle < 0) {
                angle += 180.0;
            }

            if ((angle >= 0.0 && angle < 22.5) || (angle >= 157.5 && angle < 180.0)) {
                sector.map[r][c] = 0;
            }
            else if (angle >= 22.5 && angle < 67.5) {
                sector.map[r][c] = 1;
            }
            else if (angle >= 67.5 && angle < 112.5) {
                sector.map[r][c] = 2;
            }
            else {
                sector.map[r][c] = 3;
            }
        }
    }

    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            double mag = magnitude.map[r][c];
            int s = (int)sector.map[r][c];

            double neighbor1 = 0.0;
            double neighbor2 = 0.0;

            if (s == 0) {
                if (c - 1 >= 0) {
                    neighbor1 = magnitude.map[r][c - 1];
                }
                if (c + 1 < img.width) {
                    neighbor2 = magnitude.map[r][c + 1];
                }
            }
            else if (s == 1) {
                if (r - 1 >= 0 && c + 1 < img.width) {
                    neighbor1 = magnitude.map[r - 1][c + 1];
                }
                if (r + 1 < img.height && c - 1 >= 0) {
                    neighbor2 = magnitude.map[r + 1][c - 1];
                }
            }
            else if (s == 2) {
                if (r - 1 >= 0) {
                    neighbor1 = magnitude.map[r - 1][c];
                }
                if (r + 1 < img.height) {
                    neighbor2 = magnitude.map[r + 1][c];
                }
            }
            else {
                if (r - 1 >= 0 && c - 1 >= 0) {
                    neighbor1 = magnitude.map[r - 1][c - 1];
                }
                if (r + 1 < img.height && c + 1 < img.width) {
                    neighbor2 = magnitude.map[r + 1][c + 1];
                }
            }

            if (mag >= neighbor1 && mag >= neighbor2) {
                suppressedMatrix.map[r][c] = mag;
            }
            else {
                suppressedMatrix.map[r][c] = 0.0;
            }
        }
    }

    double maxMag = 0;
    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            if (suppressedMatrix.map[r][c] > maxMag) {
                maxMag = suppressedMatrix.map[r][c];
            }
        }
    }

    double lowThreshold = CANNY_LOW_RATIO * maxMag;
    double highThreshold = CANNY_HIGH_RATIO * maxMag;

    Matrix edges = createMatrix(img.height, img.width);
    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            double val = suppressedMatrix.map[r][c];

            if (val >= highThreshold) {
                edges.map[r][c] = 255;
            }
            else if (val <= lowThreshold) {
                edges.map[r][c] = 0;
            }
            else {
                edges.map[r][c] = 25;
            }
        }
    }

    int changed = 1;
    while (changed) {
        changed = 0;

        for (int r = 0; r < img.height; r++) {
            for (int c = 0; c < img.width; c++) {
                if (edges.map[r][c] == 25) {
                    int connectedToStrong = 0;

                    for (int dr = -1; dr <= 1; dr++) {
                        for (int dc = -1; dc <= 1; dc++) {
                            if (dr == 0 && dc == 0) {
                                continue;
                            }

                            int nr = r + dr;
                            int nc = c + dc;

                            if (nr >= 0 && nr < img.height && nc >= 0 && nc < img.width) {
                                if (edges.map[nr][nc] == 255) {
                                    connectedToStrong = 1;
                                }
                            }
                        }
                    }

                    if (connectedToStrong) {
                        edges.map[r][c] = 255;
                        changed = 1;
                    }
                }
            }
        }
    }

    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            if (edges.map[r][c] == 25) {
                edges.map[r][c] = 0;
            }
        }
    }

    deleteMatrix(matrixImg);
    deleteMatrix(gaussianFilter);
    deleteMatrix(smoothedMatrix);
    deleteMatrix(horizontalFilter);
    deleteMatrix(verticalFilter);
    deleteMatrix(gX);
    deleteMatrix(gY);
    deleteMatrix(magnitude);
    deleteMatrix(direction);
    deleteMatrix(sector);
    deleteMatrix(suppressedMatrix);

    return edges;
}

Matrix combineSobelCanny(Matrix sobelMatrix, Matrix cannyMatrix) {
    Matrix combined = createMatrix(sobelMatrix.height, sobelMatrix.width);

    for (int r = 0; r < sobelMatrix.height; r++) {
        for (int c = 0; c < sobelMatrix.width; c++) {
            if (cannyMatrix.map[r][c] > 0.0 && sobelMatrix.map[r][c] >= SOBEL_THRESHOLD) {
                combined.map[r][c] = 255.0;
            }
            else {
                combined.map[r][c] = 0.0;
            }
        }
    }

    return combined;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s input.pgm output.pgm\n", argv[0]);
        return 1;
    }

    Image img = readImage(argv[1]);

    Matrix sobelMatrix = sobel(img);

    Matrix cannyMatrix = canny(img);

    Matrix combinedMatrix = combineSobelCanny(sobelMatrix, cannyMatrix);

    Image result = matrix2Image(combinedMatrix, 0, 1.0);
    writeImage(result, argv[2]);

    deleteImage(img);
    deleteImage(result);
    deleteMatrix(sobelMatrix);
    deleteMatrix(cannyMatrix);
    deleteMatrix(combinedMatrix);

    return 0;
}
