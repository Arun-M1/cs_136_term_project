#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../netpbm.h"

Matrix smoothing_filter(Matrix m1, Matrix m2);
Matrix createGaussianFilter();
Matrix convolve(Matrix m1, Matrix m2);
Image canny(Image img);

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

Image canny(Image img) {
    //smooth with gaussian filter
    Matrix matrixImg = image2Matrix(img);
    Matrix gaussianFilter = createGaussianFilter();
    Matrix smoothedMatrix = smoothing_filter(matrixImg, gaussianFilter);

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


    //find magnitude and orientiation of gradient
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

    //non-maxima suppression
    Matrix sector = createMatrix(img.height, img.width);
    Matrix suppressedMatrix = createMatrix(img.height, img.width);
    
    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            double angle = direction.map[r][c] * 180.0 / M_PI;
            
            //check half the circle, get other half in negative direction
            if (angle < 0) {
                angle += 180.0;
            }

            //sector 0 check
            if ((angle >= 0.0 && angle < 22.5) || (angle >= 157.5 && angle < 180.0)) {
                sector.map[r][c] = 0;
            }//sector 1 check
            else if (angle >= 22.5 && angle < 67.5) {
                sector.map[r][c] = 1;
            }//sector 2 check
            else if (angle >= 67.5 && angle < 112.5) {
                sector.map[r][c] = 2;
            } 
            else {//sector 3
                sector.map[r][c] = 3;
            }
        }
    }

    //use sector matrix to get neighbors and compare
    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            double mag = magnitude.map[r][c];
            int s = (int)sector.map[r][c];

            double neighbor1 = 0.0;
            double neighbor2 = 0.0;

            if (s == 0) {
                //check in between bounds
                if (c - 1 >= 0) {
                    neighbor1 = magnitude.map[r][c - 1];
                } else {
                    neighbor1 = 0.0;
                }

                if (c + 1 < img.width) {
                    neighbor2 = magnitude.map[r][c + 1];
                } else {
                    neighbor2 = 0.0;
                }
            }
            else if (s == 1) {
                if (r - 1 >= 0 && c + 1 < img.width) {
                    neighbor1 = magnitude.map[r - 1][c + 1];
                } else {
                    neighbor1 = 0.0;
                }

                if (r + 1 < img.height && c - 1 >= 0) {
                    neighbor2 = magnitude.map[r + 1][c - 1];
                } else {
                    neighbor2 = 0.0;
                }
            }
            else if (s == 2) {
                if (r - 1 >= 0) {
                    neighbor1 = magnitude.map[r - 1][c];
                } else {
                    neighbor1 = 0.0;
                }

                if (r + 1 < img.height) {
                    neighbor2 = magnitude.map[r + 1][c];
                } else {
                    neighbor2 = 0.0;
                }
            }
            else {
                if (r - 1 >= 0 && c - 1 >= 0) {
                    neighbor1 = magnitude.map[r - 1][c - 1];
                } else {
                    neighbor1 = 0.0;
                }

                if (r + 1 < img.height && c + 1 < img.width) {
                    neighbor2 = magnitude.map[r + 1][c + 1];
                } else {
                    neighbor2 = 0.0;
                }
            }

            if (mag >= neighbor1 && mag >= neighbor2) {
                suppressedMatrix.map[r][c] = mag;
            } else {
                suppressedMatrix.map[r][c] = 0.0;
            }
        }
    }

    //hysterisis thresholding
    double maxMag = 0;
    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            if (suppressedMatrix.map[r][c] > maxMag) {
                maxMag = suppressedMatrix.map[r][c];
            }
        }
    }

    double lowThreshold = 0.075 * maxMag;
    double highThreshold = 0.15 * maxMag;

    // classify pixels
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
                edges.map[r][c] = 25;//? placeholder value 
            }
        }
    }

    //weak become strong edges
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

    // Remove remaining edge candidates (not connected to strong edges)
    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            if (edges.map[r][c] == 25) {
                edges.map[r][c] = 0;
            }
        }
    }

    Image result = matrix2Image(edges, 0, 1.0);
    return result;
}


int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s input.pgm output.pgm\n", argv[0]);
        return 1;
    }

    Image img = readImage(argv[1]);
    Image result = canny(img);
    writeImage(result, argv[2]);
    deleteImage(img);
    deleteImage(result);
    return 0;
}
