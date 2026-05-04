#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../netpbm.h"

#define LAWS_SIZE 5
#define LAWS_COUNT 5

void buildLawsVectors(int laws[LAWS_COUNT][LAWS_SIZE]);
void buildLawsMasks(int laws[LAWS_COUNT][LAWS_SIZE], int masks[25][5][5]);
Matrix convolveLawsMask(Matrix input, int mask[5][5]);
Image normalizeMatrixToImage(Matrix m);
Matrix absMatrix(Matrix input);
Matrix localAverage(Matrix input, int windowSize);
double patchMeanFromMatrix(Matrix m, int r0, int c0, int patchSize);
void printFeatureVectors(double **features, int patchRows, int patchCols, int featureDim);
void normalizeFeatures(double **features, int count, int dim);
double squaredDistance(double *a, double *b, int dim);
void copyVector(double *dst, double *src, int dim);
void kmeans(double **features, int count, int dim, int k, int *labels);
Image labelsToColorImage(int *labels, int patchRows, int patchCols, int patchSize, int height, int width);
void printLabelsGrid(int *labels, int patchRows, int patchCols);
Image segmentTexture(Image inputImg, int segments);

void buildLawsVectors(int laws[LAWS_COUNT][LAWS_SIZE]) {
    int L5[5] = {1, 4, 6, 4, 1};
    int E5[5] = {-1, -2, 0, 2, 1};
    int S5[5] = {-1, 0, 2, 0, -1};
    int R5[5] = {1, -4, 6, -4, 1};
    int W5[5] = {-1, 2, 0, -2, 1};

    for (int i = 0; i < 5; i++) {
        laws[0][i] = L5[i];
        laws[1][i] = E5[i];
        laws[2][i] = S5[i];
        laws[3][i] = R5[i];
        laws[4][i] = W5[i];
    }
}

void buildLawsMasks(int laws[LAWS_COUNT][LAWS_SIZE], int masks[25][5][5]) {
    int idx = 0;

    for (int a = 0; a < LAWS_COUNT; a++) {
        for (int b = 0; b < LAWS_COUNT; b++) {
            for (int i = 0; i < LAWS_SIZE; i++) {
                for (int j = 0; j < LAWS_SIZE; j++) {
                    masks[idx][i][j] = laws[a][i] * laws[b][j];
                }
            }
            idx++;
        }
    }
}

Matrix convolveLawsMask(Matrix input, int mask[5][5]) {
    Matrix output = createMatrix(input.height, input.width);

    for (int r = 0; r < input.height; r++) {
        for (int c = 0; c < input.width; c++) {
            double sum = 0.0;

            for (int i = -2; i <= 2; i++) {
                for (int j = -2; j <= 2; j++) {
                    int rr = r + i;
                    int cc = c + j;

                    if (rr >= 0 && rr < input.height && cc >= 0 && cc < input.width) {
                        sum += input.map[rr][cc] * mask[i + 2][j + 2];
                    }
                }
            }

            output.map[r][c] = sum;
        }
    }

    return output;
}

Image normalizeMatrixToImage(Matrix m) {
    Image out = createImage(m.height, m.width);

    double minVal = m.map[0][0];
    double maxVal = m.map[0][0];

    for (int r = 0; r < m.height; r++) {
        for (int c = 0; c < m.width; c++) {
            if (m.map[r][c] < minVal) minVal = m.map[r][c];
            if (m.map[r][c] > maxVal) maxVal = m.map[r][c];
        }
    }

    double range = maxVal - minVal;
    if (range == 0.0) range = 1.0;

    for (int r = 0; r < m.height; r++) {
        for (int c = 0; c < m.width; c++) {
            int value = (int)(255.0 * (m.map[r][c] - minVal) / range);

            if (value < 0) value = 0;
            if (value > 255) value = 255;

            out.map[r][c].i = value;
            out.map[r][c].r = value;
            out.map[r][c].g = value;
            out.map[r][c].b = value;
        }
    }

    return out;
}

Matrix absMatrix(Matrix input) {
    Matrix output = createMatrix(input.height, input.width);

    for (int r = 0; r < input.height; r++) {
        for (int c = 0; c < input.width; c++) {
            output.map[r][c] = fabs(input.map[r][c]);
        }
    }

    return output;
}

Matrix localAverage(Matrix input, int windowSize) {
    Matrix output = createMatrix(input.height, input.width);
    int half = windowSize / 2;

    for (int r = 0; r < input.height; r++) {
        for (int c = 0; c < input.width; c++) {
            double sum = 0.0;
            int count = 0;

            for (int dr = -half; dr <= half; dr++) {
                for (int dc = -half; dc <= half; dc++) {
                    int rr = r + dr;
                    int cc = c + dc;

                    if (rr >= 0 && rr < input.height && cc >= 0 && cc < input.width) {
                        sum += input.map[rr][cc];
                        count++;
                    }
                }
            }

            if (count > 0) {
                output.map[r][c] = sum / count;
            } else {
                output.map[r][c] = 0.0;
            }
        }
    }

    return output;
}

double patchMeanFromMatrix(Matrix m, int r0, int c0, int patchSize) {
    double sum = 0.0;
    int count = 0;

    for (int r = r0; r < r0 + patchSize && r < m.height; r++) {
        for (int c = c0; c < c0 + patchSize && c < m.width; c++) {
            sum += m.map[r][c];
            count++;
        }
    }

    if (count == 0) return 0.0;
    return sum / count;
}

void printFeatureVectors(double **features, int patchRows, int patchCols, int featureDim) {
    int total = patchRows * patchCols;

    printf("\nFeature vectors (showing first few patches):\n");
    for (int i = 0; i < total && i < 10; i++) {
        printf("Patch %d: [", i);
        for (int d = 0; d < featureDim; d++) {
            printf("%.3f", features[i][d]);
            if (d < featureDim - 1) printf(", ");
        }
        printf("]\n");
    }
}

void normalizeFeatures(double **features, int count, int dim) {
    for (int d = 0; d < dim; d++) {
        double minVal = features[0][d];
        double maxVal = features[0][d];

        for (int i = 1; i < count; i++) {
            if (features[i][d] < minVal) minVal = features[i][d];
            if (features[i][d] > maxVal) maxVal = features[i][d];
        }

        double range = maxVal - minVal;
        if (range == 0.0) range = 1.0;

        for (int i = 0; i < count; i++) {
            features[i][d] = (features[i][d] - minVal) / range;
        }
    }
}

double squaredDistance(double *a, double *b, int dim) {
    double sum = 0.0;

    for (int i = 0; i < dim; i++) {
        double diff = a[i] - b[i];
        sum += diff * diff;
    }

    return sum;
}

void copyVector(double *dst, double *src, int dim) {
    for (int i = 0; i < dim; i++) {
        dst[i] = src[i];
    }
}

void kmeans(double **features, int count, int dim, int k, int *labels) {
    double **centroids = (double **)malloc(k * sizeof(double *));
    double **newCentroids = (double **)malloc(k * sizeof(double *));
    int *clusterSizes = (int *)malloc(k * sizeof(int));

    int *chosen = (int *)malloc(k * sizeof(int));
    for (int c = 0; c < k; c++) {
        centroids[c] = (double *)malloc(dim * sizeof(double));
        newCentroids[c] = (double *)malloc(dim * sizeof(double));
        int pick, duplicate;
        do {
            pick = rand() % count;
            duplicate = 0;
            for (int prev = 0; prev < c; prev++) {
                if (chosen[prev] == pick) { duplicate = 1; break; }
            }
        } while (duplicate);
        chosen[c] = pick;
        copyVector(centroids[c], features[pick], dim);
    }
    free(chosen);

    for (int i = 0; i < count; i++) {
        labels[i] = -1;
    }

    int changed = 1;
    int iterations = 0;
    int maxIterations = 100;

    while (changed && iterations < maxIterations) {
        changed = 0;
        iterations++;

        for (int i = 0; i < count; i++) {
            double bestDist = squaredDistance(features[i], centroids[0], dim);
            int bestCluster = 0;

            for (int c = 1; c < k; c++) {
                double dist = squaredDistance(features[i], centroids[c], dim);
                if (dist < bestDist) {
                    bestDist = dist;
                    bestCluster = c;
                }
            }

            if (labels[i] != bestCluster) {
                labels[i] = bestCluster;
                changed = 1;
            }
        }

        for (int c = 0; c < k; c++) {
            clusterSizes[c] = 0;
            for (int d = 0; d < dim; d++) {
                newCentroids[c][d] = 0.0;
            }
        }

        for (int i = 0; i < count; i++) {
            int cluster = labels[i];
            clusterSizes[cluster]++;

            for (int d = 0; d < dim; d++) {
                newCentroids[cluster][d] += features[i][d];
            }
        }

        for (int c = 0; c < k; c++) {
            if (clusterSizes[c] > 0) {
                for (int d = 0; d < dim; d++) {
                    centroids[c][d] = newCentroids[c][d] / clusterSizes[c];
                }
            }
        }
    }

    for (int c = 0; c < k; c++) {
        free(centroids[c]);
        free(newCentroids[c]);
    }
    free(centroids);
    free(newCentroids);
    free(clusterSizes);
}

Image labelsToColorImage(int *labels, int patchRows, int patchCols, int patchSize, int height, int width) {
    Image result = createImage(height, width);

    Pixel palette[10];
    palette[0].r = 255; palette[0].g = 0;   palette[0].b = 0;   palette[0].i = 255;
    palette[1].r = 0;   palette[1].g = 255; palette[1].b = 0;   palette[1].i = 255;
    palette[2].r = 0;   palette[2].g = 0;   palette[2].b = 255; palette[2].i = 255;
    palette[3].r = 255; palette[3].g = 255; palette[3].b = 0;   palette[3].i = 255;
    palette[4].r = 255; palette[4].g = 0;   palette[4].b = 255; palette[4].i = 255;
    palette[5].r = 0;   palette[5].g = 255; palette[5].b = 255; palette[5].i = 255;
    palette[6].r = 255; palette[6].g = 128; palette[6].b = 0;   palette[6].i = 255;
    palette[7].r = 128; palette[7].g = 0;   palette[7].b = 255; palette[7].i = 255;
    palette[8].r = 0;   palette[8].g = 128; palette[8].b = 255; palette[8].i = 255;
    palette[9].r = 180; palette[9].g = 180; palette[9].b = 180; palette[9].i = 255;

    int idx = 0;
    for (int pr = 0; pr < patchRows; pr++) {
        for (int pc = 0; pc < patchCols; pc++) {
            int label = labels[idx] % 10;
            Pixel color = palette[label];

            int r0 = pr * patchSize;
            int c0 = pc * patchSize;
            int r1 = r0 + patchSize;
            int c1 = c0 + patchSize;

            if (r1 > height) r1 = height;
            if (c1 > width) c1 = width;

            for (int r = r0; r < r1; r++) {
                for (int c = c0; c < c1; c++) {
                    result.map[r][c].r = color.r;
                    result.map[r][c].g = color.g;
                    result.map[r][c].b = color.b;
                    result.map[r][c].i = color.i;
                }
            }

            idx++;
        }
    }

    return result;
}

void printLabelsGrid(int *labels, int patchRows, int patchCols) {
    printf("\nCluster labels by patch:\n");
    for (int pr = 0; pr < patchRows; pr++) {
        for (int pc = 0; pc < patchCols; pc++) {
            printf("%2d ", labels[pr * patchCols + pc]);
        }
        printf("\n");
    }
}

Image segmentTexture(Image inputImg, int segments) {
    int laws[LAWS_COUNT][LAWS_SIZE];
    int masks[25][LAWS_SIZE][LAWS_SIZE];

    buildLawsVectors(laws);
    buildLawsMasks(laws, masks);

    Matrix gray = image2Matrix(inputImg);

    // 9 Law's masks chosen to capture level, edge, spot, ripple, wave energy
    int selectedMasks[9] = {1, 5, 2, 10, 3, 15, 4, 20, 12};
    int featureDim = 9;
    int patchSize  = 8;

    // --- Step 1: compute texture energy maps ---
    Matrix energyMaps[9];
    for (int k = 0; k < featureDim; k++) {
        Matrix filtered    = convolveLawsMask(gray, masks[selectedMasks[k]]);
        Matrix absFiltered = absMatrix(filtered);
        energyMaps[k]      = localAverage(absFiltered, 9);
        deleteMatrix(filtered);
        deleteMatrix(absFiltered);
    }

    // --- Step 2: build patch feature vectors ---
    int patchRows  = (gray.height + patchSize - 1) / patchSize;
    int patchCols  = (gray.width  + patchSize - 1) / patchSize;
    int numPatches = patchRows * patchCols;

    double **features = (double **)malloc(numPatches * sizeof(double *));
    for (int i = 0; i < numPatches; i++)
        features[i] = (double *)malloc(featureDim * sizeof(double));

    int idx = 0;
    for (int pr = 0; pr < patchRows; pr++) {
        for (int pc = 0; pc < patchCols; pc++) {
            int r0 = pr * patchSize;
            int c0 = pc * patchSize;
            for (int d = 0; d < featureDim; d++)
                features[idx][d] = patchMeanFromMatrix(energyMaps[d], r0, c0, patchSize);
            idx++;
        }
    }

    // --- Step 3: normalize and cluster ---
    normalizeFeatures(features, numPatches, featureDim);

    int *labels = (int *)malloc(numPatches * sizeof(int));
    kmeans(features, numPatches, featureDim, segments, labels);

    // --- Step 4: paint each patch with its cluster color ---
    Image result = labelsToColorImage(labels, patchRows, patchCols, patchSize, inputImg.height, inputImg.width);

    // --- Cleanup ---
    for (int d = 0; d < featureDim; d++)
        deleteMatrix(energyMaps[d]);
    for (int i = 0; i < numPatches; i++)
        free(features[i]);
    free(features);
    free(labels);
    deleteMatrix(gray);

    return result;
}



int main(int argc, char *argv[]) {
    Image img;
    Image result;
    int segments;

    if (argc != 4) {
        printf("Usage: %s input.pgm output.ppm segments\n", argv[0]);
        return 1;
    }

    img = readImage(argv[1]);
    segments = atoi(argv[3]);
    result = segmentTexture(img, segments);
    writeImage(result, argv[2]);
    deleteImage(img);
    deleteImage(result);
    return 0;
}
