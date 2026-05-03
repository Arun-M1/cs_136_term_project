// fom.c — Pratt's Figure of Merit evaluation
// Usage: ./fom predicted.pgm groundtruth.pgm
//
// Uses the Felzenszwalb-Huttenlocher separable Distance Transform (O(W*H))
// to compute exact squared Euclidean distances from each pixel to the
// nearest ground-truth edge pixel.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "netpbm.h"

#define ALPHA (1.0 / 9.0)

// --- Felzenszwalb-Huttenlocher 1-D distance transform helpers ---

// Compute 1D squared-distance transform of a row/column of length n.
// f[i] holds the "source" value: 0 if it is a GT edge pixel, INF otherwise.
// dt[i] = min_k { (i-k)^2 + f[k] }
static void dt1d(const double *f, double *dt, int *v, double *z, int n) {
    int k = 0;
    v[0] = 0;
    z[0] = -1e38;
    z[1] =  1e38;

    for (int q = 1; q < n; q++) {
        double s;
        while (1) {
            // intersection of parabolas at v[k] and q
            s = ((f[q] + q * q) - (f[v[k]] + v[k] * v[k])) / (2.0 * q - 2.0 * v[k]);
            if (s > z[k])
                break;
            k--;
        }
        k++;
        v[k]     = q;
        z[k]     = s;
        z[k + 1] = 1e38;
    }

    k = 0;
    for (int q = 0; q < n; q++) {
        while (z[k + 1] < q)
            k++;
        int d = q - v[k];
        dt[q] = d * d + f[v[k]];
    }
}

// Compute 2-D squared Euclidean DT of a binary map (1 = GT edge, 0 = background).
// Returns a height×width array of squared distances (allocated, caller frees).
static double **compute_dt(int **gt, int height, int width) {
    const double INF = 1e30;
    int maxdim = (height > width) ? height : width;

    double  *f  = malloc(maxdim * sizeof(double));
    double  *dt = malloc(maxdim * sizeof(double));
    int     *v  = malloc(maxdim * sizeof(int));
    double  *z  = malloc((maxdim + 1) * sizeof(double));

    // Allocate 2-D result
    double **D = malloc(height * sizeof(double *));
    for (int r = 0; r < height; r++)
        D[r] = malloc(width * sizeof(double));

    // Pass 1: per-column DT
    for (int c = 0; c < width; c++) {
        for (int r = 0; r < height; r++)
            f[r] = gt[r][c] ? 0.0 : INF;
        dt1d(f, dt, v, z, height);
        for (int r = 0; r < height; r++)
            D[r][c] = dt[r];
    }

    // Pass 2: per-row DT using column results
    for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++)
            f[c] = D[r][c];
        dt1d(f, dt, v, z, width);
        for (int c = 0; c < width; c++)
            D[r][c] = dt[c];
    }

    free(f); free(dt); free(v); free(z);
    return D;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s predicted.pgm groundtruth.pgm\n", argv[0]);
        return 1;
    }

    Image pred = readImage(argv[1]);
    Image gt   = readImage(argv[2]);

    if (pred.height != gt.height || pred.width != gt.width) {
        fprintf(stderr, "Error: image dimensions do not match (%dx%d vs %dx%d)\n",
                pred.width, pred.height, gt.width, gt.height);
        return 1;
    }

    int H = pred.height, W = pred.width;

    // Build binary int maps
    int **pred_bin = malloc(H * sizeof(int *));
    int **gt_bin   = malloc(H * sizeof(int *));
    for (int r = 0; r < H; r++) {
        pred_bin[r] = malloc(W * sizeof(int));
        gt_bin[r]   = malloc(W * sizeof(int));
        for (int c = 0; c < W; c++) {
            pred_bin[r][c] = (pred.map[r][c].i > 0) ? 1 : 0;
            gt_bin[r][c]   = (gt.map[r][c].i   > 0) ? 1 : 0;
        }
    }

    // Count N_A and N_I
    long N_A = 0, N_I = 0;
    for (int r = 0; r < H; r++)
        for (int c = 0; c < W; c++) {
            if (pred_bin[r][c]) N_A++;
            if (gt_bin[r][c])   N_I++;
        }

    // Compute squared DT of GT map
    double **D = compute_dt(gt_bin, H, W);

    // Sum Pratt's FOM numerator
    double sum = 0.0;
    for (int r = 0; r < H; r++)
        for (int c = 0; c < W; c++)
            if (pred_bin[r][c])
                sum += 1.0 / (1.0 + ALPHA * D[r][c]);

    long denom = (N_A > N_I) ? N_A : N_I;
    double fom = (denom > 0) ? sum / denom : 0.0;

    printf("FOM: %.6f  (N_detected=%ld, N_gt=%ld)\n", fom, N_A, N_I);

    // Cleanup
    for (int r = 0; r < H; r++) {
        free(pred_bin[r]);
        free(gt_bin[r]);
        free(D[r]);
    }
    free(pred_bin); free(gt_bin); free(D);
    deleteImage(pred); deleteImage(gt);

    return 0;
}
