#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "netpbm.h"

#define EDGE_THRESHOLD 40.0
#define MIN_RADIUS 14
#define MAX_RADIUS 65
#define RADIUS_STEP 2
#define ANGLE_STEP_DEG 6
#define NUMBER_OF_CIRCLES 10
#define MIN_CENTER_SEPARATION 22.0

typedef struct
{
	int height, width, depth;
	double ***map;
} HoughSpace;


// Create a 3D Hough space of size height x width x depth and initialize to zero.
HoughSpace createHoughSpace(int height, int width, int depth)
{
	int m, n, r;
	HoughSpace hs;

	hs.height = height;
	hs.width = width;
	hs.depth = depth;
	hs.map = (double ***) malloc(sizeof(double **) * height);

	for (m = 0; m < height; m++)
	{
		hs.map[m] = (double **) malloc(sizeof(double *) * width);
		for (n = 0; n < width; n++)
		{
			hs.map[m][n] = (double *) malloc(sizeof(double) * depth);
			for (r = 0; r < depth; r++)
				hs.map[m][n][r] = 0.0;
		}
	}

	return hs;
}


// Delete a previously created 3D Hough space.
void deleteHoughSpace(HoughSpace hs)
{
	int m, n;

	for (m = 0; m < hs.height; m++)
	{
		for (n = 0; n < hs.width; n++)
			free(hs.map[m][n]);
		free(hs.map[m]);
	}
	free(hs.map);
}


// Project the 3D Hough space to a 2D matrix by taking the maximum vote over all radii
// for each candidate circle center.
Matrix projectHoughSpace(HoughSpace hs)
{
	int m, n, r;
	double best;
	Matrix mx = createMatrix(hs.height, hs.width);

	for (m = 0; m < hs.height; m++)
		for (n = 0; n < hs.width; n++)
		{
			best = 0.0;
			for (r = 0; r < hs.depth; r++)
				if (hs.map[m][n][r] > best)
					best = hs.map[m][n][r];
			mx.map[m][n] = best;
		}

	return mx;
}


// Simple Sobel-based edge detector.
// Produces a binary edge image with values 255.0 for edges and 0.0 otherwise.
Matrix detectEdges(Matrix inputMatrix)
{
	int m, n;
	double gx, gy, mag;
	Matrix edgeMatrix = createMatrix(inputMatrix.height, inputMatrix.width);

	for (m = 1; m < inputMatrix.height - 1; m++)
		for (n = 1; n < inputMatrix.width - 1; n++)
		{
			gx = -inputMatrix.map[m - 1][n - 1] + inputMatrix.map[m - 1][n + 1]
			   - 2.0 * inputMatrix.map[m][n - 1] + 2.0 * inputMatrix.map[m][n + 1]
			   - inputMatrix.map[m + 1][n - 1] + inputMatrix.map[m + 1][n + 1];

			gy = -inputMatrix.map[m - 1][n - 1] - 2.0 * inputMatrix.map[m - 1][n] - inputMatrix.map[m - 1][n + 1]
			   + inputMatrix.map[m + 1][n - 1] + 2.0 * inputMatrix.map[m + 1][n] + inputMatrix.map[m + 1][n + 1];

			mag = sqrt(gx * gx + gy * gy);
			edgeMatrix.map[m][n] = (mag >= EDGE_THRESHOLD) ? 255.0 : 0.0;
		}

	return edgeMatrix;
}


// builds a 3D Hough parameter space for circles.
// Dimensions are center row, center column, and radius index.
HoughSpace houghTransformLines(Matrix mxSpatial, int minRadius, int maxRadius, int radiusStep)
{
	int m, n, a, radius, cy, cx, radiusIndex;
	int depth = (maxRadius - minRadius) / radiusStep + 1;
	int angleCount = 360 / ANGLE_STEP_DEG;
	double alpha;
	double *sinTable;
	double *cosTable;
	HoughSpace hs = createHoughSpace(mxSpatial.height, mxSpatial.width, depth);

	sinTable = (double *) malloc(sizeof(double) * angleCount);
	cosTable = (double *) malloc(sizeof(double) * angleCount);

	for (a = 0; a < angleCount; a++)
	{
		alpha = 2.0 * PI * (double) (a * ANGLE_STEP_DEG) / 360.0;
		sinTable[a] = sin(alpha);
		cosTable[a] = cos(alpha);
	}

	for (m = 0; m < mxSpatial.height; m++)
		for (n = 0; n < mxSpatial.width; n++)
			if (mxSpatial.map[m][n] > 200.0)
				for (radiusIndex = 0; radiusIndex < depth; radiusIndex++)
				{
					radius = minRadius + radiusIndex * radiusStep;
					for (a = 0; a < angleCount; a++)
					{
						cy = (int) (m - radius * sinTable[a] + 0.5);
						cx = (int) (n - radius * cosTable[a] + 0.5);

						if (cy >= 0 && cy < hs.height && cx >= 0 && cx < hs.width)
							hs.map[cy][cx][radiusIndex] += 1.0;
					}
				}

	free(sinTable);
	free(cosTable);

	return hs;
}


// Test whether entry (m, n) in matrix mx is a local maximum, i.e., is not exceeded by any
// of its maximally 8 neighbors. Return 1 if true, 0 otherwise.
int isLocalMaximum(Matrix mx, int m, int n)
{
	double strength = mx.map[m][n];
	int i, j;
	int iMin = (m == 0) ? 0 : (m - 1);
	int iMax = (m == mx.height - 1) ? m : (m + 1);
	int jMin = (n == 0) ? 0 : (n - 1);
	int jMax = (n == mx.width - 1) ? n : (n + 1);

	for (i = iMin; i <= iMax; i++)
		for (j = jMin; j <= jMax; j++)
			if (mx.map[i][j] > strength)
				return 0;

	return 1;
}


// Insert a new entry, consisting of vPos, hPos, and strength, into the list of maxima mx.
void insertMaxEntry(Matrix mx, int vPos, int hPos, double strength)
{
	int m, n = mx.width - 1;

	while (n > 0 && mx.map[2][n - 1] < strength)
	{
		for (m = 0; m < 3; m++)
			mx.map[m][n] = mx.map[m][n - 1];
		n--;
	}

	mx.map[0][n] = (double) vPos;
	mx.map[1][n] = (double) hPos;
	mx.map[2][n] = strength;
}


// Delete entry number i from the list of maxima mx.
void deleteMaxEntry(Matrix mx, int i)
{
	int m, n;

	for (n = i; n < mx.width - 1; n++)
		for (m = 0; m < 3; m++)
			mx.map[m][n] = mx.map[m][n + 1];

	mx.map[2][mx.width - 1] = -1.0;
}


// Find the <number> highest maxima in a 2D Hough projection that are separated by a
// Euclidean distance of at least <minSeparation>. The result is a three-row matrix
// with each column representing the row, the column, and the strength of one maximum,
// in descending order of strength.
Matrix findHoughMaxima(Matrix mx, int number, double minSeparation)
{
	int j, m, n, k;
	double minSepSquare = SQR(minSeparation);
	double strength;
	Matrix maxima = createMatrix(3, number);

	for (j = 0; j < number; j++)
		maxima.map[2][j] = -1.0;

	for (m = 0; m < mx.height; m++)
		for (n = 0; n < mx.width; n++)
			if (mx.map[m][n] > 0.0 && isLocalMaximum(mx, m, n))
			{
				strength = mx.map[m][n];

				if (strength > maxima.map[2][number - 1])
				{
					insertMaxEntry(maxima, m, n, strength);

					for (j = 0; j < number; j++)
					{
						if (maxima.map[2][j] < 0.0)
							continue;

						for (k = j + 1; k < number; )
						{
							if (maxima.map[2][k] < 0.0)
							{
								k++;
								continue;
							}

							if (SQR(maxima.map[0][j] - maxima.map[0][k]) +
								SQR(maxima.map[1][j] - maxima.map[1][k]) < minSepSquare)
								deleteMaxEntry(maxima, k);
							else
								k++;
						}
					}
				}
			}

	return maxima;
}


// For a selected center (m, n), find a good radius.
// This prefers the largest radius whose vote is at least 85% of the best vote,
// which helps favor the outer boundary over smaller inner circular details.
int bestRadiusAtCenter(HoughSpace hs, int m, int n, int minRadius, int radiusStep)
{
	int r, bestIndex = 0, chosenIndex = 0;
	double bestVote = hs.map[m][n][0];

	for (r = 1; r < hs.depth; r++)
		if (hs.map[m][n][r] > bestVote)
		{
			bestVote = hs.map[m][n][r];
			bestIndex = r;
		}

	chosenIndex = bestIndex;
	for (r = bestIndex; r < hs.depth; r++)
		if (hs.map[m][n][r] >= 0.85 * bestVote)
			chosenIndex = r;

	return minRadius + chosenIndex * radiusStep;
}


// Read image and write Hough transform related output images.

int main(int argc, char *argv[])
{
	int i;
	int radius;
	double minStrength;
	Image inputImage;
	Matrix inputMatrix;
	Matrix edgeMatrix, houghMatrix, maxMatrix;
	HoughSpace houghSpace;
	Image edgeImage, houghImage;

	if (argc != 5)
	{
		printf("Usage: %s input.pgm edge_output.pgm hough_output.pgm circles_output.ppm\n", argv[0]);
		return 1;
	}

	inputImage = readImage(argv[1]);
	inputMatrix = image2Matrix(inputImage);

	edgeMatrix = detectEdges(inputMatrix);
	edgeImage = matrix2Image(edgeMatrix, 1, 1.0);
	writeImage(edgeImage, argv[2]);

	houghSpace = houghTransformLines(edgeMatrix, MIN_RADIUS, MAX_RADIUS, RADIUS_STEP);
	houghMatrix = projectHoughSpace(houghSpace);
	houghImage = matrix2Image(houghMatrix, 1, 1.0);
	writeImage(houghImage, argv[3]);

	maxMatrix = findHoughMaxima(houghMatrix, NUMBER_OF_CIRCLES, MIN_CENTER_SEPARATION);

	for (i = 0; i < NUMBER_OF_CIRCLES; i++)
		if (maxMatrix.map[2][i] >= 0.0)
			ellipse(houghImage,
				(int) maxMatrix.map[0][i],
				(int) maxMatrix.map[1][i],
				8, 8,
				2, 10, 7,
				255, 255, 255, 0);

	minStrength = 0.3 * maxMatrix.map[2][0];

	for (i = 0; i < NUMBER_OF_CIRCLES; i++)
		if (maxMatrix.map[2][i] >= 0.0)
		{
			if (maxMatrix.map[2][i] < minStrength)
				continue;

			if (maxMatrix.map[0][i] < MAX_RADIUS || maxMatrix.map[0][i] >= inputImage.height - MAX_RADIUS ||
				maxMatrix.map[1][i] < MAX_RADIUS || maxMatrix.map[1][i] >= inputImage.width - MAX_RADIUS)
				continue;

			radius = bestRadiusAtCenter(houghSpace,
				(int) maxMatrix.map[0][i],
				(int) maxMatrix.map[1][i],
				MIN_RADIUS, RADIUS_STEP);

			if (radius < 24)
				continue;

			ellipse(inputImage,
				(int) maxMatrix.map[0][i],
				(int) maxMatrix.map[1][i],
				radius, radius,
				2, 18, 10,
				255, 0, 0, 0);
		}

	writeImage(inputImage, argv[4]);

	deleteMatrix(inputMatrix);
	deleteMatrix(edgeMatrix);
	deleteHoughSpace(houghSpace);
	deleteMatrix(houghMatrix);
	deleteMatrix(maxMatrix);
	deleteImage(inputImage);
	deleteImage(edgeImage);
	deleteImage(houghImage);

	return 0;
}
