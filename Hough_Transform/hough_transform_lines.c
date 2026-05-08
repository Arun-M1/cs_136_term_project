#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "netpbm.h"

#define EDGE_THRESHOLD 40.0
#define ALPHA_STEP_DEG 1
#define NUMBER_OF_LINES 10
#define MIN_LINE_SEPARATION 20.0

typedef struct
{
	int height, width;
	double **map;
} HoughSpace;


HoughSpace createHoughSpace(int height, int width)
{
	int m, n;
	HoughSpace hs;

	hs.height = height;
	hs.width = width;
	hs.map = (double **) malloc(sizeof(double *) * height);

	for (m = 0; m < height; m++)
	{
		hs.map[m] = (double *) malloc(sizeof(double) * width);
		for (n = 0; n < width; n++)
			hs.map[m][n] = 0.0;
	}

	return hs;
}


void deleteHoughSpace(HoughSpace hs)
{
	int m;

	for (m = 0; m < hs.height; m++)
		free(hs.map[m]);
	free(hs.map);
}


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


HoughSpace houghTransformLines(Matrix mxSpatial, int alphaStepDeg)
{
	int m, n, a;
	int numAlpha = 360 / alphaStepDeg;
	int maxD = (int)(sqrt(mxSpatial.height * mxSpatial.height + 
	                      mxSpatial.width * mxSpatial.width) + 0.5);
	double alpha, d;
	double *sinTable, *cosTable;
	HoughSpace hs = createHoughSpace(numAlpha, maxD);

	sinTable = (double *) malloc(sizeof(double) * numAlpha);
	cosTable = (double *) malloc(sizeof(double) * numAlpha);

	for (a = 0; a < numAlpha; a++)
	{
		alpha = 2.0 * PI * (double)(a * alphaStepDeg) / 360.0;
		sinTable[a] = sin(alpha);
		cosTable[a] = cos(alpha);
	}

	for (m = 0; m < mxSpatial.height; m++)
	{
		for (n = 0; n < mxSpatial.width; n++)
		{
			if (mxSpatial.map[m][n] > 200.0)
			{
				for (a = 0; a < numAlpha; a++)
				{
					d = m * cosTable[a] + n * sinTable[a];

					int dIndex = (int)(d + 0.5);
					if (dIndex >= 0 && dIndex < maxD)
					{
						hs.map[a][dIndex] += 1.0;
					}
				}
			}
		}
	}

	free(sinTable);
	free(cosTable);

	return hs;
}


Matrix houghSpace2Matrix(HoughSpace hs)
{
	int m, n;
	Matrix mx = createMatrix(hs.height, hs.width);

	for (m = 0; m < hs.height; m++)
		for (n = 0; n < hs.width; n++)
			mx.map[m][n] = hs.map[m][n];

	return mx;
}


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


void deleteMaxEntry(Matrix mx, int i)
{
	int m, n;

	for (n = i; n < mx.width - 1; n++)
		for (m = 0; m < 3; m++)
			mx.map[m][n] = mx.map[m][n + 1];

	mx.map[2][mx.width - 1] = -1.0;
}


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


void drawLine(Image img, int alphaIndex, int dValue, int alphaStepDeg)
{
	double alpha = 2.0 * PI * (double)(alphaIndex * alphaStepDeg) / 360.0;
	double cosAlpha = cos(alpha);
	double sinAlpha = sin(alpha);
	int p, i, j;
	int lastI = -1000, lastJ = -1000;
	
	for (p = -2000; p < 2000; p++)
	{
		i = (int)(dValue * cosAlpha + p * sinAlpha + 0.5);
		j = (int)(dValue * sinAlpha - p * cosAlpha + 0.5);
		
		if (i >= 0 && i < img.height && j >= 0 && j < img.width)
		{
			if (i != lastI || j != lastJ)
			{
				setPixel(img, i, j, 255, 0, 0, NO_CHANGE);
				lastI = i;
				lastJ = j;
			}
		}
	}
}


int main(int argc, char *argv[])
{
	int i;
	char outputFile[512];
	
	if (argc != 3)
	{
		printf("Usage: %s input.pgm output_prefix\n", argv[0]);
		printf("Example: %s coastline.pgm results/coastline\n", argv[0]);
		printf("Output files: output_prefix_edges.pgm, output_prefix_hough.pgm, output_prefix_lines.ppm\n");
		return 1;
	}
	
	Image inputImage = readImage(argv[1]);
	Matrix inputMatrix = image2Matrix(inputImage);
	Matrix edgeMatrix, houghMatrix, maxMatrix;
	HoughSpace houghSpace;
	Image edgeImage, houghImage, resultImage;

	// Edge detection
	edgeMatrix = detectEdges(inputMatrix);
	edgeImage = matrix2Image(edgeMatrix, 0, 1.0);
	sprintf(outputFile, "%s_edges.pgm", argv[2]);
	writeImage(edgeImage, outputFile);

	// Hough transform
	houghSpace = houghTransformLines(edgeMatrix, ALPHA_STEP_DEG);
	houghMatrix = houghSpace2Matrix(houghSpace);
	houghImage = matrix2Image(houghMatrix, 1, 1.0);
	sprintf(outputFile, "%s_hough.pgm", argv[2]);
	writeImage(houghImage, outputFile);

	// Find maxima
	maxMatrix = findHoughMaxima(houghMatrix, NUMBER_OF_LINES, MIN_LINE_SEPARATION);

	// Mark maxima in Hough space
	for (i = 0; i < NUMBER_OF_LINES; i++)
	{
		if (maxMatrix.map[2][i] >= 0.0)
		{
			int alphaIdx = (int)maxMatrix.map[0][i];
			int dIdx = (int)maxMatrix.map[1][i];
			filledEllipse(houghImage, alphaIdx, dIdx, 5, 5, 255, 255, 0, NO_CHANGE);
		}
	}
	sprintf(outputFile, "%s_hough_max.ppm", argv[2]);
	writeImage(houghImage, outputFile);

	// Draw detected lines
	resultImage = readImage(argv[1]);
	double minStrength = 0.3 * maxMatrix.map[2][0];

	for (i = 0; i < NUMBER_OF_LINES; i++)
	{
		if (maxMatrix.map[2][i] >= minStrength)
		{
			int alphaIdx = (int)maxMatrix.map[0][i];
			int dIdx = (int)maxMatrix.map[1][i];
			drawLine(resultImage, alphaIdx, dIdx, ALPHA_STEP_DEG);
		}
	}

	sprintf(outputFile, "%s_lines.ppm", argv[2]);
	writeImage(resultImage, outputFile);

	// Cleanup
	deleteMatrix(inputMatrix);
	deleteMatrix(edgeMatrix);
	deleteHoughSpace(houghSpace);
	deleteMatrix(houghMatrix);
	deleteMatrix(maxMatrix);
	deleteImage(inputImage);
	deleteImage(edgeImage);
	deleteImage(houghImage);
	deleteImage(resultImage);

	return 0;
}