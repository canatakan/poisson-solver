//============================================================================
// CMPE 478 - Homework 1
// Name: Can Atakan Ugur
// Student ID: 2017400057
//============================================================================

#include <iostream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#define N 80
#define TOLERANCE 0.00001

float exact(int x, int y, int z) { // exact solution is xyz
	return x * y * z;
}

float f(int x, int y, int z) { // after taking the second derivatives f(x, y, z) = 0 + 0 + 0
	return 0;
}

int main(void) {

	double u[2][N][N][N]; // at 0: previous, at 1: current

	// initialize the array
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			for (int k = 0; k < N; ++k) {
				if (i == 0 || j == 0 || k == 0 || i == N - 1 || j == N - 1
						|| k == N - 1) { // call exact() on the boundaries
					u[0][i][j][k] = exact(i, j, k);
					u[1][i][j][k] = exact(i, j, k);
				} else { // fill with 0 everywhere else
					u[0][i][j][k] = 0;
					u[1][i][j][k] = 0;
				}
			}
		}
	}

	int t = 1;
	double maxError = 0.0;

	double itime, ftime, exec_time;
	itime = omp_get_wtime();

	do {

		maxError = 0.0; // reset the maximum error

		double difference = 0.0; // hold the maximum of all the differences in this variable

#pragma omp parallel for reduction(max:difference) // maximum difference will be reduced
		for (int i = 1; i < N - 1; ++i) { // skip the boundaries
			for (int j = 1; j < N - 1; ++j) {
				for (int k = 1; k < N - 1; ++k) {

					double neighborSum = u[(t - 1) % 2][i + 1][j][k] // compute the neighbor sum from the previous array
					+ u[(t - 1) % 2][i - 1][j][k] + u[(t - 1) % 2][i][j + 1][k]
							+ u[(t - 1) % 2][i][j - 1][k]
							+ u[(t - 1) % 2][i][j][k + 1]
							+ u[(t - 1) % 2][i][j][k - 1];

					u[t % 2][i][j][k] = (neighborSum / 6) // assign the value computed to the current array
					- (1 / (6 * N * N)) * f(i, j, k);

					difference = fabs(
							u[t % 2][i][j][k] - u[(t - 1) % 2][i][j][k]); // calculate the difference

					if (difference > maxError)
						maxError = difference;

				}
			}
		}
		t = t + 1;
	} while (maxError >= TOLERANCE);

	ftime = omp_get_wtime();
	exec_time = ftime - itime;

	cout << "For N = " << N << ", tolerance = " << TOLERANCE << endl;
	cout << "-------------------------------" << endl;
	cout << "Iteration count (t): " << t << endl;
	cout << "Maximum error: " << maxError << endl;
	cout << "Execution time: " << exec_time << endl;

}
