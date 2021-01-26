#ifndef VERIFY_H
#define VERIFY_H

#include "readMatrix.h"
#include <stdlib.h>

// Verify if L*x = b. In case of descrepancy, compare with the results generated from naive lsolve
void verify(int* &col, int* &row, double* &val, int n, double* &x, std::string rhsFName) {
    double *res = new double[n]();

    // L*x sparse-matrix vector multiplcation
    for (int i = 0; i < n; ++i) {
        for (int j = col[i]; j < col[i+1]; ++j) {
            res[row[j]] += val[j] * x[i];
        }
    }

    double *x1, *x2;
    if (!readRHSMatrix(rhsFName, x1, n)) {
        std::cout << "Error while reading RHS matrix." << std::endl;
        return;
    }

    x2 = new double[n]();
    for (int i = 0; i < n; ++i) {
        x2[i] = x1[i];
    }

    // Naive solve on x2. Used to verify in-case of descrepancy.
    for (int j = 0; j < n; ++j) {
        x2[j] /= val[col[j]];
        int p = col[j] + 1;
        for (; p < col[j+1]; ++p) {
            x2[row[p]] -= val[p] * x2[j];
        }
    }

    int count = 0;
    for (int i = 0; i < n; ++i) {
        double diff = abs(x1[i] - res[i]);

        // Comparing the non-zeros
        if (x1[i] != 0) {
            if (diff > 0.001) {
                // Non-zeros dont match. Does the result match with naive solve?
                diff = abs(x[i] - x2[i]);
                if (diff > 0.001) {
                    std::cout << "Verification failed!" << i << " " << x[i] << " " << x2[i] << std::endl;
                }
            }
        } else {
            if (diff > 0.001) {
                ++count;
            }
        }
    }
    std::cout << "Result verified with " << count << " values of matrix product L*x not matching the zeros in RHS." << std::endl;
}

#endif
