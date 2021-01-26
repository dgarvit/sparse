#include <omp.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <chrono>
#include <sstream>
#include <set>
#include <vector>
#include <assert.h>
#include "readMatrix.h"
#include "makeSet.h"

// Verify if L*x = b. In case of descrepancy, compare with the results generated from naive lsolve
int verify(int* &col, int* &row, double* &val, int n, double* &x, std::string rhsFName) {
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
        return -1;
    }

    x2 = new double[n]();
    for (int i = 0; i < n; ++i) {
        x2[i] = x1[i];
    }

    int count = 0;
    for (int i = 0; i < n; ++i) {
        double diff = x1[i] - res[i];
        if (diff < 0) {
            diff *= -1;
        }
        if (x1[i] != 0) {
            if (diff > 0.001) {
                
            }
        }
        // if (diff > 0.001) {
        //     ++count;
        //     std::cout << i << " " << res[i] << " " << x1[i] << std::endl;
        // }
    }
    return count;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "Usage: ./<bin> /path/to/L.mtx /path/to/RHS.mtx" << std::endl;
        return -1;
    }
    std::string fName = argv[1];
    std::string rhsFName = argv[2];
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;
    // int chunk = atoi(argv[1]);
    int *col, *row;
    double  *y, *val, *x, *x1, *x2;

    // printLower(fName);
    // exit(0);

    int n, nnz;
    // std::cout << "here";
    if (!readSparseMatrix(fName,n,nnz,col,row,val)) {
        std::cout << "Error while reading the sparse matrix. Make sure that it is a lower-triangular matrix in matrix market format." << std::endl;
        return -1;
    }
    if (!readRHSMatrix(rhsFName, x, n)) {
        std::cout << "Error while reading RHS matrix." << std::endl;
        return -1;
    }
    x1 = new double[n]();
    for (int i = 0; i < n; ++i) {
        x1[i] = x[i];
    }

    // Making Reachset (First constructing the DAG)
    bool *visited = new bool[n]();
    std::vector<std::set<int>> dag = makeDAG(col, row, val, n, visited);
    std::set<int> reachset = makeReachset(dag, visited, n, x);

#if 0
    // Reach set made. Now Solve
    start = std::chrono::system_clock::now();
    for (auto j : reachset) {
        // std::cout << f + 1 << ",";
        x[j] /= val[col[j]];
        int p = col[j] + 1;
        for (; p < col[j+1]; ++p) {
            x[row[p]] -= val[p] * x[j];
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout << "Optimised lsolve: " << elapsed_seconds.count() << std::endl;
    //*******************************************************************
    // Serial Lsolve
    start = std::chrono::system_clock::now();
    for (int j = 0; j < n; ++j) {
        // std::cout << f + 1 << ",";
        x1[j] /= val[col[j]];
        int p = col[j] + 1;
        for (; p < col[j+1]; ++p) {
            x1[row[p]] -= val[p] * x1[j];
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout << "Naive lsolve: " << elapsed_seconds.count() << std::endl;
    //***************************************************************
    // Verify
    for (int i = 0; i < n; ++i) {
        double diff = x[i] - x1[i];
        if (diff < 0) {
            diff *= -1;
        }
        assert(diff < 0.001);
    }
    //****************************************************************
#endif
#if 1
    std::vector< std::vector<int> > levels = makeLevelSet(dag, reachset, n); // Store levels of parallelism. Tasks inside same level can be performed parallely; each 2 levels need to sync
    
    // Levels created, now parallel solve
    // start = std::chrono::system_clock::now();
    // for (int i = 0; i < levels.size(); ++i) {
    //     // std::cout << "Level " << i << std::endl;
    //     #pragma omp parallel for
    //     for (int j = 0; j < levels[i].size(); ++j) {
    //         // std::cout << omp_get_thread_num() << " performing " << levels[i][j] + 1 << std::endl;
    //         int idx = levels[i][j];
    //         x[idx] /= val[col[idx]];
    //         int p = col[idx] + 1;
    //         for (; p < col[idx+1]; ++p) {
    //             x[row[p]] -= val[p] * x[idx];
    //         }
    //     }
    // }

    int currLevel = 0;
    start = std::chrono::system_clock::now();
    #pragma omp parallel
    {
        while (currLevel < levels.size()) {
            #pragma omp for schedule(dynamic, 5)
            for (int j = 0; j < levels[currLevel].size(); ++j) {
                int idx = levels[currLevel][j];
                x[idx] /= val[col[idx]];
                int p = col[idx] + 1;
                for (; p < col[idx+1]; ++p) {
                    x[row[p]] -= val[p] * x[idx];
                }
            }
            #pragma omp master
            ++currLevel;
            #pragma omp barrier
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    // std::cout << "Parallel lsolve: " << elapsed_seconds.count() << std::endl;
    std::cout << elapsed_seconds.count() << std::endl;

    // for (int j = 0; j < n; ++j) {
    //     x1[j] /= val[col[j]];
    //     int p = col[j] + 1;
    //     for (; p < col[j+1]; ++p) {
    //         x1[row[p]] -= val[p] * x1[j];
    //     }
    // }
    // std::cout << verify(col, row, val, n, x, rhsFName) << std::endl;
    // VERIFY
    // double *res = new double[n]();
    // for (int i = 0; i < n; ++i) {
    //     for (int j = col[i]; j < col[i+1]; ++j) {
    //         res[row[j]] += val[j] * x[i];
    //     }
    // }

    // for (int i = 0; i < n; ++i) {
    //     // std::cout << res[i] << " " << x1[i] << std::endl;
    //     double diff = x1[i] - res[i];
    //     if (diff < 0) {
    //         diff *= -1;
    //     }
    //     if (diff > 0.001) {
    //         std::cout << i << " " << res[i] << " " << x1[i] << std::endl;
    //     }
    // }

    // int count = 0;
    // for (int i = 0; i < n; ++i) {
    //     double diff = x[i] - x1[i];
    //     if (diff < 0) {
    //         diff *= -1;
    //     }
    //     if (diff > 0.1) {
    //         if (diff / x1[i] > 0.01)
    //         std::cout << i << " " << x1[i] << " " << x[i] << " " << diff << " " << diff / x1[i] << std::endl;
    //     }
    //     // if (x1[i] != x[i]) {
    //     //     std::cout << x[i] << " " << x1[i] << std::endl;
    //     // }
    // }
    // std::cout << count;

    // for (int i = 0; i < n; ++i) {
    //     std::cout << x[i] << " " << x1[i] << std::endl;
    // }
#endif
}