#include <omp.h>
#include <iostream>
#include <chrono>
#include <set>
#include <vector>
#include "readMatrix.h"
#include "makeSet.h"
#include "verify.h"

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cout << "Usage: ./main /path/to/L.mtx /path/to/RHS.mtx" << std::endl;
        return -1;
    }
    std::string fName = argv[1];
    std::string rhsFName = argv[2];
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsed_seconds;
    int *col, *row;
    double  *y, *val, *x, *x1;


    int n, nnz;
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
    //*******************************************************************
    // NAIVE LSOLVE
    start = std::chrono::system_clock::now();
    for (int j = 0; j < n; ++j) {
        x[j] /= val[col[j]];
        int p = col[j] + 1;
        for (; p < col[j+1]; ++p) {
            x[row[p]] -= val[p] * x[j];
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout << "Naive lsolve: " << elapsed_seconds.count() << std::endl;
    verify(col, row, val, n, x, rhsFName);
    //********************************************************************
    // SINGLE-THREADED OPTIMIZED SOLVE
    for (int i = 0; i < n; ++i) { // Reconstructing the RHS matrix
        x[i] = x1[i];
    }
    // Making Reachset (First constructing the DAG)
    bool *visited = new bool[n]();
    std::vector<std::set<int>> dag = makeDAG(col, row, val, n, visited);
    std::set<int> reachset = makeReachset(dag, visited, n, x);

    // Reach set made. Now Solve
    start = std::chrono::system_clock::now();
    for (auto j : reachset) {
        x[j] /= val[col[j]];
        int p = col[j] + 1;
        for (; p < col[j+1]; ++p) {
            x[row[p]] -= val[p] * x[j];
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout << "Optimised lsolve: " << elapsed_seconds.count() << std::endl;
    verify(col, row, val, n, x, rhsFName);
    //*********************************************************************
    //PARALLEL ALGORITHM 1
    for (int i = 0; i < n; ++i) {
        x[i] = x1[i];
    }
    // Store levels of parallelism. Tasks inside same level can be performed parallely; each 2 levels need to sync
    std::vector< std::vector<int> > levels = makeLevelSet(dag, reachset, n);
    
    // Levels created, now parallel solve
    start = std::chrono::system_clock::now();
    for (int i = 0; i < levels.size(); ++i) {
        #pragma omp parallel for
        for (int j = 0; j < levels[i].size(); ++j) { // This loop is parallelized
            int idx = levels[i][j];
            x[idx] /= val[col[idx]];
            int p = col[idx] + 1;
            for (; p < col[idx+1]; ++p) {
                x[row[p]] -= val[p] * x[idx];
            }
        }
    }
    end = std::chrono::system_clock::now();
    elapsed_seconds = end-start;
    std::cout << "Parallel Algorithm 1 lsolve: " << elapsed_seconds.count() << std::endl;
    verify(col, row, val, n, x, rhsFName);
    //*******************************************************************************
    // PARALLEL ALGORITHM 2
    for (int i = 0; i < n; ++i) {
        x[i] = x1[i];
    }
    int currLevel = 0; // store current level
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
    std::cout << "Parallel Algorithm 2 lsolve: " << elapsed_seconds.count() << std::endl;
    verify(col, row, val, n, x, rhsFName);
    //*********************************************************************************
}
