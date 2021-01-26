#include <omp.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <chrono>
#include <sstream>
#include <set>
#include <vector>
#include <assert.h>

/*
 * reading a CSC matrix from a coordinate file, stored col-ordered
 */
bool readSparseMatrix(std::string fName, int &n, int &NNZ, int* &col,//FIXME change col type to size_t
                int* &row, double* &val){
 /*This function reads the input matrix from "fName" file and
  * allocate memory for matrix A, L and U.
  * - The input file is a coordinate version and e
  * ach row of the file shows (col, row, nnz)
  * - The matrices are zero-indexed
  */

 std::ifstream inFile;
 inFile.open(fName);
 std::string line,banner, mtx, crd, arith, sym;
 /*  File format:
  *    %%MatrixMarket matrix coordinate real general/symmetric/...
  *    % ...
  *    % (optional comments)
  *    % ...
  *    #rows    #non-zero
  *    Triplet in the rest of lines: row    col    value
  */
 std::getline(inFile,line);
 for (unsigned i=0; i<line.length(); line[i]=tolower(line[i]),i++);
 std::istringstream iss(line);
 if (!(iss >> banner >> mtx >> crd >> arith >> sym)){
  std::cout<<"Invalid header (first line does not contain 5 tokens)" << std::endl;
  return false;
 }

 if(banner.compare("%%matrixmarket")) {
  std::cout<<"Invalid header (first token is not \"%%%%MatrixMarket\")" << std::endl;
  return false;
 }
 if(mtx.compare("matrix")) {
  std::cout<<"Not a matrix; this driver cannot handle that." << std::endl;
  return false;
 }
 if(crd.compare("coordinate")) {
  std::cout<<"Not in coordinate format; this driver cannot handle that." << std::endl;
  return false;
 }
 if(arith.compare("real")) {
  if(!arith.compare("complex")) {
   std::cout<<"Complex matrix; use zreadMM instead!" << std::endl;
   return false;
  }
  else if(!arith.compare("pattern")) {
   std::cout<<"Pattern matrix; values are needed!" << std::endl;
   return false;
  }
  else {
   std::cout<<"Unknown arithmetic" << std::endl;
   return false;
  }
 }
 while (!line.compare(0,1,"%"))
 {
  std::getline(inFile, line);
 }
 std::istringstream issDim(line);
 if (!(issDim >> n >> n >> NNZ)){
  std::cout<<"The matrix dimension is missing" << std::endl;
  return false;
 }
 if(n <= 0 || NNZ <= 0)
  return false;
 col = new int[n + 1]();
 // colL = new int[n + 1]; colU = new int[n + 1];
 row = new int[NNZ];
 // rowL = new int[factorSize]; rowU = new int[factorSize];
 val = new double[NNZ];
 // valL = new double[factorSize]; valU = new double[factorSize];
 if(!val || !col || !row)
  return false;
 //Initializing the result vector
 int y, x, colCnt=0, nnzCnt=0;
 double value;

 col[0]=0;
 for (int i = 0; i < NNZ; ++i) {
    inFile>>x;//zero indexing
    inFile>>y;
    inFile>>value;
    if(y > n || y > x || x > n) {
        std::cout << x << " " << y << " " << value ;
        return false;
    }
    val[i] = value;
    row[i] = x - 1;
    ++col[y];
 }
 for (int i = 2; i <= n; ++i)
    col[i] += col[i-1];
 return true;
}

bool readRHSMatrix(std::string fName, double* &val, const int size) {

 std::ifstream inFile;
 inFile.open(fName);
 std::string line,banner, mtx, crd, arith, sym;
 /*  File format:
  *    %%MatrixMarket matrix coordinate real general/symmetric/...
  *    % ...
  *    % (optional comments)
  *    % ...
  *    #rows    #non-zero
  *    Triplet in the rest of lines: row    col    value
  */
 std::getline(inFile,line);
 for (unsigned i=0; i<line.length(); line[i]=tolower(line[i]),i++);
 std::istringstream iss(line);
 if (!(iss >> banner >> mtx >> crd >> arith >> sym)){
  std::cout<<"Invalid header (first line does not contain 5 tokens)" << std::endl;
  return false;
 }

 if(banner.compare("%%matrixmarket")) {
  std::cout<<"Invalid header (first token is not \"%%%%MatrixMarket\")" << std::endl;
  return false;
 }
 if(mtx.compare("matrix")) {
  std::cout<<"Not a matrix; this driver cannot handle that." << std::endl;
  return false;
 }
 if(crd.compare("coordinate")) {
  std::cout<<"Not in coordinate format; this driver cannot handle that." << std::endl;
  return false;
 }
 if(arith.compare("real")) {
  if(!arith.compare("complex")) {
   std::cout<<"Complex matrix; use zreadMM instead!" << std::endl;
   return false;
  }
  else if(!arith.compare("pattern")) {
   std::cout<<"Pattern matrix; values are needed!" << std::endl;
   return false;
  }
  else {
   std::cout<<"Unknown arithmetic" << std::endl;
   return false;
  }
 }
 while (!line.compare(0,1,"%"))
 {
  std::getline(inFile, line);
 }
 int n, cols, NNZ;
 std::istringstream issDim(line);
 if (!(issDim >> n >> cols >> NNZ)){
  std::cout<<"The matrix dimension is missing" << std::endl;
  return false;
 }
 if(n <= 0 || NNZ <= 0)
  return false;
 if (n != size) {
     std::cout << "Dimensions of L and RHS do not match" << std::endl;
     std::cout << size << " " << n << std::endl;
     return false;
 }
 if (cols != 1) {
     std::cout << "RHS matrix has incorrect dimensions" << std::endl;
     return false;
 }
 val = new double[n];
 // valL = new double[factorSize]; valU = new double[factorSize];
 if(!val)
  return false;
 //Initializing the result vector
 int y, x;
 double value;

 for (int i = 0; i < NNZ; ++i) {
    inFile>>x;//zero indexing
    inFile>>y;
    inFile>>value;
    if(y > 1 || x > n) {
        std::cout << x << " " << y << " " << value ;
        return false;
    }
    val[x-1] = value;
 }
 return true;
}

bool printLower(std::string fName){
 /*This function reads the input matrix from "fName" file and
  * allocate memory for matrix A, L and U.
  * - The input file is a coordinate version and e
  * ach row of the file shows (col, row, nnz)
  * - The matrices are zero-indexed
  */

 std::ifstream inFile;
 double tol=0.1;
 size_t n, NNZ;
 inFile.open(fName);
 std::string line,banner, mtx, crd, arith, sym;
 /*  File format:
  *    %%MatrixMarket matrix coordinate real general/symmetric/...
  *    % ...
  *    % (optional comments)
  *    % ...
  *    #rows    #non-zero
  *    Triplet in the rest of lines: row    col    value
  */
 std::getline(inFile,line);
 //std::cout<<line<<std::endl;
 std::cout<<"%%MatrixMarket matrix coordinate real symmetric" << std::endl;
 for (unsigned i=0; i<line.length(); line[i]=tolower(line[i]),i++);
 std::istringstream iss(line);
 if (!(iss >> banner >> mtx >> crd >> arith >> sym)){
  std::cout<<"Invalid header (first line does not contain 5 tokens)" << std::endl;
  return false;
 }

 if(banner.compare("%%matrixmarket")) {
  std::cout<<"Invalid header (first token is not \"%%%%MatrixMarket\")" << std::endl;
  return false;
 }
 if(mtx.compare("matrix")) {
  std::cout<<"Not a matrix; this driver cannot handle that." << std::endl;
  return false;
 }
 if(crd.compare("coordinate")) {
  std::cout<<"Not in coordinate format; this driver cannot handle that." << std::endl;
  return false;
 }
 if(arith.compare("real")) {
  if(!arith.compare("complex")) {
   std::cout<<"Complex matrix; use zreadMM instead!" << std::endl;
   return false;
  }
  else if(!arith.compare("pattern")) {
   std::cout<<"Pattern matrix; values are needed!" << std::endl;
   return false;
  }
  else {
   std::cout<<"Unknown arithmetic" << std::endl;
   return false;
  }
 }
 while (!line.compare(0,1,"%"))
 {
  std::getline(inFile, line);
  //std::cout<<line<<std::endl;
 }
 std::istringstream issDim(line);
 if (!(issDim >> n >> n >> NNZ)){
  std::cout<<"The matrix dimension is missing" << std::endl;
  return false;
 }
 if(n <= 0 || NNZ <= 0)
  return false;
 //Initializing the result vector
 int y, x, colCnt=0, nnzCnt=0;
 std::string value;
 std::vector<std::pair<std::pair<int, int>, std::string>> contents; // storing row, col and val
 for (int i = 0; i < NNZ; ++i) {
    inFile>>x;
    inFile>>y;
    inFile>>value;
    if (x >= y) {
        contents.push_back(std::make_pair(std::make_pair(x, y), value));
        ++nnzCnt;
    }
 }
 std::cout << n << " " << n << " " << nnzCnt << std::endl;
 for (auto f : contents) {
    std::cout << f.first.first << " " << f.first.second << " " << f.second << std::endl;
 }
 return true;
}

void visit(std::vector<std::set<int>> const &dag, int idx, bool visited[], std::set<int>& reachset) {
    reachset.insert(idx);
    visited[idx] = true;
    for (auto f : dag[idx]) {
        if (!visited[f])
            visit(dag, f, visited, reachset);
    }
}

void insertToLevel(std::vector< std::vector<int> >& levels, int const node, int const currLevel) {
    if (levels.size() <= currLevel) {
        levels.emplace_back(); // construct a new level
    }
    levels[currLevel].push_back(node);
}

// Function that returns the DAG. visited matrix is later used to generate reachset
std::vector<std::set<int>> makeDAG(int* &col, int* &row, double* &val, int n, bool* &visited) {
    std::vector<std::set<int>> dag(n);
    for (int i = 0; i < n; ++i) {
        visited[i] = false;
    }
    

    // Making the dag
    for (int i = 1; i < n; ++i) {
        for (int j = col[i-1]; j < col[i]; ++j) {
            if (row[j] != i - 1) // Do not create a self-loop
                dag[i-1].insert(row[j]);
        }
    }
    return dag;
}

std::set<int> makeReachset(std::vector<std::set<int>> const &dag, bool* &visited, int n, double* &x) {
    std::set<int> reachset;

    // Creating dependencies using DFS
    for (int i = 0; i < n; ++i) {
        if (x[i] != 0) {
            if (!visited[i]) {
                visit(dag, i, visited, reachset);
            }
        }
    }
    return reachset;
}

// Verify if L*x = b
int verify(int* &col, int* &row, double* &val, int n, double* &x, std::string rhsFName) {
    double *res = new double[n]();

    // L*x sparse-matrix vector multiplcation
    for (int i = 0; i < n; ++i) {
        for (int j = col[i]; j < col[i+1]; ++j) {
            res[row[j]] += val[j] * x[i];
        }
    }

    double *x1;
    if (!readRHSMatrix(rhsFName, x1, n)) {
        std::cout << "Error while reading RHS matrix." << std::endl;
        return -1;
    }

    int count = 0;
    for (int i = 0; i < n; ++i) {
        double diff = x1[i] - res[i];
        if (diff < 0) {
            diff *= -1;
        }
        if (diff > 0.001) {
            ++count;
            std::cout << i << " " << res[i] << " " << x1[i] << std::endl;
        }
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
    

    // for (auto f : reachset) {
    //     std::cout << f+1 << ": ";
    //     for (auto j : dag[f]) {
    //         std::cout << j+1 << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << "here";
    // exit(0);

#if 1
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
#if 0
    // Creating levels and tasks
    std::vector<int> indegree(n); // Keep a count of incoming edges to record root nodes
    std::vector< std::vector<int> > levels; // Store levels of parallelism. Tasks inside same level can be performed parallely; each 2 levels need to sync
    int currLevel = 0;
    // levels.emplace_back();

    for (auto f : reachset) {
        for (auto j : dag[f]) {
            ++indegree[j];
        }
    }

    while (!reachset.empty()) {
        for (auto f : reachset) {
            if (indegree[f] == 0) {
                insertToLevel(levels, f, currLevel);
            }
        }

        for (auto f : levels[currLevel]) {
            reachset.erase(f);
            for (auto j : dag[f]) {
                --indegree[j];
            }
        }
        ++currLevel;
    }
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

    currLevel = 0;
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
    std::cout << "Parallel lsolve: " << elapsed_seconds.count() << std::endl;

    for (int j = 0; j < n; ++j) {
        // std::cout << f + 1 << ",";
        x1[j] /= val[col[j]];
        int p = col[j] + 1;
        for (; p < col[j+1]; ++p) {
            x1[row[p]] -= val[p] * x1[j];
        }
    }
    std::cout << verify(col, row, val, n, x, rhsFName) << std::endl;
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


    // for (int i = 0; i < levels.size(); ++i) {
    //     std::cout << i+1 << ": ";
    //     for (auto j : levels[i]) {
    //         std::cout << j+1 << ", ";
    //     }
    //     std::cout << std::endl;
    // }
}