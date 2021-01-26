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
bool readSparseMatrix(std::string fName, int &n, int &NNZ, int* &col, int* &row, double* &val){
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

/*
 * reading a dense vector from a coordinate file, stored col-ordered
 */
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
