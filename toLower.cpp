#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>

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

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cout << "Usage: ./toLower /path/to/input.mtx > /path/to/output.mtx" << std::endl;
        return -1;
    }
    std::string fName = argv[1];
    printLower(fName);
}
