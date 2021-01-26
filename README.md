# Optimizing Sparse Matrix Kernels
This repository contains optimized single-threaded and multithreaded sparse-matrix triangular solve. To build the project, run:
```
make
```

To run the project, first ensure that the input matrix is a lower triangular matrix. The project contains a program to print the lower half of a matrix in the following way:
```
./toLower /path/to/input.mtx > /path/to/output.mtx
```

Finally, to run the project:
```
./main /path/to/L.mtx /path/to/RHS.mtx
```
