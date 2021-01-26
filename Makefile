CC=g++
CFLAGS=-std=c++11 -O3 -fopenmp

all: main toLower

main: main.cpp makeSet.h readMatrix.h verify.h
	$(CC) $(CFLAGS) main.cpp -o main

toLower: toLower.cpp
	$(CC) toLower.cpp -o toLower
