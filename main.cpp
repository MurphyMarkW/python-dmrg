#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "dmtk/dmtk.h"

// Debug statements?
//#define DEBUG

// Debug subroutine.
#ifdef DEBUG
  #define Debug(x) std::cout << x << std::endl
#else
  #define Debug(x)
#endif//DEBUG

// Use dmrg or array?
#define USE_DMRG
//#define USE_ARRAY

void dmrgMat(unsigned int size) {
  dmtk::Matrix<double> mat(size,size);
  for(size_t i = 0; i < mat.size1(); ++i) {
    for(size_t j = 0; j < mat.size2(); ++j) {
      mat(i,j) = std::rand();
      Debug(mat(j,i));
    }
  }
  mat = ((mat + 5.0) * 10.0);
};

void arrayMat(unsigned int size) {
  vector<double> mat(size*size,0.0);
  for(size_t i = 0; i < size; ++i) {
    for(size_t j = 0; j < size; ++j) {
      mat[i*size+j] = std::rand();
      Debug(mat[i*size+j]);
    }
  }
  
  for(size_t i = 0; i < size; ++i) {
    for(size_t j = 0; j < size; ++j) {
      mat[i*size+j] = (mat[i*size+j] + 5.0) * 10.0;
    }
  }
  
};

int main() {
  std::srand(std::time(NULL));
  unsigned int size = 16000;

  #ifdef USE_DMRG
  dmrgMat(size);
  #endif//USE_DMRG

  #ifdef USE_ARRAY
  arrayMat(size);
  #endif//USE_ARRAY
}
