#ifndef _GRIDFUNCTION_H
#define _GRIDFUNCTION_H
#include "tnt/tnt.h"
#include "TNT2.h"
#include <fstream>

using namespace std;

template <class T>
class GridFunction
{
 public:
  //See .cpp file for explanations of functions
  GridFunction(int vecSize, int arraySize,T initvalue);
  GridFunction(int vecSize, int arraySize);
  void set(int vcoord, TNT::Array1D<T>); 
  void set(int vcoord, int acoord,T value); 
  void append(TNT::Array1D<T> array);
  T get(int vcoord, int acoord);
  TNT::Array1D<T> get(int vcoord);
  int gridDim(); //Dimension of vector
  int pointsDim(); //Dimension of Array1D
  void save(string filename);
  
  // private:
 private:
  vector<TNT::Array1D<T>> data;
  int GFvectorDim; 
  int GFarrayDim;
};

template <typename T>
GridFunction<T> operator+(GridFunction<T> gf1,GridFunction<T> gf2);

template <typename T>
GridFunction<T> operator*(T A,GridFunction<T> gf2);

template <typename T>
GridFunction<complex<T>> operator*(T A,GridFunction<complex<T>> gf2);

#include "GridFunction.tpp"

#endif
