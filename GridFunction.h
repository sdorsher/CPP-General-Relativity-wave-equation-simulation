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
  GridFunction(int vecSize, int arraySize,T initvalue);
  GridFunction(int vecSize, int arraySize);
  void set(int vcoord, TNT::Array1D<T>);
  void set(int vcoord, int acoord,T value);
  void append(TNT::Array1D<T> array);
  T get(int vcoord, int acoord);
  TNT::Array1D<T> get(int vcoord);
  int gridDim();
  int pointsDim();
  void save(string filename);
  
  // private:
 private:
  vector<TNT::Array1D<T>> data;
  int GFvectorDim;
  int GFarrayDim;

  

//  void initFromFile(string filename);






//  GridFunction(GridFunction&&);
//  GridFunction& operator=(GridFunction&&);
  //need move constructors for addition of RHS and boundary conditions
//  GridFunction(const GridFunction&);
//  GridFunction& operator=(const GridFunction&);
  //need copy constructors for intermediate steps of time evolution
 
};

template <typename T>
GridFunction<T> operator+(GridFunction<T> gf1,GridFunction<T> gf2);
template <typename T>
GridFunction<T> operator*(T A,GridFunction<T> gf2);

#include "GridFunction.tpp"

#endif
