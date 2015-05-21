#ifndef VECTORGRIDFUNCTION_H
#define VECTORGRIDFUNCTION_H

#include "TNT2.h"
#include "GridFunction.h"
#include <fstream>


template <class T>
class VectorGridFunction
{
 public:
  vector<GridFunction<T>> data;
  int VGFvectorDim;
  int GFarrayDim;
  int GFvectorDim;

 public:
  VectorGridFunction(int VGFvecSize, int GFvecSize, int GFarraySize, T initvalue);
  VectorGridFunction(int VGFvecSize, int GFvecSize, int GFarraySize);
  T get(int VGFvcoord, int GFvcoord, int GFacoord);
  Array1D<T> get(int VGFvcoord,int GFvcoord);
  GridFunction<T> get(int VGFvcoord);
  void set(int VGFvcoord, int GFvcoord, int GFacoord,T value);
  void set(int VGFcoord,GridFunction<T> value);
  void set(int VGFcoord,int GFcoord,TNT::Array1D<T> arr);
  void setVector(int GFcoord, int GFacoord, vector<T> vec);
  int vectorDim();//the dimension of the external vector
  int gridDim();//the dimension of the vector within the GridFunction
  int pointsDim();//the dimension of the array within the GridFunction
  void append(GridFunction<T> gf);
  void save(string filenames); //vector of filenames to print to
  //  void initFromFile(string filename);
  vector<T> getVector(int GFvcoord, int GFacoord);
  
};

template <typename T>
VectorGridFunction<T> operator+(VectorGridFunction<T>, VectorGridFunction<T>);
  //need addition operator for addition of RHS and boundary conditions

template <typename T>
VectorGridFunction<T> operator*(T, VectorGridFunction<T>);
//for easy multiplication in rk4 routine

#include "VectorGridFunction.tpp"

#endif
