#ifndef VECTORGRIDFUNCTION_H
#define VECTORGRIDFUNCTION_H

#include "TNT2.h"
#include "GridFunction.h"
#include <fstream>

//See .cpp file for explanations of functions
template <class T>
class VectorGridFunction
{
 public:
  vector<GridFunction<T>> data;
  int VGFvectorDim; //outer vector dimension
  int GFarrayDim; //array dimension
  int GFvectorDim; //inner vector dimension

 public:
  VectorGridFunction(int VGFvecSize, int GFvecSize, int GFarraySize, 
                     T initvalue);
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
  vector<T> getVector(int GFvcoord, int GFacoord);
  Array1D<T> getVectorAsArray1D(int GFvcoord, int GFacoord,int vmin, int vmax);
  Array2D<T> getVectorNodeArray2D(int GFcoord,int startvec, int stopvec);
};

template <typename T>
VectorGridFunction<T> operator+(VectorGridFunction<T>, VectorGridFunction<T>);

template <typename T>
VectorGridFunction<T> operator*(T, VectorGridFunction<T>);

#include "VectorGridFunction.tpp"

#endif
