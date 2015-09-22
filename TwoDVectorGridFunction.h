#ifndef TWODVECTORGRIDFUNCTION_H
#define TWODVECTORGRIDFUNCTION_H

#include "TNT2.h"
#include "VectorGridFunction.h"
#include <fstream>
#include <complex>

//See .cpp file for explanations of functions
template <class T>
class TwoDVectorGridFunction
{
 public:
  vector<VectorGridFunction<T>> data;
  int TDVGFvectorDim; //outermost vector dimension
  int VGFvectorDim; //middle vector dimension
  int GFvectorDim; //inner vector dimension
  int GFarrayDim; //array dimension

 public:
  TwoDVectorGridFunction(int TDVGFvecSize, int VGFvecSize, int GFvecSize, 
                         int GFarraySize, T initvalue);
  TwoDVectorGridFunction(int TDVGFSize, int VGFvecSize, int GFvecSize, 
                         int GFarraySize);
  VectorGridFunction<T> get(int TDVGFvcoord);
  T get(int TDVGFvcoord, int VGFvcoord, int GFvcoord, int GFacoord);
  Array1D<T> get(int TDVGFvcoord, int VGFvcoord,int GFvcoord);
  GridFunction<T> get(int TDVGFvcoord, int VGFvcoord);
  void set(int TDGVFvcoord, int VGFvcoord, int GFvcoord, int GFacoord,T value);
  void set(int TDVGFvcoord, VectorGridFunction<T> vgf);
  void set(int TDVGFvcoord, int VGFcoord,GridFunction<T> value);
  void set(int TDVGFvcoord, int VGFcoord,int GFcoord,TNT::Array1D<T> arr);
  void setVector(int TDVGFvcoord, int GFcoord, int GFacoord, vector<T> vec);
  int TDVGFdim(); //the dimension of the outermost vector
  int VGFdim(); // the dimension of the middle vector
  int GFvecDim(); //the dimension of the vector within the GridFunction
  int GFarrDim(); //the dimension of the array within the GridFunction
  void append(VectorGridFunction<T> gf);
  vector<T> getVector(int TDVGFvcoord, int GFvcoord, int GFacoord);
  Array1D<T> getVectorAsArray1D(int vectorCoord, int GFvcoord, int GFacoord,
                                int vmin, int vmax, int dimension);
  Array2D<T> getVectorNodeArray2D(int vectorCoord, int GFcoord,int startvec, 
                                  int stopvec, int dimension);
};


//Adds two TwoDVectorGridFunctions and returns a third
template <typename T>
TwoDVectorGridFunction<T> operator+(TwoDVectorGridFunction<T>, TwoDVectorGridFunction<T>);

//Multiplies a real and a TwoDVectorGridFunctions and returns a TwoDVectorGridFunction
template <typename T>
TwoDVectorGridFunction<T> operator*(T, TwoDVectorGridFunction<T>);

//Multiplies a real and a complex TwoDVectorGridFunctions and returns a complex TwoDVectorGridFunction
template <typename T>
TwoDVectorGridFunction<complex<T>> operator*(T, TwoDVectorGridFunction<complex<T>>);

#include "TwoDVectorGridFunction.tpp"

#endif
