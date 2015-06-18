#ifndef TWODVECTORGRIDFUNCTION_H
#define TWODVECTORGRIDFUNCTION_H

#include "TNT2.h"
#include "VectorGridFunction.h"
#include <fstream>

//See .cpp file for explanations of functions
template <class T>
class TwoDVectorGridFunction
{
 public:
  vector<VectorGridFunction<T>> data;
  int TDVGFvectorDim; //outermost vector dimension
  int VGFvectorDim; //middle vector dimension
  int GFarrayDim; //array dimension
  int GFvectorDim; //inner vector dimension

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
  int modesDim(); //the dimension of the outermost vector
  int vectorDim();//the dimension of the middle vector
  int gridDim();//the dimension of the vector within the GridFunction
  int pointsDim();//the dimension of the array within the GridFunction
  void append(VectorGridFunction<T> gf);
  vector<T> getVector(int TDVGFvcoord, int GFvcoord, int GFacoord);
  Array1D<T> getVectorAsArray1D(int vectorCoord, int GFvcoord, int GFacoord,
                                int vmin, int vmax, int dimension);
  Array2D<T> getVectorNodeArray2D(int vectorCoord, int GFcoord,int startvec, 
                                  int stopvec, int dimension);
};

template <typename T>
VectorGridFunction<T> operator+(VectorGridFunction<T>, VectorGridFunction<T>);

template <typename T>
VectorGridFunction<T> operator*(T, VectorGridFunction<T>);

#include "TwoDVectorGridFunction.tpp"

#endif
