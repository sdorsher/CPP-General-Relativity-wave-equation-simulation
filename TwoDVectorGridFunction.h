#ifndef TWODVECTORGRIDFUNCTION_H
#define TWODVECTORGRIDFUNCTION_H

#include "VectorGridFunction.h"
#include <fstream>
#include <complex>
#include <iostream>
#include <vector>

using namespace std;

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
  //constructor that takes an initial value

  TwoDVectorGridFunction(int TDVGFSize, int VGFvecSize, int GFvecSize, 
                         int GFarraySize);
  //constructor for empty arrays
  
  inline VectorGridFunction<T> get(int TDVGFvcoord);
  // get the VectorGridFunction at the outer TDVGFcoord specified
  
  inline T get(int TDVGFvcoord, int VGFvcoord, int GFvcoord, int GFacoord);
  // get the value at the full set of four coordinates specified
  
  inline vector<T> get(int TDVGFvcoord, int VGFvcoord,int GFvcoord);
  // get the vector at the outer three vector coordinates specified
  
  inline GridFunction<T> get(int TDVGFvcoord, int VGFvcoord);
  // Get the grid function at the outer two vector coordinates specified
  
  void set(int TDGVFvcoord, int VGFvcoord, int GFvcoord, int GFacoord,T value);
  // Set the TDVGF at the full set of four coordinates specified to the value specified
  
  void set(int TDVGFvcoord, VectorGridFunction<T> vgf);
  //Set the TDVGF at the outer TDVGF vector coordinate to the Vector Grid Function specified.
  
  void set(int TDVGFvcoord, int VGFcoord,GridFunction<T> value);
  // Set the TDVGF at the outer two vector coordinates to the Grid Function specified. 
  
  void set(int TDVGFvcoord, int VGFcoord,int GFcoord,std::vector<T> arr);
  // Set the TDGVGF at the outer three vector coordinates the the vector specified
  
  void setVector(int TDVGFvcoord, int GFcoord, int GFacoord, vector<T> vec);
  // Set the TDVGF using a for loop along the VGF vector dimension to the vector specified at
  // the other three coordinates specified.

  int TDVGFdim(); //the dimension of the outermost vector
  int VGFdim(); // the dimension of the middle vector
  int GFvecDim(); //the dimension of the vector within the GridFunction
  int GFarrDim(); //the dimension of the array within the GridFunction
  void append(VectorGridFunction<T> gf);

  //  vector<T> getVector(int TDVGFvcoord, int GFvcoord, int GFacoord);
  //Using a for loop, get the vector along the VGFvector dimension.

  vector<T> getVectorRange(int vectorCoord, int GFvcoord, int GFacoord,
                                int vmin, int vmax, int dimension);
  //Using for loops, get that same vector as an vector.

  //was Array2D
   vector<T> getVectorNode2D(int vectorCoord, int GFcoord,int startvec, 
			     int stopvec, int dimension, int* dim1, int* dim2);
  //Using for loops, get an array along the either the outermost or second outermost
  //vector dimension and the GFvector dimensions at the other
  // coordinates specified. Choose which of the two dimensions using "dimension". Choose
  // which elements to take from that dimension using startvec and stopvec. 

};


//Adds two TwoDVectorGridFunctions and returns a third
template <typename T>
inline TwoDVectorGridFunction<T> operator+(TwoDVectorGridFunction<T>, TwoDVectorGridFunction<T>);

//Multiplies a real and a TwoDVectorGridFunctions and returns a TwoDVectorGridFunction
template <typename T>
inline TwoDVectorGridFunction<T> operator*(T, TwoDVectorGridFunction<T>);

//Multiplies a real and a complex TwoDVectorGridFunctions and returns a complex
//TwoDVectorGridFunction
template <typename T>
inline TwoDVectorGridFunction<complex<T>> operator*(T, TwoDVectorGridFunction<complex<T>>);


//Get dimension of outer vector.
template <class T>
inline int TwoDVectorGridFunction<T>::TDVGFdim()
{
  return TDVGFvectorDim;
}

//Get dimension of middle vector.
template <class T>
inline int TwoDVectorGridFunction<T>::VGFdim()
{
  return VGFvectorDim;
}

//Get dimension of inner vector.
template <class T>
inline int TwoDVectorGridFunction<T>::GFvecDim()
{
  return GFvectorDim;
}

//Get dimension of array.
template <class T>
inline int TwoDVectorGridFunction<T>::GFarrDim()
{
  return GFarrayDim;
}




#include "TwoDVectorGridFunction.tpp"

#endif
