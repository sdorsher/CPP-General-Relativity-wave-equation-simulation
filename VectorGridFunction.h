#ifndef VECTORGRIDFUNCTION_H
#define VECTORGRIDFUNCTION_H

#include "TNT2.h"
#include "GridFunction.h"
#include <fstream>
#include <complex>

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
  //Constructor that sets values to default value. 
  
  VectorGridFunction(int VGFvecSize, int GFvecSize, int GFarraySize);
  //Empty constructor.

  inline T get(int VGFvcoord, int GFvcoord, int GFacoord);
  //Get value at full coordinate specification.
  
  inline Array1D<T> get(int VGFvcoord,int GFvcoord);
  // Get Array1D at outer two vector coordinates. 
  
  inline GridFunction<T> get(int VGFvcoord);
  //Get grid function at outermost vectorcoordinate. 
  
  void set(int VGFvcoord, int GFvcoord, int GFacoord,T value);
  //Set value at full coordinate specification
  
  void set(int VGFcoord,GridFunction<T> value);
  // Set Grid Function at outermost vector specification.
  
  void set(int VGFcoord,int GFcoord,TNT::Array1D<T> arr);
  // Set Array1D at outer two vector coordinates.
  
  void setVector(int GFcoord, int GFacoord, vector<T> vec);
  //Using a for loop, set the values along the VGFcoord vector dimension at the inner two
  //coordinates.
  
  int VGFdim();//the dimension of the external vector
  int GFvecDim();//the dimension of the vector within the GridFunction
  int GFarrDim();//the dimension of the array within the GridFunction
  void append(GridFunction<T> gf); //appends a GridFunction to the end of the VGF.
  void save(string filenames); //vector of filenames to print to
  
  vector<T> getVector(int GFvcoord, int GFacoord);
  //Using a for loop, returns a vector of the values along the VGFcoord dimension at the
  // inner two coordinates specified.
  
  Array1D<T> getVectorAsArray1D(int GFvcoord, int GFacoord,int vmin, int vmax);
  // Does the same, but returns an Array1D.

  Array2D<T> getVectorNodeArray2D(int GFcoord,int startvec, int stopvec);
  // Returns an Array2D of outer vector (VGFcoord dimension) values and inner coordinate
  //(GFacoord dimension) values. VGFcoord runs from startvec to stopvec. 

};

//Adds two VGFs.
template <typename T>
VectorGridFunction<T> operator+(VectorGridFunction<T>, VectorGridFunction<T>);

//Multiplies a real scalar by a real VGF.
template <typename T>
VectorGridFunction<T> operator*(T, VectorGridFunction<T>);

//Multiplies a real scalar by a complex VGF to return a complex VGF.
template <typename T>
VectorGridFunction<complex<T>> operator*(T, VectorGridFunction<complex<T>>);


//Get dimension of outer vector.
template <class T>
inline int VectorGridFunction<T>::VGFdim()
{
  return VGFvectorDim;
}

//Get dimension of inner vector.
template <class T>
inline int VectorGridFunction<T>::GFvecDim()
{
  return GFvectorDim;
}

//Get dimension of array.
template <class T>
inline int VectorGridFunction<T>::GFarrDim()
{
  return GFarrayDim;
}




#include "VectorGridFunction.tpp"

#endif
