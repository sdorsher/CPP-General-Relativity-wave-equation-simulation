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
  inline T get(int vcoord, int acoord);
  inline TNT::Array1D<T> get(int vcoord);
  inline int GFvecDim(); //Dimension of vector
  inline int GFarrDim(); //Dimension of Array1D
  void save(string filename);
  
  // private:
 private:
  vector<TNT::Array1D<T>> data;
  int GFvectorDim; 
  int GFarrayDim;
};

//Adds a grid function to a grid function
template <typename T>
GridFunction<T> operator+(GridFunction<T> gf1,GridFunction<T> gf2);


//Multiplies a real by a real grid function
template <typename T>
GridFunction<T> operator*(T A,GridFunction<T> gf2);


//Multiplies a real by a complex grid function
template <typename T>
GridFunction<complex<T>> operator*(T A,GridFunction<complex<T>> gf2);

//Get a value at a vector and array coordinate.
template <class T>
T GridFunction<T>::get(int vcoord, int acoord)
{
  if((0>vcoord) || (vcoord>=GFvectorDim)) {
    throw out_of_range("Grid coordinate out of range in get(int,int)");
  } else if((0>acoord) || (acoord>=GFarrayDim)) {
    throw out_of_range("Grid function coordinate out of range in get(int,int)");
  } else {
    return data.at(vcoord)[acoord];
  }
}


//Get an array at a vector coordinate.
template <class T>
inline TNT::Array1D<T> GridFunction<T>::get(int vcoord)
{
  if((0>vcoord) || (vcoord>=GFvectorDim)) {
    throw out_of_range("Grid coordinate out of range in get(int)");
  } else {
    return data.at(vcoord);
  }
}


//Get the dimension of the vector.
template <class T>
inline int GridFunction<T>::GFvecDim()
{
  return GFvectorDim;
}


//Get the dimension of the array.
template <class T>
inline int GridFunction<T>::GFarrDim()
{
  return GFarrayDim;
}





#include "GridFunction.tpp"

#endif
