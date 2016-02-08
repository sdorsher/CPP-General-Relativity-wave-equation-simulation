#ifndef _GRIDFUNCTION_H
#define _GRIDFUNCTION_H
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <iostream>

using namespace std;


template <class T>
class GridFunction
{
 public:
  //See .cpp file for explanations of functions
  GridFunction(int vecSize, int arraySize,T initvalue);
  GridFunction(int vecSize, int arraySize);
  ~GridFunction();

  void set(int vcoord, std::vector<T>); 
  void set(int vcoord, int acoord,T value); 
  void append(std::vector<T> array);
  inline T get(int vcoord, int acoord);
  inline std::vector<T> get(int vcoord);
  inline int GFvecDim(); //Dimension of std::vector
  inline int GFarrDim(); //Dimension of Array1D
  void save(string filename);
  
  // private:
 private:
  std::vector<std::vector<T>> data;
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
template <typename T>
T GridFunction<T>::get(int vcoord, int acoord)
{
  if((0>vcoord) || (vcoord>=GFvectorDim)) {
    cout <<"Grid coordinate out of range in get(int,int)";
  } else if((0>acoord) || (acoord>=GFarrayDim)) {
    cout << "Grid function coordinate out of range in get(int,int)";
  } else {
    return data.at(vcoord).at(acoord);
  }
}


//Get an array at a vector coordinate.
template <typename T>
inline std::vector<T> GridFunction<T>::get(int vcoord)
{
  if((0>vcoord) || (vcoord>=GFvectorDim)) {
    cout <<"Grid coordinate out of range in get(int)";
  } else {
    return data.at(vcoord);
  }
}


//Get the dimension of the vector.
template <typename T>
inline int GridFunction<T>::GFvecDim()
{
  return GFvectorDim;
}


//Get the dimension of the array.
template <typename T>
inline int GridFunction<T>::GFarrDim()
{
  return GFarrayDim;
}





#include "GridFunction.tpp"

#endif
