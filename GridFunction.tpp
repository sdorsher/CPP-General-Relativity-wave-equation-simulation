#include "GridFunction.h"


//Constructor. Initializes GridFunction to store Array1Ds consisting
//of initial value of specified size. If vecSize=0 and arraySize is not
//equal to zero, data with the same arraySize can be appended later.
template <class T>
GridFunction<T>::GridFunction(int vecSize, int arraySize,
                                  T initvalue)
                                  :GFvectorDim(vecSize), GFarrayDim(arraySize)
{
  if((vecSize<0) || (arraySize<0)) {
    cout << "Negative grid function dimensions at GridFunction constructor.";
  }
 
  for(int v=0; v < GFvectorDim; v++) {
    std::vector<T> temp(GFarrayDim, initvalue);
    data.push_back(temp);
  }
}


//Constructor with uninitialized values, has same property for vecSize=0 
//as above.
template <class T>
GridFunction<T>::GridFunction(int vecSize, int arraySize)
  :GFvectorDim(vecSize), GFarrayDim(arraySize)
{
  if((vecSize<0) || (arraySize<0)) {
    cout << "Negative grid function dimensions at GridFunction constructor.";
  }

  for(int v=0; v < GFvectorDim; v++) {
    std::vector<T> temp(GFarrayDim);
    data.push_back(temp);
  }
}

template <class T>
GridFunction<T>::~GridFunction(){
}

//Set value at vector coordinate to an array
template <class T>
void GridFunction<T>::set(int vcoord, std::vector<T> arraydata)
{
  if((0>vcoord) || (vcoord>=GFvectorDim)) {
    cout << "Grid coordinate out of range in set(int,Array1D)";
  } else if (arraydata.size()!=GFarrayDim) {
    cout << "Grid function data size does not match.";
  } else {
    data.at(vcoord) = arraydata; 
  }
}

//Set vector coordinate and array coordinate to a value.
template <class T>
void GridFunction<T>::set(int vcoord, int acoord, T value)
{
  if((0>vcoord) || (vcoord>=GFvectorDim)) {
    cout << "Grid coordinate out of range in set(int,int,T)";
  } else if((0>acoord) || (acoord>=GFarrayDim)) {
    cout << "Grid function coordinate out of range in set(int, int,T)";
  } else {
    data.at(vcoord).at(acoord) = value;
  }
}

//Append an array of the same size as the arrays already present.
template <class T>
void GridFunction<T>::append(std::vector<T> array)
{
  if(array.size() != GFarrayDim) {
    cout << GFarrayDim << endl;
    cout << "Incorrect array dimensions in GridFunction::append";
  }

  data.push_back(array);
  GFvectorDim++;
}


//-----------------------------------------
// Not in class

//Addition operator for GridFunctions.
template <typename T>
GridFunction<T> operator+(GridFunction<T> gf1,GridFunction<T> gf2)
{
  if((gf1.GFvecDim() != gf2.GFvecDim()) 
      || (gf1.GFarrDim() != gf2.GFarrDim())) {
    cout << "Grid function dimension mismatch in + operator";
  } else {
    GridFunction<T> gfout(gf1.GFvecDim(), gf1.GFarrDim());
    for(int i = 0; i < gf1.GFvecDim(); i++) {
      for(int j=0; j<gf1.GFarrDim(); j++){
        gfout.set(i,j,gf1.get(i,j) + gf2.get(i,j));
      }
    }
    return gfout;
  }
}

//Multiplication operator for a grid function and a scalar.
template <typename T>
GridFunction<T> operator*(T A, GridFunction<T> gf)
{
  GridFunction<T> gfout(gf.GFvecDim(), gf.GFarrDim());
  for(int i = 0; i < gf.GFvecDim(); i++) {
    for(int j=0; j< gf.GFarrDim(); j++){
      gfout.set(i,j,A * gf.get(i,j));	   
    }
 }
  return gfout;
}

//Multiplication operator for a complex real scalar and a complex GridFunction
template <typename T>
GridFunction<complex<T>> operator*(T A, GridFunction<complex<T>> gf)
{
  GridFunction<complex<T>> gfout(gf.GFvecDim(), gf.GFarrDim());
  for(int i = 0; i < gf.GFvecDim(); i++) {
    for(int j=0; j< gf.GFarrDim(); j++) {
      gfout.set(i,j,A * gf.get(i,j));
    }
  }
  return gfout;
}
