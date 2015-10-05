#include "GridFunction.h"

using namespace std;


//Constructor. Initializes GridFunction to store Array1Ds consisting
//of initial value of specified size. If vecSize=0 and arraySize is not
//equal to zero, data with the same arraySize can be appended later.
template <class T>
GridFunction<T>::GridFunction(int vecSize, int arraySize,
                                  T initvalue)
                                  :GFvectorDim(vecSize), GFarrayDim(arraySize)
{
  if((vecSize<0) || (arraySize<0)) {
    throw invalid_argument("Negative grid function dimensions at GridFunction constructor.");
  }
 
  for(int v=0; v < GFvectorDim; v++) {
    TNT::Array1D<T> temp(GFarrayDim, initvalue);
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
    throw invalid_argument("Negative grid function dimensions at GridFunction constructor.");
  }

  for(int v=0; v < GFvectorDim; v++) {
    TNT::Array1D<T> temp(GFarrayDim);
    data.push_back(temp);
  }
}

//Set value at vector coordinate to an array
template <class T>
void GridFunction<T>::set(int vcoord, TNT::Array1D<T> arraydata)
{
  if((0>vcoord) || (vcoord>=GFvectorDim)) {
    throw out_of_range("Grid coordinate out of range in set(int,Array1D)");
  } else if (arraydata.dim()!=GFarrayDim) {
    throw invalid_argument("Grid function data size does not match.");
  } else {
    data.at(vcoord) = arraydata.copy(); 
    //copy to avoid shallow copy memory problem in TNT Arrays
  }
}

//Set vector coordinate and array coordinate to a value.
template <class T>
void GridFunction<T>::set(int vcoord, int acoord, T value)
{
  if((0>vcoord) || (vcoord>=GFvectorDim)) {
    throw out_of_range("Grid coordinate out of range in set(int,int,T");
  } else if((0>acoord) || (acoord>=GFarrayDim)) {
    throw out_of_range("Grid function coordinate out of range in set(int, int,T)");
  } else {
    data.at(vcoord)[acoord] = value;
  }
}

//Append an array of the same size as the arrays already present.
template <class T>
void GridFunction<T>::append(TNT::Array1D<T> array)
{
  if(array.dim() != GFarrayDim) {
    cout << GFarrayDim << endl;
    throw invalid_argument("Incorrect array dimensions in GridFunction::append");
  }

  data.push_back(array.copy());
  GFvectorDim++;
}


//May be obsolete.
template <class T>
void GridFunction<T>::save(string filename)
{
  ofstream fs;
  fs.open(filename);
  for(int i=0;i<GFvectorDim;i++){
    for(int j=0;j<GFarrayDim;j++){
      fs << data.at(i)[j] <<endl;
    }
  }
  fs.close();
}

//-----------------------------------------
// Not in class

//Addition operator for GridFunctions.
template <typename T>
GridFunction<T> operator+(GridFunction<T> gf1,GridFunction<T> gf2)
{
  if((gf1.GFvecDim() != gf2.GFvecDim()) 
      || (gf1.GFarrDim() != gf2.GFarrDim())) {
    throw invalid_argument("Grid function dimension mismatch in + operator");
  } else {
    GridFunction<T> gfout(0, gf1.GFarrDim());
    for(int i = 0; i < gf1.GFvecDim(); i++) {
      gfout.append(gf1.get(i) + gf2.get(i));
    }
    return gfout;
  }
}

//Multiplication operator for a grid function and a scalar.
template <typename T>
GridFunction<T> operator*(T A, GridFunction<T> gf)
{
  GridFunction<T> gfout(0, gf.GFarrDim());
  for(int i = 0; i < gf.GFvecDim(); i++) {
    gfout.append(A * gf.get(i));
  }
  return gfout;
}

//Multiplication operator for a complex real scalar and a complex GridFunction
template <typename T>
GridFunction<complex<T>> operator*(T A, GridFunction<complex<T>> gf)
{
  GridFunction<complex<T>> gfout(0, gf.GFarrDim());
  for(int i = 0; i < gf.GFvecDim(); i++) {
    gfout.append(A * gf.get(i));
  }
  return gfout;
}
