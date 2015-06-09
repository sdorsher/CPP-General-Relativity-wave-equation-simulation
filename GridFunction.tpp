#include "GridFunction.h"

using namespace std;

template <class T>
GridFunction<T>::GridFunction(int vecSize, int arraySize,
                                  T initvalue)
                                  :GFvectorDim(vecSize),GFarrayDim(arraySize)
{
  if((vecSize<0)||(arraySize<0)){
    throw invalid_argument("Negative grid function dimensions at GridFunction constructor.");
  }

  for(int v=0; v<GFvectorDim; v++){
    TNT::Array1D<T> temp(GFarrayDim,initvalue);
    data.push_back(temp);
  }
}

template <class T>
GridFunction<T>::GridFunction(int vecSize, int arraySize)
  :GFvectorDim(vecSize),GFarrayDim(arraySize)
{
  if((vecSize<0)||(arraySize<0)){
    throw invalid_argument("Negative grid function dimensions at GridFunction constructor.");
  }

  for(int v=0; v<GFvectorDim; v++){
    TNT::Array1D<T> temp(GFarrayDim);
    data.push_back(temp);
  }
}

template <class T>
void GridFunction<T>::set(int vcoord, TNT::Array1D<T> arraydata)
{
  if((0>vcoord)||(vcoord>=GFvectorDim)){
    throw out_of_range("Grid coordinate out of range in set(int,Array1D)");
  } else if (arraydata.dim()!=GFarrayDim) {
    throw invalid_argument("Grid function data size does not match.");
  } else {
    data.at(vcoord)=arraydata.copy();
  }
}

template <class T>
void GridFunction<T>::set(int vcoord, int acoord,T value)
{
  if((0>vcoord)||(vcoord>=GFvectorDim)){
    throw out_of_range("Grid coordinate out of range in set(int,int,T");
  } else if((0>acoord)||(acoord>=GFarrayDim)) {
    throw out_of_range("Grid function coordinate out of range in set(int, int,T)");
  } else {
    data.at(vcoord)[acoord]=value;
  }
}

template <class T>
void GridFunction<T>::append(TNT::Array1D<T> array)
{
  if(array.dim()!=GFarrayDim){
    cout << GFarrayDim << endl;
    throw invalid_argument("Incorrect array dimensions in GridFunction::append");
  }

  data.push_back(array.copy());
  GFvectorDim++;
}

template <class T>
T GridFunction<T>::get(int vcoord, int acoord)
{
  if((0>vcoord)||(vcoord>=GFvectorDim)){
    throw out_of_range("Grid coordinate out of range in get(int,int)");
  } else if((0>acoord)||(acoord>=GFarrayDim)) {
    throw out_of_range("Grid function coordinate out of range in get(int,int)");
  } else {
    return data.at(vcoord)[acoord];
  }
}

template <class T>
TNT::Array1D<T> GridFunction<T>::get(int vcoord)
{
  if((0>vcoord)||(vcoord>=GFvectorDim)){
    throw out_of_range("Grid coordinate out of range in get(int)");
  } else {
    return data.at(vcoord);
  }
}

template <class T>
int GridFunction<T>::gridDim()
{
  return GFvectorDim;
}

template <class T>
int GridFunction<T>::pointsDim()
{
  return GFarrayDim;
}

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

template <typename T>
GridFunction<T> operator+(GridFunction<T> gf1,GridFunction<T> gf2)
{
  if((gf1.gridDim()!=gf2.gridDim())||(gf1.pointsDim()!=gf2.pointsDim())){
    throw invalid_argument("Grid function dimension mismatch in + operator");
  } else {
    GridFunction<T> gfout(0,gf1.pointsDim());
    for(int i=0; i<gf1.gridDim(); i++) {
      gfout.append(gf1.get(i)+gf2.get(i));
    }
    return gfout;
  }
}

template <typename T>
GridFunction<T> operator*(T A, GridFunction<T> gf)
{
  GridFunction<T> gfout(0,gf.pointsDim());
  for(int i=0; i<gf.gridDim(); i++) {
    gfout.append(A*gf.get(i));
  }
  return gfout;
}
