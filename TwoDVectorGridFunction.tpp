#include "TwoDVectorGridFunction.h"

//Constructor that initializes to a specific value. If inner vector size and 
//array size are fixed, outer vector size may be zero and additional 
//GridFunctions may be appended later.
template <class T>
TwoDVectorGridFunction<T>::TwoDVectorGridFunction(int outerVecSize, 
                                                  int middleVecSize, 
                                                  int innerVecSize, 
                                                  int arraySize, 
                                                  T initvalue):
  TDVGFvectorDim(outerVecSize), 
  VGFvectorDim(middleVecSize), 
  GFvectorDim(innerVecSize), 
  GFarrayDim(arraySize)
{
  if(outerVecSize<0) {
    throw invalid_argument("Negative Vector dimensions at TwoDVectorGridFunction constructor.");
  }
  
  for(int i=0; i<outerVecSize; i++) {
    VectorGridFunction<T> temp(middleVecSize,innerVecSize,arraySize,initvalue);
    data.push_back(temp);
  }
}

//Uninitialized constructor. Same property for outerVecSize=0 as above.
template <class T>
TwoDVectorGridFunction<T>::TwoDVectorGridFunction(int outerVecSize, 
                                                  int middleVecSize,
                                                  int innerVecSize, 
                                                  int arraySize):
  TDVGFvectorDim(outerVecSize),  
  VGFvectorDim(middleVecSize), 
  GFvectorDim(innerVecSize), 
  GFarrayDim(arraySize)
{
  if(outerVecSize<0){
    throw invalid_argument("Negative Vector dimensions at TwoDVectorGridFunction constructor.");
  }
  
  for(int i = 0; i < outerVecSize; i++){
    VectorGridFunction<T> temp(middleVecSize,innerVecSize,arraySize);
    data.push_back(temp);
  }
}

//Get dimension of outer vector.
template <class T>
int TwoDVectorGridFunction<T>::modesDim()
{
  return TDVGFvectorDim;
}

//Get dimension of middle vector.
template <class T>
int TwoDVectorGridFunction<T>::vectorDim()
{
  return VGFvectorDim;
}

//Get dimension of inner vector.
template <class T>
int TwoDVectorGridFunction<T>::gridDim()
{
  return GFvectorDim;
}

//Get dimension of array.
template <class T>
int TwoDVectorGridFunction<T>::pointsDim()
{
  return GFarrayDim;
}

//Set the value of a specific point. 
template <class T>
void TwoDVectorGridFunction<T>::set(int TDVGFvcoord, int VGFvcoord, 
                                    int GFvcoord, int GFacoord,
                                    T value)
{
  if((0 > TDVGFvcoord) || (TDVGFvcoord >= TDVGFvectorDim)){
    throw invalid_argument("2D vector index out of range in set");
  }
  
  data.at(TDVGFvcoord).set(VGFvcoord, GFvcoord, GFacoord, value);
}

//Set an outer and middle vector coordinate to a grid function.
template <class T>
void TwoDVectorGridFunction<T>::set(int TDVGFvcoord, int VGFcoord, 
                                GridFunction<T> gf)
{
  if((0 > TDVGFvcoord) || (TDVGFvcoord >= TDVGFvectorDim)){
    throw invalid_argument("2D vector index out of range in set");
  }

  data.at(TDVGFvcoord).set(VGFcoord,gf);
}

//Set an outer, middle, and inner vector coordinate to an Array1D.
template <class T>
void TwoDVectorGridFunction<T>::set(int TDVGFvcoord, int VGFcoord, 
                                    int GFcoord, TNT::Array1D<T> arr)
{

  if((0 >TDVGFvcoord) || (TDVGFvcoord >= TDVGFvectorDim)){
    throw invalid_argument("2D vector index out of range in set");
  }

  data.at(TDVGFvcoord).set(VGFcoord,GFcoord,arr);
}

//Set outer vector coordinate and inner vector coordinate and an array 
//coordinate to a vector of values.
template <class T>
void TwoDVectorGridFunction<T>::setVector(int TDVGFvcoord, int GFcoord, 
                                          int GFacoord, vector<T> vec)
{
  if((GFcoord < 0) || (GFacoord < 0) || (GFacoord > GFarrayDim)
     || (TDVGFvcoord < 0 ) || (GFcoord > GFvectorDim)){
    throw invalid_argument("Coordinates in setVector out of range.");
  }
  for(int i = 0; i < VGFvectorDim; i++){
    data.at(TDVGFvcoord).set(i,GFcoord,GFacoord,vec[i]);
  }
}

//Get the value of a specific point.
template <class T>
T TwoDVectorGridFunction<T>::get(int TDVGFvcoord, int VGFvcoord, 
                                 int GFvcoord, int GFacoord)
{
  if((0 > TDVGFvcoord) || (TDVGFvcoord >= TDVGFvectorDim)){
    throw invalid_argument("2D vector index out of range in get");
  }
  return data.at(TDVGFvcoord).get(VGFvcoord, GFvcoord, GFacoord);
}

//Get an array from an inner, middle, and outer vector coordinate.
template <class T>
Array1D<T> TwoDVectorGridFunction<T>::get(int TDVGFvcoord, int VGFvcoord,
                                          int GFvcoord)
{
  if((0 > TDVGFvcoord) || (TDVGFvcoord >= TDVGFvectorDim)){
    throw invalid_argument("2D vector index out of range in get");
  }
  return data.at(TDVGFvcoord).get(VGFvcoord, GFvcoord);
}
 

//Get a GridFunction from an outer and middle vector coordinate.
template <class T>
GridFunction<T> TwoDVectorGridFunction<T>::get(int TDVGFvcoord, int VGFvcoord)
{
  if((0 > TDVGFvcoord) || (TDVGFvcoord >= TDVGFvectorDim)){
    throw invalid_argument("2D vector index out of range in get");
  }
  return data.at(TDVGFvcoord).get(VGFvcoord);
}


//Append a VectorGridFunction to the end of the TwoDVectorGridFunction.
template <class T>
void TwoDVectorGridFunction<T>::append(VectorGridFunction<T> vgf)
{
  TDVGFvectorDim++;
  data.push_back(vgf);
}

//Get a vector from an outer vector coordinate and an inner vector 
//coordinate and an array coordinate.
template <class T>              
vector<T> TwoDVectorGridFunction<T>::getVector(int TDVGFvcoord, 
                                               int GFvcoord, int GFacoord)
{
  if((GFvcoord < 0) || (GFvcoord >= GFvectorDim) || (GFacoord < 0)
     || (GFacoord > GFarrayDim) || (TDVGFvcoord < 0) || 
     (TDVGFvcoord >= TDVGFvectorDim)) {
    throw invalid_argument("Get indices out of range in TwoDVectorGridFunction.");
  }
  vector<T> outputvec;
  for(int j = 0; j < TDVGFvectorDim; j++) {
    for(int i = 0; i < VGFvectorDim; i++) {
      outputvec.push_back(get(j, i, GFvcoord, GFacoord));
    }
  }
  return outputvec;
}
    
//HERE
//Get an Array1D from an outer vector coordinate, a middle vector coordinate,
//and inner vector coordinate 
//and an array coordinate, 
//ranging from middle (dimension 0) or outer (dimension 1) vector index 
//vmin to vmax with coordinate in the other dimension vectorCoord.
template <class T>
Array1D<T> TwoDVectorGridFunction<T>::getVectorAsArray1D(int vectorCoord,
                                                         int GFvcoord, 
                                                         int GFacoord, 
                                                         int vmin, 
                                                         int vmax,
                                                         int dimension)
{
  if((GFvcoord < 0) || (GFvcoord >= GFvectorDim) || (GFacoord < 0)
     || (GFacoord > GFarrayDim)) {
    throw invalid_argument("Get indices out of range in TwoDVectorGridFunction");
  }
  
  if((vmax >= VGFvectorDim) || (vmin < 0)){
    throw invalid_argument("Get max or min vector indices out of range in TwoDVectorGridFunction");
  }
  
  Array1D<T> outputvec(vmax - vmin + 1);
  if(dimension==0){
    for(int i = vmin; i <= vmax; i++){
      outputvec[i - vmin] =get(vectorCoord, i, GFvcoord, GFacoord);
    }
  } else if (dimension == 1) {
    for(int j = vmin; j <= vmax; j++) {
      outputvec[j - vmin] = get(j, vectorCoord, GFvcoord, GFacoord);
    }
  } else {
    throw invalid_argument("Dimension out of range in TwoDVectorGridFunction::getVectorAsArray1D.");
  }
  return outputvec;
}

//Get an Array2D containing the values of the VGF at an inner vector
//coordinate, with the middle (dimension 0) or outer (dimension 1)
// vector values running along the rows
//and the array values running along the columns. Begins selecting
//vector values at startvec and ends at stopvec. The unvaried vector dimension
//is selected at vectorCoord.
template <class T>
Array2D<T> TwoDVectorGridFunction<T>::getVectorNodeArray2D(int vectorCoord,
                                                           int GFcoord,
                                                           int startvec, 
                                                           int stopvec,
                                                           int dimension)
{
  if((GFcoord < 0) || (GFcoord > GFvectorDim)) {
    throw invalid_argument("Get indices out of range");
  }
  
  if((startvec < 0) || (stopvec >= VGFvectorDim)) {
    throw invalid_argument("Endpoints of vector requested are out of range");
  }
  
  Array2D<T> output(GFarrayDim, stopvec - startvec + 1);
  for(int i = startvec; i <= stopvec; i++){
    for(int j = 0; j < GFarrayDim; j++){
      if(dimension==0){
        output[j][i - startvec]=get(vectorCoord,i, GFcoord, j);
      } else if (dimension == 1) {
        output[j][i - startvec] = get(i, vectorCoord, GFcoord, j);
      } else {
        throw invalid_argument("Dimension out of range in TwoDVectorGridFunction::getVectorNodeArray2D");
      }
    }
  }
  return output;
}
  
//-----------------------------
//not in class

//Addition operator for VectorGridFunctions.
template <typename T>
TwoDVectorGridFunction<T> operator+(TwoDVectorGridFunction<T> tdvgf1, 
                                TwoDVectorGridFunction<T> tdvgf2)
{
  if(tdvgf1.vectorDim() != tdvgf2.vectorDim()){
    throw invalid_argument("Second dimension vector dimension mismatch in + operation");
  }
  TwoDVectorGridFunction<T> tdvgfsum(0, tdvgf1.vectorDim(), tdvgf1.gridDim(), 
                                 tdvgf1.pointsDim());
  for(int i = 0; i < tdvgf1.modesDim(); i++){
    tdvgfsum.append(tdvgf1.get(i) + tdvgf2.get(i));
  }
  return tdvgfsum;
}

//Scalar multiplication operator for TwoDVectorGridFunctions.
template <typename T>
TwoDVectorGridFunction<T> operator*(T A, TwoDVectorGridFunction<T> tdvgf)
//for easy multiplication in rk4 routine
{
  TwoDVectorGridFunction<T> tdvgfprod(0, tdvgf.vectorDim(), tdvgf.gridDim(), 
                                      tdvgf.pointsDim());
  for(int i = 0; i < tdvgf.modesDim(); i++){
   tdvgfprod.append(A * tdvgf.get(i));
  }
  return tdvgfprod;
}

  

