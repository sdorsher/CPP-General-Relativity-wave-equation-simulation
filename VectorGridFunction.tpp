#include "VectorGridFunction.h"


template <class T>
VectorGridFunction<T>::VectorGridFunction(int outerVecSize, 
                           int innerVecSize, int arraySize, T initvalue):
                             VGFvectorDim(outerVecSize),GFvectorDim(innerVecSize),GFarrayDim(arraySize)
{
  if(outerVecSize<0)
    {
      throw invalid_argument("Negative Vector dimensions at VectorGridFunction constructor.");
    }
  
  for(int i=0; i<outerVecSize; i++)
    {
      GridFunction<T> temp(innerVecSize,arraySize,initvalue);
      data.push_back(temp);
    }
  
}
template <class T>
VectorGridFunction<T>::VectorGridFunction(int outerVecSize, 
                           int innerVecSize, int arraySize):
                             VGFvectorDim(outerVecSize),GFvectorDim(innerVecSize),GFarrayDim(arraySize)
{
  if(outerVecSize<0)
    {
      throw invalid_argument("Negative Vector dimensions at VectorGridFunction constructor.");
    }
  
  for(int i=0; i<outerVecSize; i++)
    {
      GridFunction<T> temp(innerVecSize,arraySize);
      data.push_back(temp);
    }
  
}

template <class T>
int VectorGridFunction<T>::vectorDim()
{
  return VGFvectorDim;
}
template <class T>
int VectorGridFunction<T>::gridDim()
{
  return GFvectorDim;
}
template <class T>
int VectorGridFunction<T>::pointsDim()
{
  return GFarrayDim;
}

template <class T>
void VectorGridFunction<T>::set(int VGFvcoord, int GFvcoord, int GFacoord,T value)
{
  if((0>VGFvcoord)||(VGFvcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in set");
    }
  
  data.at(VGFvcoord).set(GFvcoord,GFacoord,value);
}
template <class T>
void VectorGridFunction<T>::set(int VGFcoord,GridFunction<T> gf)
{

  if((0>VGFcoord)||(VGFcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in set");
    }

  data.at(VGFcoord)=gf;
}

template <class T>
void VectorGridFunction<T>::set(int VGFcoord,int GFcoord,TNT::Array1D<T> arr)
{

  if((0>VGFcoord)||(VGFcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in set");
    }

  data.at(VGFcoord).set(GFcoord,arr);

}
template <class T>
void VectorGridFunction<T>::setVector(int GFcoord, int GFacoord, vector<T> vec)
{
   if((GFcoord<0)||(GFacoord<0)||(GFacoord>GFarrayDim))
    {
      throw invalid_argument("Coordinates in setVector out of range.");
    }
  for(int i=0; i<VGFvectorDim; i++)
    {
      data.at(i).set(GFcoord,GFacoord,vec[i]);
    }
    
    
}

template <class T>
T VectorGridFunction<T>::get(int VGFvcoord, int GFvcoord, int GFacoord)
{
  if((0>VGFvcoord)||(VGFvcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in get");
    }
  return data.at(VGFvcoord).get(GFvcoord,GFacoord);
}
template <class T>
Array1D<T> VectorGridFunction<T>::get(int VGFvcoord,int GFvcoord)
{
  if((0>VGFvcoord)||(VGFvcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in get");
    }
  return data.at(VGFvcoord).get(GFvcoord);
}
 
template <class T>
GridFunction<T> VectorGridFunction<T>::get(int VGFvcoord)
{
  if((0>VGFvcoord)||(VGFvcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in get");
    }
  return data.at(VGFvcoord);

}
template <class T>
void VectorGridFunction<T>::append(GridFunction<T> gf)
{
  VGFvectorDim++;
  data.push_back(gf);
}

template <class T>
void VectorGridFunction<T>::save(string filename)
{
  ofstream fs;
  fs.open(filename);
   
  for(int j=0; j<GFvectorDim; j++)
    {
      for(int k = 0; k<GFarrayDim; k++)
        {
          for(int i=0; i<VGFvectorDim; i++)
            {
              fs << data.at(i).get(j,k) << " ";
            }
          fs << endl;
        }
      
    }
  fs.close();
}

template <class T>              
vector<T> VectorGridFunction<T>::getVector(int GFvcoord, int GFacoord)
{
  if((GFvcoord<0)||(GFvcoord>=GFvectorDim)||(GFacoord<0)||(GFacoord>GFarrayDim))
    {
      throw invalid_argument("Get indices out of range.");
    }
  
  vector<T> outputvec;
  for(int i=0; i<VGFvectorDim; i++)
    {
      outputvec.push_back(get(i,GFvcoord,GFacoord));
    }
  return outputvec;

}
    

template <class T>
Array1D<T> VectorGridFunction<T>::getVectorAsArray1D(int GFvcoord, int GFacoord,int vmin, int vmax)

{
  if((GFvcoord<0)||(GFvcoord>=GFvectorDim)||(GFacoord<0)||(GFacoord>GFarrayDim))
    {
      throw invalid_argument("Get indices out of range.");
    }
  
  if((vmax>VGFvectorDim)||(vmin<0))
    {
      throw invalid_argument("Get max or min vector indices out of range");
    }

  Array1D<T> outputvec(vmax-vmin+1);
  for(int i=vmin; i<=vmax; i++)
    {
      outputvec[i-vmin] =get(i,GFvcoord,GFacoord);
    }
  return outputvec;

}

template <class T>
Array2D<T> VectorGridFunction<T>::getVectorNodeArray2D(int GFcoord,int startvec, int stopvec)
{
  if((GFcoord<0)||(GFcoord>GFvectorDim))
    {
      throw invalid_argument("Get indices out of range");
    }
  if((startvec<0)||(stopvec>=VGFvectorDim))
    {
      throw invalid_argument("Endpoints of vector requested are out of range");
    }

  Array2D<T> output(GFarrayDim,stopvec-startvec+1);
  for(int i=startvec; i<=stopvec; i++)
    {
      for(int j=0; j<GFarrayDim; j++)
        {
          output[j][i-startvec]=get(i,GFcoord,j);
        }
    }
  return output;
    
}



//-----------------------------
//not in class

template <typename T>
VectorGridFunction<T> operator+(VectorGridFunction<T> vgf1, VectorGridFunction<T> vgf2)
{
  if(vgf1.vectorDim()!=vgf2.vectorDim())
    {
      throw invalid_argument("Vector dimension mismatch in + operation");
    }
    VectorGridFunction<T> vgfsum(0,vgf1.gridDim(),vgf1.pointsDim());
    for(int i=0; i<vgf1.vectorDim(); i++)
      {
        vgfsum.append(vgf1.get(i)+vgf2.get(i));
      }
    return vgfsum;
    
}

template <typename T>
VectorGridFunction<T> operator*(T A, VectorGridFunction<T> vgf)
//for easy multiplication in rk4 routine
{
  VectorGridFunction<T> vgfprod(0,vgf.gridDim(),vgf.pointsDim());
  for(int i=0; i<vgf.vectorDim(); i++)
    {
      vgfprod.append(A*vgf.get(i));
    }
  return vgfprod;
}

  

