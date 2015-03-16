#include "GridFunction.h"

using namespace std;
GridFunction::GridFunction(int vecSize, int arraySize,bool initZeros):GFvectorDim(vecSize),GFarrayDim(arraySize)
{
  for(int v=0; v<GFvectorDim; v++)
    {
      if(initZeros)
        {
          TNT::Array1D<double> temp(GFarrayDim,0.0);
          data.push_back(temp);
        }
      else
        {
          TNT::Array1D<double> temp(GFarrayDim);
          data.push_back(temp);
        }
    }
}

void GridFunction::set(int vcoord, TNT::Array1D<double> arraydata)
{
  if((0>vcoord)||(vcoord>=GFvectorDim))
    {
      throw out_of_range("Grid coordinate out of range");
    }
  else if (arraydata.dim()!=GFarrayDim)
    {
      throw invalid_argument("Grid function data size does not match.");
    }
  else
    {
      data.at(vcoord)=arraydata.copy();
    }
}

void GridFunction::set(int vcoord, int acoord,double value)
{
  if((0>vcoord)||(vcoord>=GFvectorDim))
    {
      throw out_of_range("Grid coordinate out of range");
    }
  else if((0>acoord)||(acoord>=GFarrayDim))
    {
      throw out_of_range("Grid function coordinate out of range");
    }
  else
    {
      data.at(vcoord)[acoord]=value;
    }
}

void GridFunction::append(TNT::Array1D<double> array)
{
  data.push_back(array.copy());
  GFvectorDim++;
}

double GridFunction::get(int vcoord, int acoord)
{
  if((0>vcoord)||(vcoord>=GFvectorDim))
    {
      throw out_of_range("Grid coordinate out of range");
    }
  else if((0>acoord)||(acoord>=GFarrayDim))
    {
      throw out_of_range("Grid function coordinate out of range");
    }
  else
    {
      return data.at(vcoord)[acoord];
    }


}

int GridFunction::gridDim()
{
  return GFvectorDim;
}

int GridFunction::functionDim()
{
  return GFarrayDim;
}

//void GridFunction::operator+(GridFunction gf1, GridFunction gf2)

  
  

