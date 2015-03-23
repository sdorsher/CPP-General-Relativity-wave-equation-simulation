#include "GridFunction.h"

using namespace std;
GridFunction::GridFunction(int vecSize, int arraySize,bool initZeros):GFvectorDim(vecSize),GFarrayDim(arraySize)
{
  if((vecSize<0)||(arraySize<0))
    {
      throw invalid_argument("Negative grid function dimensions at GridFunction constructor.");
    }

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
  if(array.dim()!=GFarrayDim)
    {
     throw invalid_argument("Incorrect array dimensions");
    }

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

TNT::Array1D<double> GridFunction::get(int vcoord)
{
  if((0>vcoord)||(vcoord>=GFvectorDim))
    {
      throw out_of_range("Grid coordinate out of range");
    }
  else
    {
      return data.at(vcoord);
    }
}

int GridFunction::gridDim()
{
  return GFvectorDim;
}

int GridFunction::pointsDim()
{
  return GFarrayDim;
}


void GridFunction::save(string filename)
{

  ofstream fs;
  fs.open(filename);
  for(int i=0;i<GFvectorDim;i++)
    {
      for(int j=0;j<GFarrayDim;j++)
        {
          fs << data.at(i)[j] <<endl;
        }
    }
  fs.close();
}


//-----------------------------------------
// Not in class

GridFunction operator+(GridFunction gf1,GridFunction gf2)
{
  if((gf1.gridDim()!=gf2.gridDim())||(gf1.pointsDim()!=gf2.pointsDim()))
    {
      throw invalid_argument("Grid function dimension mismatch in + operator");
    }
  else
    {
      GridFunction gfout(0,gf1.pointsDim(),false);
      for(int i=0; i<gf1.gridDim(); i++)
        {
          gfout.append(gf1.get(i)+gf2.get(i));
        }
      return gfout;
    }
}
