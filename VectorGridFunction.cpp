#include "VectorGridFunction.h"



VectorGridFunction::VectorGridFunction(int outerVecSize, int innerVecSize, int arraySize, bool initZero):VGFvectorDim(outerVecSize),GFvectorDim(innerVecSize),GFarrayDim(arraySize)
{
  if(outerVecSize<0)
    {
      throw invalid_argument("Negative Vector dimensions at VectorGridFunction constructor.");
    }
  
  for(int i=0; i<outerVecSize; i++)
    {
      GridFunction temp(innerVecSize,arraySize,initZero);
      data.push_back(temp);
    }
  
}

int VectorGridFunction::vectorDim()
{
  return VGFvectorDim;
}

int VectorGridFunction::gridDim()
{
  return GFvectorDim;
}

int VectorGridFunction::pointsDim()
{
  return GFarrayDim;
}


void VectorGridFunction::set(int VGFvcoord, int GFvcoord, int GFacoord,double value)
{
  if((0>VGFvcoord)||(VGFvcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in set");
    }
  
  data.at(VGFvcoord).set(GFvcoord,GFacoord,value);
}

void VectorGridFunction::set(int VGFcoord,GridFunction gf)
{

  if((0>VGFcoord)||(VGFcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in set");
    }

  data.at(VGFcoord)=gf;
}


void VectorGridFunction::set(int VGFcoord,int GFcoord,TNT::Array1D<double> arr)
{

  if((0>VGFcoord)||(VGFcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in set");
    }

  data.at(VGFcoord).set(GFcoord,arr);

}


double VectorGridFunction::get(int VGFvcoord, int GFvcoord, int GFacoord)
{
  if((0>VGFvcoord)||(VGFvcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in get");
    }
  return data.at(VGFvcoord).get(GFvcoord,GFacoord);
}

Array1D<double> VectorGridFunction::get(int VGFvcoord,int GFvcoord)
{
  if((0>VGFvcoord)||(VGFvcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in get");
    }
  return data.at(VGFvcoord).get(GFvcoord);
}
 

GridFunction VectorGridFunction::get(int VGFvcoord)
{
  if((0>VGFvcoord)||(VGFvcoord>=VGFvectorDim))
    {
      throw invalid_argument("Vector index out of range in get");
    }
  return data.at(VGFvcoord);

}

void VectorGridFunction::append(GridFunction gf)
{
  VGFvectorDim++;
  data.push_back(gf);
}


void VectorGridFunction::save(string filename)
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

              
vector<double> VectorGridFunction::getVector(int GFvcoord, int GFacoord)
{
  if((GFvcoord<0)||(GFvcoord>=GFvectorDim)||(GFacoord<0)||(GFacoord>GFarrayDim))
    {
      throw invalid_argument("Get indices out of range.");
    }
  
  vector<double> outputvec;
  for(int i=0; i<VGFvectorDim; i++)
    {
      outputvec.push_back(get(i,GFvcoord,GFacoord));
    }
  return outputvec;

}
    


//-----------------------------
//not in class

VectorGridFunction operator+(VectorGridFunction vgf1, VectorGridFunction vgf2)
{
  if(vgf1.vectorDim()!=vgf2.vectorDim())
    {
      throw invalid_argument("Vector dimension mismatch in + operation");
    }
    VectorGridFunction vgfsum(0,vgf1.gridDim(),vgf1.pointsDim(),false);
    for(int i=0; i<vgf1.vectorDim(); i++)
      {
        vgfsum.append(vgf1.get(i)+vgf2.get(i));
      }
    return vgfsum;
    
}


//looper function
//all functions need to inherrit from function class 
//(or this needs to be a template function)



void loop(GridFunction& grid, VectorGridFunction& input, VectorGridFunction& output, FUNCTYPE func)
{
  //func should be function not class

  if (input.vectorDim()!=output.vectorDim())
    {
      throw invalid_argument("uh and RHS dimensions do not agree in loop.");
    }


  for(int vecindex=0; vecindex<output.vectorDim(); vecindex++)
    {
      for(int gridindex =0; gridindex<output.gridDim(); gridindex++)
        {
          (func)(grid, input, output, vecindex, gridindex);
          //is this a general enough functional form?
        }
    }
}

void testfunc(GridFunction& grid, VectorGridFunction& uh, VectorGridFunction& RHS, int vecindex, int gridindex)
{
  for(int arrayindex=0; arrayindex<RHS.pointsDim(); arrayindex++)
    {
      RHS.set(vecindex,gridindex,arrayindex,grid.get(gridindex,arrayindex)*vecindex);
    }
}
