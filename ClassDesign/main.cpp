#include "tnt_array1d.h"
#include "tnt_array2d."
#include "tnt_array1d_utils.h"
#include "tnt_array2d_utils.h"
#include "DifferentialEquation.h"
#include "TimeEvolution.h"
#include "GridFunction.h"
#include "VectorGridFunction.h"

int main()
{
  int timesteps=10000;
  int PDEnum = 1; //number of independent PDEs. 
  int Np=21; //order of element

  //initialization of reference element
  ReferenceElement refelem();
  refelem::initialize(Np);
  
  //initialization of grid
  Grid thegrid(fileElemBoundaries, refelem);
  int NumElem=theGrid.getN();

  //declaration of calculation variables and 
  //initialization to either zero or value read from file
  VectorGridFunction(PDEnum,N,Np,false) uh; 
  //solution to PDE, possibly a vector 
  VectorGridFunction(PDEnum,N,Np,true) RHS; //right hand side of PDE
  uh.initFromFile(scalarfilename);

  //evolve in time and print out at select time steps
  TimeEvolution te();
  
  for(i=0;i<timesteps;i++){
    TimeEvolution::rk4lowstorage(thegrid.getNodeLocations(), uh, RHS);
    if(outputcondition) uh.print(filenames);
  }

}
