#include "tnt.h"
#include "TNT2.h"
#include "Grid.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include "VectorGridFunction.h"
#include "Evolution.h"
#include <fstream>

double analyticsoln(double);

int main()
{
  int PDEnum = 1; //number of independent PDEs. 
  int Np=3; //order of element
  int NumElems = 1;
  string fileElemBoundaries= "oneElem.txt";
  //vector<string> uhfilename
  string uhfilename="fourthorderODE.txt";

  //initialization of grid and calculation of reference element
  Grid thegrid(fileElemBoundaries, Np, NumElems);

  //declaration of calculation variables and 
  //initialization to either zero or value read from file
  VectorGridFunction uh(PDEnum,NumElems,Np,true); 
  //solution to PDE, possibly a vector 
  VectorGridFunction RHSvgf(PDEnum,NumElems,Np,true); //right hand side of PDE
  
  GridFunction nodes = thegrid.gridNodeLocations();

  //initial conditions
  //uh.initFromFile(scalarfilename);

  //print out at select time steps
  
  ofstream fs;
  fs.open(uhfilename);

  double t0=0.0;
  double deltat=0.01;
  double tmax=10.0;
  for(double t=t0;t<tmax;t+=deltat){
    fs << t << " " << uh.get(0,0,0) <<" " << analyticsoln(t) <<endl;
    rk4lowStorage(nodes,uh,RHSvgf,t,deltat);

    //need to do something about numerical fluxes
    //    if(outputcondition) uh.save(uhfilename);
  }
  fs.close();

}

double analyticsoln(double t)
{
  return 0.2*pow(t,5.0);
}
