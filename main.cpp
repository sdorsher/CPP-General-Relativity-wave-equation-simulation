#include "tnt.h"
#include "TNT2.h"
#include "Grid.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include "VectorGridFunction.h"
#include "Evolution.h"
#include <fstream>

double analyticsoln(double);
void initialconditions(VectorGridFunction& uh, Grid grd);
void Linferror(double nominal,double theoretical, double&);

int main()
{
  int PDEnum = 2; //number of independent PDEs. 
  int Np=10; //order of element
  int NumElems = 5;
  string fileElemBoundaries= "fiveElems.txt";
  //vector<string> uhfilename
  // string uhfilename="fourthorderODE.txt";
  string uhfilename="waveequation.txt";
  
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

  initialconditions(uh,thegrid);

  //print out at select time steps
  
  ofstream fs;
  fs.open(uhfilename);

  double t0=0.0;
  double tmax=10.0;
  double deltat=0.01;
  for(double t=t0; t<tmax; t+=deltat)
    {
      rk4lowStorage(thegrid,uh,RHSvgf,t,deltat);
    }
  
  //initial conditions, numerical fluxes, boundary conditions handled inside 
  //Evolution.cpp, in RHS.


  
}

/*
//analytic solution to du=t^4 with t=0 as initial condition
double analyticsoln(double t)
{
return 0.2*pow(t,5.0);
}
*/
void initialconditions(VectorGridFunction& uh, Grid grd){
   double sigma = 0.1;
   double amplitude = 1.0;
   double position = 5.0;
   GridFunction nodes(uh.gridDim(),uh.pointsDim(),false);
   nodes=grd.gridNodeLocations();
   for(int i=0; i<uh.gridDim(); i++)
     {
       for(int j=0; j<uh.pointsDim(); j++)
         {
           double gauss=amplitude*exp(-pow((nodes.get(i,j)-position),2.0)
                                      /2.0/pow(sigma,2.0));
           uh.set(0,i,j,gauss);
           uh.set(1,i,j,0.0);
         }
     }
}

void Linferror(double nominal,double theoretical,double& maxerror)
 {
   
   double newerror= fabs(nominal-theoretical);
   maxerror = newerror>maxerror ? newerror : maxerror;
 }

 //analytic solution to du=-omega^2 u with initial conditions of A=2.0
 //and omega = 1.0

 /* double analyticsoln(double t)
    {
    double A=2.0;
    double omega=1.0;
    return A*pow(omega,2.0)*cos(omega*t);
    }*/

//initial conditions for SHO
 /*void initialconditions(VectorGridFunction& uh) {
   vector<double> A{2.0,0.0};
   for(int i=0; i<A.size();i++){
    for(int j=0; j<uh.gridDim(); j++)
    {
    for(int k=0; k<uh.pointsDim(); k++)
    {
    uh.set(i,j,k,A[i]);
    }
    }
  }
  
  }*/


