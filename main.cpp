#include "tnt.h"
#include "TNT2.h"
#include "Grid.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include "VectorGridFunction.h"
#include "Evolution.h"
#include <fstream>

double analyticsoln(double);
void initialconditions(VectorGridFunction& uh);
double Linferror(double nominal,double theoretical);

int main()
{
  int PDEnum = 2; //number of independent PDEs. 
  int Np=3; //order of element
  int NumElems = 1;
  string fileElemBoundaries= "oneElem.txt";
  //vector<string> uhfilename
  // string uhfilename="fourthorderODE.txt";
  string uhfilename="shoODE.txt";
  
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
  double tmax=10.0;
  double deltatinit=1.0;
  for(int i=0; i<7; i++)
    {
      cout << i<< endl;
      initialconditions(uh);
      double deltat=deltatinit*pow(0.1,i);
  
      for(double t=t0;t<tmax;t+=deltat){
        rk4lowStorage(nodes,uh,RHSvgf,t,deltat);
        
        //need to do something about numerical fluxes
        //    if(outputcondition) uh.save(uhfilename);
      }
      double nominal=uh.get(0,0,0);
      double theoretical=analyticsoln(tmax);
      fs << deltat << " " <<nominal << " " <<theoretical << " " 
         << Linferror(nominal, theoretical) << endl;
    }
  fs.close();

}

/*
//analytic solution to du=t^4 with t=0 as initial condition
double analyticsoln(double t)
{
  return 0.2*pow(t,5.0);
}
*/

 //analytic solution to du=-omega^2 u with initial conditions of A=2.0
 //and omega = 1.0

 double analyticsoln(double t)
{
  double A=2.0;
  double omega=1.0;
  return A*pow(omega,2.0)*cos(omega*t);
}

void initialconditions(VectorGridFunction& uh) {
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

}

double Linferror(double nominal,double theoretical)
 {
   return fabs(nominal-theoretical);
 }
