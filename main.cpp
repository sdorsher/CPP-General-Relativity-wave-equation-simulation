#include "tnt.h"
#include "TNT2.h"
#include "Grid.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include "VectorGridFunction.h"
#include "Evolution.h"
#include "globals.h"
#include <cmath>
#include <fstream>

double analyticsoln(double);
void initialconditions(VectorGridFunction& uh, Grid grd);
void Linferror(double nominal,double theoretical, double&);

int main()
{
  int PDEnum = 3; //number of independent PDEs. 
  int NumElems=20;
  string fileElemBoundaries= "elemBoundaries.txt";
  //vector<string> uhfilename
  // string uhfilename="fourthorderODE.txt";
  string uhfilename="waveequation.txt";
  
  //initialization of grid and calculation of reference element
  Grid thegrid(fileElemBoundaries, ELEMORDER, NumElems);

  //declaration of calculation variables and 
  //initialization to either zero or value read from file
  VectorGridFunction uh(PDEnum,NumElems,ELEMORDER+1,true); 
  VectorGridFunction uh0(PDEnum,NumElems,ELEMORDER+1,true);
  //solution to PDE, possibly a vector 
  VectorGridFunction RHSvgf(PDEnum,NumElems,ELEMORDER+1,true); //right hand side of PDE
 


  GridFunction nodes = thegrid.gridNodeLocations();

  //initial conditions
  //uh.initFromFile(scalarfilename);

  //get problem in initial conditions
  initialconditions(uh,thegrid);
  uh0=uh;
  

  //print out at select time steps
  
  double dt0=nodes.get(0,1)-nodes.get(0,0);
  //sets to smallest grid spacing


  ofstream fs;
  fs.open(uhfilename);

  double t0=0.0;
  double tmax=20.0;
  // double deltat=0.002;
  
  int nt=ceil(tmax/dt0);
  double deltat=tmax/nt;

  cout << dt0 << " " << deltat << endl;

  double outputinterval=1.0;
  double output=deltat/2.0;
  int outputcount =0;

  for(double t=t0; t<tmax+deltat; t+=deltat)
    {
      
      if(output>0.0){
        fs<<endl <<endl;
        fs<< " #time = " << t << endl;
        for (int i=0; i<uh.gridDim(); i++)
          {
            for(int j=0; j<uh.pointsDim(); j++)
              {
                fs << thegrid.gridNodeLocations().get(i,j) << " " << uh.get(0,i,j) << endl;
              }
            
          }
        if(outputcount==20)
          {
            ofstream fs2;
            fs2.open("diffGaus.txt");
            for(int i=0; i<uh.gridDim(); i++)
              {
                for (int j=0; j<uh.pointsDim(); j++)
                  {

                    fs2 << thegrid.gridNodeLocations().get(i,j) << " " << uh.get(0,i,j) - uh0.get(0,i,j) <<endl;
                  }}
            fs2.close();}
        output-=outputinterval; 
        outputcount++;
      }
   
      rk4lowStorage(thegrid,uh,RHSvgf,t,deltat);
      output+=deltat;
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
   double sigma = 1.0;
   double amplitude = 1.0;
   double position = 10.1;
   GridFunction nodes(uh.gridDim(),uh.pointsDim(),false);
   nodes=grd.gridNodeLocations();

   // std::ofstream fs;
   //fs.open("initialconditions.txt");
   

   for(int i=0; i<uh.gridDim(); i++)
     {
       for(int j=0; j<uh.pointsDim(); j++)
         {
           double gauss=amplitude*exp(-pow((nodes.get(i,j)-position),2.0)
                                      /2.0/pow(sigma,2.0));
           double dgauss =-(nodes.get(i,j)-position)/pow(sigma,2.0)*gauss;
           uh.set(0,i,j,gauss);
           uh.set(1,i,j,0.0);
           //uh.set(1,i,j,-speed*dgauss);
           uh.set(2,i,j,dgauss);
           // fs << nodes.get(i,j) << " " << gauss << " " << speed*dgauss << " " << -dgauss << endl;
         }
     }
   //   fs.close();

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


