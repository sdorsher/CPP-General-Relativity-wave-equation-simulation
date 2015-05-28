#include "TNT2.h"
#include "Grid.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include "VectorGridFunction.h"
#include "Evolution.h"
#include "globals.h"
#include <cmath>
#include <fstream>
#include "ConfigParams.h"


double analyticsoln(double);
void initialGaussian(VectorGridFunction<double>& uh, Grid grd);
void initialSinusoid(VectorGridFunction<double>& uh, Grid grd);
void Linferror(double nominal,double theoretical, double&);
double LTwoerror(Grid thegrid, VectorGridFunction<double>& uh0, 
                 VectorGridFunction<double>& uhend);

int main()
{
  //string filename="elemBoundaries.txt";

  Grid thegrid(params.grid.elemorder,params.grid.numelems,params.grid.lowerlim, params.grid.upperlim);

  //declaration of calculation variables and 
  //initialization to either zero or value read from file
  VectorGridFunction<double> uh(params.waveeq.pdenum,
                                params.grid.numelems,
                                params.grid.elemorder+1,
                                0.0); 
  VectorGridFunction<double> uh0(params.waveeq.pdenum,
                                 params.grid.numelems,
                                 params.grid.elemorder+1,
                                 0.0);
  //solution to PDE, possibly a vector 
  VectorGridFunction<double> RHSvgf(params.waveeq.pdenum,
                                    params.grid.numelems,
                                    params.grid.elemorder+1,
                                    0.0); //right hand side of PDE
 



  GridFunction<double> nodes = thegrid.gridNodeLocations();



  //initial conditions
  //uh.initFromFile(scalarfilename);

  //get problem in initial conditions

  if(params.waveeq.issinusoid)
    {
      initialSinusoid(uh,thegrid);
    }else if(params.waveeq.isgaussian)
    {
      initialGaussian(uh,thegrid);
    }
  uh0=uh;
  
  //  double integral= LTwoerror(thegrid,uh0,uh);
  //cout << integral << endl;
  
  //print out at select time steps
  
  double dt0=nodes.get(0,1)-nodes.get(0,0);
  //sets to smallest grid spacing


  ofstream fs;
  fs.open(params.file.pdesolution);

  //  double t0=0.0;
  //double tmax=20.0;
  //double deltat=0.001;
  
  //  double courantfac=0.25;
  int nt=ceil(params.time.tmax/params.time.courantfac/dt0);

  double deltat;
  if(params.time.usefixedtimestep)
    {
      deltat=params.time.dt;
    }
  else
    {
      deltat=params.time.tmax/nt;
    }
  cout << dt0 << " " << deltat << endl;

  //double outputinterval=1.0;
  double output=deltat/2.0;
  int outputcount =0;

  for(double t=params.time.t0; t<params.time.tmax+deltat; t+=deltat)
    {
      
      if(output>0.0){
        fs<<endl <<endl;
        fs<< " #time = " << t << endl;
        for (int i=0; i<uh.gridDim(); i++)
          {
            for(int j=0; j<uh.pointsDim(); j++)
              {
                fs << thegrid.gridNodeLocations().get(i,j) << " " << uh.get(0,i,j) << " " << uh.get(1,i,j) <<" " << uh.get(2,i,j)<< endl;
              }
            
          }
        if(outputcount==params.time.comparisoncount)
          {
            cout << t << endl;
            ofstream fs2;
            fs2.open(params.file.oneperioderror);
            for(int i=0; i<uh.gridDim(); i++)
              {
                for (int j=0; j<uh.pointsDim(); j++)
                  {

                    fs2 << thegrid.gridNodeLocations().get(i,j) << " " << uh.get(0,i,j) - uh0.get(0,i,j) <<endl;
                  }}
            fs2.close();

            ofstream fsconvergence;
            fsconvergence.open(params.file.L2error,ios::app);
            

            double L2;
            L2=LTwoerror(thegrid,uh0,uh);
            cout << "Order, deltat, num elems, L2 norm" << endl;
            cout<<params.grid.elemorder << " " << deltat << " " << params.grid.numelems << " " << L2 << endl;
            fsconvergence<<params.grid.elemorder << " " << deltat << " " << params.grid.numelems << " " << L2 << endl;
            fsconvergence.close();



}


        output-=params.time.outputinterval; 
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

void initialSinusoid(VectorGridFunction<double>& uh, Grid grd){
    double omega=2.0*PI/params.sine.wavelength;

    GridFunction<double> nodes(uh.gridDim(),uh.pointsDim(),false);
  nodes=grd.gridNodeLocations();

  for(int i=0; i<uh.gridDim(); i++)
    {
      for (int j=0; j<uh.pointsDim(); j++)
        {
          double psi=params.sine.amp*sin(omega*nodes.get(i,j)+params.sine.phase);
          double rho=omega*params.sine.amp*cos(omega*nodes.get(i,j)+params.sine.phase);
          double pivar=-params.waveeq.speed*rho;
          uh.set(0,i,j,psi);
          uh.set(1,i,j,pivar);
          uh.set(2,i,j,rho);
        }
    }
}
 


void initialGaussian(VectorGridFunction<double>& uh, Grid grd){
   // double sigma = 1.0;
   //double amplitude = 1.0;
   //double position = 10.1;
  GridFunction<double> nodes(uh.gridDim(),uh.pointsDim(),false);
   nodes=grd.gridNodeLocations();

   // std::ofstream fs;
   //fs.open(params.file.initialconditions);
   
   for(int i=0; i<uh.gridDim(); i++)
     {
       for(int j=0; j<uh.pointsDim(); j++)
         {
           double gaussian=params.gauss.amp*exp(-pow((nodes.get(i,j)
                                                      -params.gauss.mu),2.0)
                                      /2.0/pow(params.gauss.sigma,2.0));
           double dgauss =-(nodes.get(i,j)-params.gauss.mu)
             /pow(params.gauss.sigma,2.0)*gaussian;
           uh.set(0,i,j,gaussian);
           uh.set(1,i,j,0.0);
           //uh.set(1,i,j,-params.waveeq.speed*dgauss);
           uh.set(2,i,j,dgauss);
           //fs << nodes.get(i,j) << " " << gaussian << " " 
           //<< params.waveeq.speed*dgauss << " " << -dgauss << endl;
         }
     }
   //   fs.close();

}

/*void Linferror(double nominal,double theoretical,double& maxerror)
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

double LTwoerror(Grid thegrid, VectorGridFunction<double>& uh0, VectorGridFunction<double>& uhend)
{
  double L2;
  L2=0.0;
  Array1D<double> weights;
  weights=thegrid.refelem.getw();
  // vector<double> rx;
  //rx=thegrid.jacobian();
  GridFunction<double> nodes(uh0.gridDim(), uh0.pointsDim(),false);
  nodes=thegrid.gridNodeLocations();
  for(int i =0; i<uh0.gridDim(); i++)
    {
      for(int j=0; j<uh0.pointsDim(); j++)
        {
     
          //          cout << nodes.get(i,j) << " " << weights[j] << " " << rx[i] << endl;
          double added= weights[j]*pow(uh0.get(0,i,j)
                                       -uhend.get(0,i,j),2.0)
            /thegrid.jacobian(i);
          L2+=added;
          //cout << L2 <<" "<< added << endl; 
            //<< weights[j]*pow(uh0.get(0,i,j)-uhend.get(0,i,j),2.0)/rx[i] <<endl;
        }

    }
  L2=sqrt(L2); 
  return L2;
}

