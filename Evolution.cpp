#include "Evolution.h"

//Fourth order low storage Runga Kutta method for time integration

//coefficients copied and pasted from Fortran code by Peter Diener.
//He thinks the coefficients were copied and pasted from the Matlab 
//code by Hesthaven and Warburn.
//See page 63 of Hesthaven and Warburn for this routine

void rk4lowStorage(Grid thegrid, DiffEq theequation, 
                   TwoDVectorGridFunction<complex<double>>& uh, 
                   TwoDVectorGridFunction<complex<double>>& RHStdvgf, 
                   double t, double deltat)
{

  /*
  vector<double> rk4a{0.0, 
      -567301805773.0/1357537059087.0,
      -2404267990393.0/2016746695238.0,
      -3550918686646.0/2091501179385.0,
      -1275806237668.0/842570457699.0};
  vector<double> rk4b{1432997174477.0/9575080441755.0,
      5161836677717.0/13612068292357.0,
      1720146321549.0/2090206949498.0,
      3134564353537.0/4481467310338.0,
      2277821191437.0/14882151754819.0};
  vector<double> rk4c{0.0,
      1432997174477.0/9575080441755.0, 
      2526269341429.0/6820363962896.0, 
      2006345519317.0/3224310063776.0, 
      2802321613138.0/2924317926251.0};
  */
  
  vector<double> rk4a;
  rk4a[0]= 0.0; 
  rk4a[1]=-567301805773.0/1357537059087.0;
  rk4a[2]=-2404267990393.0/2016746695238.0;
  rk4a[3]=-3550918686646.0/2091501179385.0;
  rk4a[4]=-1275806237668.0/842570457699.0;

  vector<double> rk4b;
  rk4b[0]=1432997174477.0/9575080441755.0;
  rk4b[0]=5161836677717.0/13612068292357.0;
  rk4b[0]=1720146321549.0/2090206949498.0;
  rk4b[0]=3134564353537.0/4481467310338.0;
  rk4b[0]=2277821191437.0/14882151754819.0;

  vector<double> rk4c;
  rk4c[0]=0.0;
  rk4c[1]=1432997174477.0/9575080441755.0; 
  rk4c[2]=2526269341429.0/6820363962896.0; 
  rk4c[3]=2006345519317.0/3224310063776.0; 
  rk4c[4]=2802321613138.0/2924317926251.0;

  int nsteps=5;
  
  TwoDVectorGridFunction<complex<double>> k(RHStdvgf.modesDim(),
                                   RHStdvgf.vectorDim(), RHStdvgf.gridDim(),
                                   RHStdvgf.pointsDim());

  //step 0
  //  vector<Array2D<double>> du;


  //step 1
  theequation.modeRHS(thegrid, uh, RHStdvgf, t,false);
  k = deltat * RHStdvgf;
  uh = uh + rk4b[0] * k;


  //steps 2-5
  for(int i=2; i<=5; i++){
    theequation.modeRHS(thegrid,uh, RHStdvgf, t + rk4c[i-1] * deltat,false);
    k = rk4a[i-1] * k + deltat * RHStdvgf;
    uh = uh + rk4b[i-1] * k;
  }

  /*  if(RHSOUTPUT){
    ofstream fs;
    fs.open("urhs.txt", ios::app);
    fs << endl << endl;
    fs << " #time = " << t << endl;
    for(int j=0; j< RHStdvgf.gridDim(); j++){
      for(int i=0; i< RHStdvgf.pointsDim(); i++){
        fs << thegrid.gridNodeLocations().get(j,i) << " ";
        for(int k = 0; k< RHStdvgf.vectorDim(); k++){
          fs << RHStdvgf.get(0,k,j,i) << " ";
        }
        fs << endl;
      }
    }
    fs.close();
    }*/
}
