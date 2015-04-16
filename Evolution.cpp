#include "Evolution.h"

void rk4lowStorage(GridFunction& nodes, VectorGridFunction& uh, 
                   VectorGridFunction& RHSvgf, 
                   double t, double deltat)
{


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

  int nsteps=5;
  
  VectorGridFunction k(RHSvgf.vectorDim(),RHSvgf.gridDim(),RHSvgf.pointsDim(),false);

  //step 0

  
  //step 1
  RHS(nodes, uh, RHSvgf, t);
  k=deltat*RHSvgf;
  uh=uh+rk4b[0]*k;

  //step2
  for(int i=2; i<=5; i++){
    RHS(nodes, uh, RHSvgf, t+rk4c[i-1]*deltat);
    k=rk4a[i-1]*k+deltat*RHSvgf;
    uh=uh+rk4b[i-1]*k;
  }

}

void RHS(GridFunction& nodes, VectorGridFunction& uh, 
             VectorGridFunction& RHSvgf, double t)
{

  for(int vecindex=0;vecindex<RHSvgf.vectorDim(); vecindex++){
    for(int gridindex=0; gridindex<RHSvgf.gridDim(); gridindex++)
      {
        for(int pointsindex=0;pointsindex<RHSvgf.pointsDim(); pointsindex++){
          RHSvgf.set(vecindex,gridindex,pointsindex,pow(t,4.0));
        }
      }
  }


}
