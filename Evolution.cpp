#include "Evolution.h"

void rk4lowStorage(Grid thegrid, VectorGridFunction& uh, 
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
  RHS(thegrid, uh, RHSvgf, t);
  k=deltat*RHSvgf;
  uh=uh+rk4b[0]*k;

  //step2
  for(int i=2; i<=5; i++){
    RHS(thegrid, uh, RHSvgf, t+rk4c[i-1]*deltat);
    k=rk4a[i-1]*k+deltat*RHSvgf;
    uh=uh+rk4b[i-1]*k;
  }

}

//---------------------

void RHS(Grid thegrid, VectorGridFunction& uh, 
         VectorGridFunction& RHSvgf, double t)
{
  if(RHSvgf.vectorDim()>2)
    {
      throw invalid_argument("Wave equation PDE should only have two separable components.");
    }

  //wave equation specific values

  //drho/dt - -c^2dPi/dx =0
  //dpi/dt - drho/dx =0

  // PI=[[0 -c^2],[-1 0]]
  //LAMBDA=S^-1 PI S=[[-c 0],[0 c]]
  // S=[v1 v2] so that eigenvalues are diagonal
  // eigenvalues are +-c
  // S=[[c -c],[1 1]]
  // 

  double speed = 1.0;
  //Array2D<double> s(2,2);
  Array2D<double> sinv(2,2);
  Array2D<double> lambda(2,2,0.0);
  // Array2D<double> A(2,2);

  
  /*
  s[0][0]=speed;
  s[0][1]=-speed;
  s[1][0]=1.0;
  s[1][1]=1.0;
  */
  sinv[0][0]=0.5/speed;
  sinv[0][1]=0.5;
  sinv[1][0]=-0.5/speed;
  sinv[1][1]=0.5;

  lambda[0][0]=-speed;
  lambda[1][1]=speed;
  
  vector<double> jacobian=thegrid.jacobian();


  //numerical flux
  
  //index
  vector<int> ind{0,RHSvgf.pointsDim()};

  for(int elemnum=0; elemnum<uh.gridDim(); elemnum++)
    {
  
      //interior,exterior for two variables of wave equation
      Array1D<double> uint1(2,0.0);
      Array1D<double> uint0(2,0.0);
      Array1D<double> uext1(2,0.0);
      Array1D<double> uext0(2,0.0);
  
      for(int i=0; i<2; i++)
        {
          uint0[i] = uh.get(0,elemnum,ind[i]);
          uint1[i]= uh.get(1,elemnum,ind[i]);
        }
      
      Array1D<double> nx(2,1.0);
      nx[0]=-1.0;
      
      if(elemnum>1)
        {
          uext0[0]=uh.get(0,elemnum-1,RHSvgf.pointsDim());
          uext1[0]=uh.get(1,elemnum-1,RHSvgf.pointsDim());
        }else //inverted reflection
        {
          uext0[0]=0.0;
          uext1[0]=0.0;
        }

      if(elemnum<RHSvgf.pointsDim()-1)
        {
          uext0[1]=uh.get(0,elemnum+1,0);
          uext1[1]=uh.get(0,elemnum+1,0);
        }
      else //inverted reflection
        {
          uext0[1]=0.0;
          uext1[1]=0.0;
        }
  
      Array1D<double> nflux0(2,0.0);
      Array1D<double> nflux1(2,0.0);
      Array1D<double> du0(2,0.0);
      Array1D<double> du1(2,0.0);
      
      for(int i=0; i<2; i++){
        Array2D<double> lambdaminus(2,2,0.0);
        Array2D<double> lambdaplus(2,2,0.0);
        for (int j=0; j<2; j++)
          {
            if (nx[i]*lambda[i][j] <= 0.0) { //might need to re-reverse this
              lambdaminus[j][j]=nx[i]*lambda[i][j];
            }else{
              lambdaplus[j][j]=nx[i]*lambda[i][j];
            }
          }
        nflux0=matmult(lambdaplus,matmult(sinv,uint0));
        nflux1=matmult(lambdaplus,matmult(sinv,uint1));
      }
      for(int i=0; i<2; i++) //A*u at boundaries
        {
          du0[i]=nx[i]*pow(speed,2.0)*uint1[1];
          du1[i]=nx[i]*uint0[1];
        }

    
    
      //rhs including numerical flux
      
      double rx = jacobian[elemnum];
  
      double omega=1.0;
      RHSvgf.set(0,elemnum,pow(speed,2.0)
                 *matmult(rx*refelem.getD(),uh.get(1,elemnum))
                 +matmult(refelem.getLift(),rx*du0));

      RHSvgf.set(1,elemnum,
                 matmult(rx*refelem.getD(),uh.get(0,elemnum))
                 +matmult(refelem.getLift(),rx*du1));
      // 0 and 1 are vectorindices, and uh.get(1,elemnum) is an Array1D
    }
   
}



//tests use old function definition of RK4 routine
// test fourth order polynomial ODE, should be exact in RK4. it was.
/*void RHS(Grid thegrid, VectorGridFunction& uh, 
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
  }*/
 //test harmonic oscillator in time with RK4
 /*void RHS(GridFunction& nodes, VectorGridFunction& uh, 
         VectorGridFunction& RHSvgf, double t)
{
  if(RHSvgf.vectorDim()>2)
    {
      throw invalid_argument("Simple harmonic oscillator ODE should have only two separable components.");
    }

  //du=omega*v
  //dv=-omega*u

  double omega=1.0;
  RHSvgf.set(0,omega*uh.get(1));
  RHSvgf.set(1,-omega*uh.get(0));
  // 0 and 1 are vectorindices, and uh.get(1) is a GridFunction.

   
  }*/

