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
  RHS(thegrid, uh, RHSvgf, t,true);
  k=deltat*RHSvgf;
  uh=uh+rk4b[0]*k;



  //step2
  for(int i=2; i<=5; i++){
    RHS(thegrid, uh, RHSvgf, t+rk4c[i-1]*deltat,false);
    k=rk4a[i-1]*k+deltat*RHSvgf;
    uh=uh+rk4b[i-1]*k;
  }

}

//---------------------

void RHS(Grid thegrid, VectorGridFunction& uh, 
         VectorGridFunction& RHSvgf, double t,bool output)
{
  if(RHSvgf.vectorDim()>3)
    {
      throw invalid_argument("Wave equation PDE should only have three separable components.");
    }

  GridFunction nodes(thegrid.gridNodeLocations().gridDim(),thegrid.gridNodeLocations().pointsDim(),false);
  nodes=thegrid.gridNodeLocations();

  // uh=(psi,rho,pi)

  //wave equation specific values

  //drho/dt - -c^2dPi/dx =0
  //dpi/dt - drho/dx =0

  // A=[[0 -c^2],[-1 0]]
  //LAMBDA=S^-1 PI S=[[-c 0],[0 c]]
  // S=[v1 v2] so that eigenvalues are diagonal
  // eigenvalues are +-c
  // S=[[c -c],[1 1]]
  // 


  

  Array2D<double> s(2,2);
  Array2D<double> sinv(2,2);
  Array2D<double> lambd(2,2);
   Array2D<double> A(2,2);

  
   A[0][0]=0.0;
   A[0][1]=-params.speed*params.speed;
   A[1][0]=-1.0;
   A[1][1]=0.0;

  
  s[0][0]=params.speed;
  s[0][1]=-params.speed;
  s[1][0]=1.0;
  s[1][1]=1.0;
  
  sinv[0][0]=0.5/params.speed;
  sinv[0][1]=0.5;
  sinv[1][0]=-0.5/params.speed;
  sinv[1][1]=0.5;

  lambd[0][0]=-params.speed;
  lambd[1][0]=-params.speed;
  lambd[0][1]=params.speed;
  lambd[1][1]=params.speed;
  //lambda is eigenvalues of characteristic matrix at element boundaries

  //   output2D(lambd);

  //  output2D(matmult(sinv,matmult(A,s)));

  vector<double> jacobian=thegrid.jacobian();

  //numerical flux
  
  //index
  vector<int> ind{0,RHSvgf.pointsDim()-1};

  for(int elemnum=0; elemnum<uh.gridDim(); elemnum++)
    {
  
      //interior,exterior for three variables of wave equation
      Array1D<double> uint1(2,0.0);
      Array1D<double> uint0(2,0.0);
      Array1D<double> uext1(2,0.0);
      Array1D<double> uext0(2,0.0);
      for(int i=0; i<2; i++)
        {
          uint0[i] = uh.get(1,elemnum,ind[i]);
          uint1[i]= uh.get(2,elemnum,ind[i]);
        }
      Array1D<double> nx(2,1.0);
      nx[0]=-1.0;
      
      if(elemnum>0)
        {
          uext0[0]=uh.get(1,elemnum-1,RHSvgf.pointsDim()-1);
          uext1[0]=uh.get(2,elemnum-1,RHSvgf.pointsDim()-1);
        }        
      else 
        {
          uext0[0]=uh.get(1,RHSvgf.gridDim()-1,RHSvgf.pointsDim()-1);//0.0;
          uext1[0]=uh.get(2,RHSvgf.gridDim()-1,RHSvgf.pointsDim()-1);//0.0;
        }

      if(elemnum<RHSvgf.gridDim()-1)
        {
          uext0[1]=uh.get(1,elemnum+1,0);
          uext1[1]=uh.get(2,elemnum+1,0);
        }
      else 
        {
          uext0[1]=uh.get(1,0,0);//0.0;
          uext1[1]=uh.get(2,0,0);//0.0;
        }
      if(elemnum==0){
        //  cout << "LHS0 uext=" << uext0[0] << " uint=" << uint0[0] <<endl;
        // cout << "RHS0 uext=" <<uext0[1] << " uint=" << uint0[1] << endl;
        // cout << "LHS1 uext=" << uext1[0] << " uint=" << uint1[0] <<endl;
        // cout << "RHS1 uext=" <<uext1[1] << " uint=" << uint1[1] << endl;
      }
      Array1D<double> nfluxL(2,0.0);
      Array1D<double> nfluxR(2,0.0);
      Array1D<double> du0(2,0.0);
      Array1D<double> du1(2,0.0);
      Array1D<double> uintL(2);
      Array1D<double> uintR(2);
      Array1D<double> uextL(2);
      Array1D<double> uextR(2);
      uintL[0]=uint0[0];
      uintL[1]=uint1[0];
      uintR[0]=uint0[1];
      uintR[1]=uint1[1];
      uextL[0]=uext0[0];
      uextL[1]=uext1[0];
      uextR[0]=uext0[1];
      uextR[1]=uext1[1];

      for(int i=0; i<2; i++){
        Array2D<double> lambdminus(2,2,0.0);
        Array2D<double> lambdplus(2,2,0.0);

        for (int j=0; j<2; j++)
          {
            if (nx[i]*lambd[i][j] <= 0.0) { //might need to re-reverse this
              lambdminus[j][j]=nx[i]*lambd[i][j];
            }else{
              lambdplus[j][j]=nx[i]*lambd[i][j];
            }
          }
        if((elemnum==5)&&(i==0))
          {
            //output2D(lambdplus);
            //output2D(lambdminus);
            //            output1D(uintL);//okay
            //output1D(uintR);//okay

          }

        if(i==0)
          {
            nfluxL=matmult(lambdplus,matmult(sinv,uintL));
            nfluxL+=matmult(lambdminus,matmult(sinv,uextL));
            nfluxL=matmult(s,nfluxL);
             if(elemnum==5)
              {
                //         output1D(nfluxL);
              }

          }else{
          nfluxR=matmult(lambdplus,matmult(sinv,uintR));
          nfluxR+=matmult(lambdminus,matmult(sinv,uextR));
          nfluxR=matmult(s,nfluxR);
        }
      }
      Array1D<double> nflux1(2);
      Array1D<double> nflux0(2);
      nflux0[0]=nfluxL[0];
      nflux0[1]=nfluxR[0];
      nflux1[0]=nfluxL[1];
      nflux1[1]=nfluxR[1];
      for(int i=0; i<2; i++) //A*u at boundaries
        {//correct. don't arbitratrily change the sign of this
          du0[i]=-nx[i]*pow(params.speed,2.0)*uint1[i]-nflux0[i];
          du1[i]=-nx[i]*uint0[i]-nflux1[i];
        }
      if(elemnum==5)
        {
          // cout << "du " << du0[0] << " " << du0[1] << " " << du1[0] << " " << du1[1] <<endl;
        }
      //rhs including numerical flux
      
      double rx = jacobian[elemnum];
  
      if(elemnum==0)
        {
          //          cout << uh.get(2,elemnum,0) <<" "<<uh.get(2,elemnum,uh.pointsDim()-1)<< endl;
          //        cout << rx << endl;
        }

      double omega=1.0;
      RHSvgf.set(0,elemnum,uh.get(1,elemnum));
      if(output)
        {
          //          std::ofstream fs;
          // fs.open("rhsoutput.txt", ios::app);
  
          Array1D<double> RHS1(RHSvgf.pointsDim());
          RHS1=pow(params.speed,2.0)
            *matmult(rx*thegrid.refelem.getD(),uh.get(2,elemnum));
          Array1D<double> RHS2(RHSvgf.pointsDim());
          RHS2=matmult(rx*thegrid.refelem.getD(),uh.get(1,elemnum));
          for(int k=0; k<RHS1.dim(); k++)
            {
              //     cout << nodes.get(elemnum,k) << " "<< RHS1[k] << " " << RHS2[k] << endl;
              // fs << nodes.get(elemnum,k) << " " << RHS1[k] << " " << RHS2[k] << endl;
            }
          //fs.close();

        }
      RHSvgf.set(1,elemnum,pow(params.speed,2.0)
                 *matmult(rx*thegrid.refelem.getD(),uh.get(2,elemnum))
                 +matmult(thegrid.refelem.getLift(),rx*du0));

      RHSvgf.set(2,elemnum,
                 matmult(rx*thegrid.refelem.getD(),uh.get(1,elemnum))
                 +matmult(thegrid.refelem.getLift(),rx*du1));

      

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

