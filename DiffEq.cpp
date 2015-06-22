#include "DiffEq.h"

/*
Flat spacetime wave equation:
drho/dt = c^2 dpi/dx
dpi/dt=drho/dx
drho/dt=psi


also dpi/dx=psi, but this is not needed for evolution
*/

// A matrix must be formatted such that zero rows are at the top

void setupABmatrices(Grid thegrid&, GridFunction<double>& nodes, 
GridFunction<Array2D<double>>& gfA, VectorGridFunction<Array2D<double>>& vgfB)
{
  double Omega, Omegap, H, Hp, eL, eLp, fT, fTp, fTpp,rm2M;
  GridFunction<Array2D<double>> gfA(nodes.gridDim(), nodes.pointsDim());
  for(int i = 0; i < nodes.gridDim(); i++){
    for(int j = 0; j < nodes.pointsDim(); j++){
      
      //regular wave equation
      //      if(params.metric.flatspacetime){
        Array2D<double> A(3, 3, 0.0);
        A[1][2] = -pow(params.waveeq.speed, 2.0);
        A[2][1] = -1.0;
        gfA.set(i, j, A);
        for(int k = 0; k < thegrid.modesDim(); k++) {
          Array2D<double> B(3, 3, 0.0);
          B[0][1] = -1.0;
          vgfB.set(k, i, j, B);
        }

        /*      } else if (params.metric.schwarzchild) {
        //scri-minus
        if(nodes.get(i, j)==params.grid.Sminus) {
          Omega = 0.0;
          Omegap = 0.0;
          eL = 1.0;
          eLp = 0.0;
          H = -1.0;
          Hp = 0.0;
          thegrid.rstar.set(i,j,empty); //WHAT IS EMPTY?
          thegrid.rschw.set(i,j,2.0*mass); //WHERE IS MASS SET?
          term1 = 0.0;
          term2 = 1.0;
          Array2D<double> A(3, 3, 0.0);
          A[1][2] = -1.0;
          A[2][1] = 0.0;
          A[2][2] = -1.0;
          gfA.set(i, j, A);
          for(int k = 0; k < thegrid.modesDim(); k++) {
            Array2D<double> B(3, 3, 0.0);
            B[0][2] = -1.0;
            vgfB.set(k, i, j, B);
          }
          //inner hyperboloidal layer
        } else if ((nodes.get(i, j) > Sminus) && (nodes.get(i,j) < Rminus)) {
          transition(thegrid.rho.get(i, j), Rminus, Sminus, fT, fTp, fTpp);
          //WRITE TRANSITION
          Omega = 1.0 - thegrid.rho.get(i, j) / Sminus * fT;
          Omegap = -(fT + thegrid.rho.get(i, j) * fTp) / Sminus;
          eL = 1.0 + pow(thegrid.rho.get(i, j), 2.0) * fTp / Sminus;
          eLp = thegrid.rho.get(i, j) 
            * (2.0 * fTp + thegrid.rho.get(i, j) * fTpp) / Sminus;
          H = -1.0 + pow(Omega, 2.0) / eL;
          Hp = (2.0 * Omega * Omegap * eL - pow(Omega, 2.0) * eLp) 
            / pow(eL, 2.0);
          thegrid.rstar.set(i, j, thegrid.rho.get(i,j) / Omega);
          //WRITE INVERT_TORTOISE
          rm2M = invert_tortoise(thegrid.rstar.get(i, j), mass);
          thegrid.rschw.set(i, j, 2.0 * mass + rm2M);
          term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i,j),3.0));
          term2 = 2.0 * mass / thegrid.rschw.get(i,j);
          Array2D<double> A(3, 3, 0.0);
          A[1][2] = -1.0;
          A[2][1] = -(1.0 + H) / (1.0 - H);
          A[2][2] = 2.0 * H / (1.0 - H);
          gfA.set(i,j,A);
          for(int k = 0; k < thegrid.modesDim(); k++) {
             Array2D<double> B(3, 3, 0.0);
             B[0][2] = -1.0;
             B[2][1] = Hp / (1.0 - H);
             B[2][2] = -Hp / (1.0 - H);
             B[2][0] = 1.0 / (1.0 - pow(H,2.0)) * pow(Omega, 2.0) * term1
               *( ll[k] * (ll[k] + 1.0) + term2);
             vgfB.set(k, i, j, B);
          }
          //central tortoise region
        } else if ((thegrid.rho.get(i,j) >= Rminus) 
                   && (thegrid.rho.get(i,j) <= Rplus)){
          Omega = 1.0;
          Omegap = 0.0;
          eL = 1.0;
          eLp = 0.0;
          H = 0.0;
          Hp = 0.0;
          thegrid.rstar.set(i, j, thegrid.rho.get(i, j));
          rm2M = invert_tortoise(thegrid.rstar.get(i, j), mass);
          thegrid.rschw.set(i, j, 2.0 * mass + rm2M);
          term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i, j), 3.0));
          term2 = 2.0 * mass / thegrid.rschw.get(i, j);
          Array2D<double> A(3, 3, 0.0);
          A[1][2] = -1.0;
          gfA.set(i, j, A);
          for(int k = 0; k < thegrid.modesDim(); k++) {
            Array2D<double> B(3, 3, 0.0);
            B[0][2] = -1.0;
            B[2][0] = 1.0 / (1.0 - pow(H, 2.0)) * pow(Omega, 2.0) * term1
              * (ll[k] * (ll[k] + 1.0) + term2);
            vgfB.set(k, i, j, B);
          }
          //outer hyperboloidal region
        } else if ((thegrid.rho.get(i, j) > Rplus) 
                   && (thegrid.rho.get(i, j) < Splus)) { 
          call transition(thegrid.rho.get(i,j), Rplus, Splus, fT, fTp, fTpp);
          Omega = 1.0 - grid.rho.get(i, j) / Splus * fT;
          Omegap = -(fT + grid.rho.get(i, j) * fTp) / Splus;
          eL = 1.0 + pow(grid.rho.get(i, j), 2.0) * fTp / Splus; 
          eLp = grid.rho.get(i, j) 
            * (2.0 * fTp + grid.rho.get(i, j) * fTpp) / Splus;
          H = 1.0 - pow(Omega, 2.0) / eL;
          Hp = -(2.0 * Omega * Omegap * eL - pow(Omega, 2.0) * eLp) 
            / pow(eL, 2.0);
          thegrid.rstar.set(i, j, thegrid.rho.get(i, j) / Omega);
          rm2M = invert_tortoise(thegrid.rstar.get(i, j), mass);
          thegrid.rschw.set(2.0 * mass + rm2M);
          term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i, j)));
          term2 = 2.0 * mass / thegrid.rschw.get(i, j);
          Array2D<double> A(3, 3, 0.0);
          A[1][2] = -1.0;
          A[2][1] = -(1.0 - H) / (1.0 + H);
          A[2][2] = 2.0 * H / (1.0 + H);
          gfA.set(i, j, A);
          for(int k = 0; k < thegrid.modesDim(); k++) {
            Array2D<double> B(3, 3, 0.0);
            B[0][2] = -1.0;
            B[2][2] = Hp / (1.0 + H);
            B[2][1] = Hp / (1.0 + H);
            B[2][0] = 1.0 / (1.0 - pow(H, 2.0)) * pow(Omega, 2.0) * term1
              * (ll[k] * (ll[k] + 1.0) + term2);
            vgfB.set(k, i, j, B);
          }
          //scri-plus
        } else if (grid.rho.get(i,j) == Splus) {
          Omega = 0.0;
          Omegap = 0.0;
          eL = 1.0;
          eLp = 0.0;
          H = 1.0;
          Hp = 0.0;
          thegrid.rstar.set(i, j, -empty); //WHAT"S EMPTY?
          }*/
    }//end inner for            
  }//end outer for
}//end function set AB

/*GridFunction<Array2D<double>> setupBmatrix(GridFunction<double>& nodes)
{
  GridFunction<Array2D<double>> gf(nodes.gridDim(), nodes.pointsDim());
  for(int i = 0; i < nodes.gridDim(); i++){
    for(int j = 0; j < nodes.pointsDim(); j++){
      Array2D<double> B(3, 3, 0.0);
      B[0][1] = -1.0;
      vgfB.set(i, j, B);
    }
  }
  return gfA;
  }*/
