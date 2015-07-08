#include "DiffEq.h"

/*
Flat spacetime wave equation:
drho/dt = c^2 dpi/dx
dpi/dt=drho/dx
drho/dt=psi


also dpi/dx=psi, but this is not needed for evolution
*/


Array2D<double> DiffEq::getA(int gridindex, int pointsindex)
{
  return Amatrices.get(gridindex, pointsindex);
}

Array2D<double> DiffEq::getB(int modesindex, int gridindex, int pointsindex)
{
  return Bmatrices.get(modesindex, gridindex, pointsindex);
}

CharacteristicFlux DiffEq::getAleft(int elemnum)
{
  return AleftBoundaries.at(elemnum);
}

CharacteristicFlux DiffEq::getAright(int elemnum)
{
  return ArightBoundaries.at(elemnum);
}

Array2D<double> DiffEq::getAtrimmed(int gridindex, int pointsindex)
{
  return trimmedAmatrices.get(gridindex, pointsindex);
}

// The A matrix must be formatted such that zero rows are at the top

DiffEq::DiffEq(Grid& thegrid, Modes& lmmodes):
  Amatrices{params.grid.numelems, params.grid.elemorder + 1},
  Bmatrices{params.waveeq.modenum, params.grid.numelems, 
      params.grid.elemorder + 1},
  trimmedAmatrices{params.grid.numelems, params.grid.elemorder + 1}
{
  //set up the A and B matrices
  setupABmatrices(thegrid, lmmodes);

  //  GridFunction<double> nodeLocs = thegrid.gridNodeLocations();

  //Get the A matrix with its zero dimensions removed for each node
  for(int i = 0; i < thegrid.numberElements(); i++){
    for(int j = 0; j <= thegrid.nodeOrder(); j++) {
      //get the trimmed A matrices at each node coordinate
      CharacteristicFlux nodechar(Amatrices.get(i,j));
      trimmedAmatrices.set(i, j, (nodechar.getAtrimmed()));
    }
  }

  //Get all characteristic equation information for each boundary node
  for (int i = 0; i < thegrid.numberElements(); i++){
    CharacteristicFlux left(Amatrices.get(i, 0));
    CharacteristicFlux right(Amatrices.get(i, thegrid.nodeOrder()));
    AleftBoundaries.push_back(left);
    ArightBoundaries.push_back(right);
  }
}

void DiffEq::setupABmatrices(Grid& thegrid, Modes& lmmodes)
{
  double Omega, Omegap, H, Hp, eL, eLp, fT, fTp, fTpp,rm2M;
  GridFunction<double> nodes = thegrid.gridNodeLocations();

  for(int i = 0; i < nodes.gridDim(); i++){
    for(int j = 0; j < nodes.pointsDim(); j++){
      //regular wave equation
      if(params.metric.flatspacetime){
        Array2D<double> A(3, 3, 0.0);
        A[1][2] = -pow(params.waveeq.speed, 2.0);
        A[2][1] = -1.0;
        Amatrices.set(i, j, A);
        //vector dimension of B is actually number of modes
        for(int k = 0; k < Bmatrices.vectorDim(); k++) {
          Array2D<double> B(3, 3, 0.0);
          B[0][1] = -1.0;
          Bmatrices.set(k, i, j, B);
        }

       } else if (params.metric.schwarzchild) {
        double Omega, Omegap, eL, eLp, H, Hp, term1, term2;
        //nodes = rho
        
        //scri-minus
        if(nodes.get(i, j)==params.hyperb.Sminus) {
          Omega = 0.0;
          Omegap = 0.0;
          eL = 1.0;
          eLp = 0.0;
          H = -1.0;
          Hp = 0.0;
          thegrid.rstar.set(i,j,DBL_MAX); 
          thegrid.rschw.set(i,j,2.0*params.schw.mass); 
          term1 = 0.0;
          term2 = 1.0;
          Array2D<double> A(3, 3, 0.0);
          A[1][2] = -1.0;
          A[2][1] = 0.0;
          A[2][2] = -1.0;
          Amatrices.set(i, j, A);
          for(int k = 0; k < lmmodes.ntotal; k++) {
            Array2D<double> B(3, 3, 0.0);
            B[0][2] = -1.0;
            Bmatrices.set(k, i, j, B);
          }
          //inner hyperboloidal layer
        } else if ((nodes.get(i, j) > params.hyperb.Sminus) 
                   && (nodes.get(i,j) < params.hyperb.Rminus)) {
          transition(nodes.get(i, j), params.hyperb.Rminus, 
                     params.hyperb.Sminus, fT, fTp, fTpp);
          Omega = 1.0 - nodes.get(i, j) / params.hyperb.Sminus * fT;
          Omegap = -(fT + nodes.get(i, j) * fTp) / params.hyperb.Sminus;
          eL = 1.0 + pow(nodes.get(i, j), 2.0) * fTp 
            / params.hyperb.Sminus;
          eLp = nodes.get(i, j) 
            * (2.0 * fTp + nodes.get(i, j) * fTpp)
            / params.hyperb.Sminus;
          H = -1.0 + pow(Omega, 2.0) / eL;
          Hp = (2.0 * Omega * Omegap * eL - pow(Omega, 2.0) * eLp) 
            / pow(eL, 2.0);
          thegrid.rstar.set(i, j, nodes.get(i,j) / Omega);
          rm2M = invert_tortoise(thegrid.rstar.get(i, j), params.schw.mass);
          thegrid.rschw.set(i, j, 2.0 * params.schw.mass + rm2M);
          term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i,j),3.0));
          term2 = 2.0 * params.schw.mass / thegrid.rschw.get(i,j);
          Array2D<double> A(3, 3, 0.0);
          A[1][2] = -1.0;
          A[2][1] = -(1.0 + H) / (1.0 - H);
          A[2][2] = 2.0 * H / (1.0 - H);
          Amatrices.set(i,j,A);
          for(int k = 0; k < lmmodes.ntotal; k++) {
             Array2D<double> B(3, 3, 0.0);
             B[0][2] = -1.0;
             B[2][1] = Hp / (1.0 - H);
             B[2][2] = -Hp / (1.0 - H);
             B[2][0] = 1.0 / (1.0 - pow(H,2.0)) * pow(Omega, 2.0) * term1
               *( lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
             Bmatrices.set(k, i, j, B);
          }
          //central tortoise region
        } else if ((nodes.get(i,j) >= params.hyperb.Rminus) 
                   && (nodes.get(i,j) <= params.hyperb.Rplus)){
          Omega = 1.0;
          Omegap = 0.0;
          eL = 1.0;
          eLp = 0.0;
          H = 0.0;
          Hp = 0.0;
          thegrid.rstar.set(i, j, nodes.get(i, j));
          rm2M = invert_tortoise(thegrid.rstar.get(i, j), params.schw.mass);
          thegrid.rschw.set(i, j, 2.0 * params.schw.mass + rm2M);
          term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i, j), 3.0));
          term2 = 2.0 * params.schw.mass / thegrid.rschw.get(i, j);
          Array2D<double> A(3, 3, 0.0);
          A[1][2] = -1.0;
          Amatrices.set(i, j, A);
          for(int k = 0; k < lmmodes.ntotal; k++) {
            Array2D<double> B(3, 3, 0.0);
            B[0][2] = -1.0;
            B[2][0] = 1.0 / (1.0 - pow(H, 2.0)) * pow(Omega, 2.0) * term1
              * (lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
            Bmatrices.set(k, i, j, B);
          }
          //outer hyperboloidal region
        } else if ((nodes.get(i, j) > params.hyperb.Rplus) 
                   && (nodes.get(i, j) < params.hyperb.Splus)) { 
          call transition(nodes.get(i,j), params.hyperb.Rplus, params.hyperb.Splus, fT, fTp, fTpp);
          Omega = 1.0 - nodes.get(i, j) / params.hyperb.Splus * fT;
          Omegap = -(fT + nodes.get(i, j) * fTp) / params.hyperb.Splus;
          eL = 1.0 + pow(nodes.get(i, j), 2.0) * fTp / params.hyperb.Splus; 
          eLp = nodes.get(i, j) 
            * (2.0 * fTp + nodes.get(i, j) * fTpp) / params.hyperb.Splus;
          H = 1.0 - pow(Omega, 2.0) / eL;
          Hp = -(2.0 * Omega * Omegap * eL - pow(Omega, 2.0) * eLp) 
            / pow(eL, 2.0);
          thegrid.rstar.set(i, j, nodes.get(i, j) / Omega);
          rm2M = invert_tortoise(thegrid.rstar.get(i, j), params.schw.mass);
          thegrid.rschw.set(2.0 * params.schw.mass + rm2M);
          term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i, j)));
          term2 = 2.0 * params.schw.mass / thegrid.rschw.get(i, j);
          Array2D<double> A(3, 3, 0.0);
          A[1][2] = -1.0;
          A[2][1] = -(1.0 - H) / (1.0 + H);
          A[2][2] = 2.0 * H / (1.0 + H);
          Amatrices.set(i, j, A);
          for(int k = 0; k < lmmodes.ntotal; k++) {
            Array2D<double> B(3, 3, 0.0);
            B[0][2] = -1.0;
            B[2][2] = Hp / (1.0 + H);
            B[2][1] = Hp / (1.0 + H);
            B[2][0] = 1.0 / (1.0 - pow(H, 2.0)) * pow(Omega, 2.0) * term1
              * (lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
            Bmatrices.set(k, i, j, B);
          }
          //scri-plus
        } else if (nodes.get(i,j) == params.hyperb.Splus) {
          Omega = 0.0;
          Omegap = 0.0;
          eL = 1.0;
          eLp = 0.0;
          H = 1.0;
          Hp = 0.0;
          thegrid.rstar.set(i, j, -DBL_MAX); 
          }*/
    }//end inner for            
  }//end outer for
}//end function set AB


vector<TNT::Array2D<double>> 
DiffEq::characteristicflux(int modenum, Grid& thegrid,
                           TwoDVectorGridFunction<double>& uh)
{
  //We can loop over characteristicFlux in an external function because
  //in general, neither the RHS of the differential equation nor du mixes
  //spherical harmonic modes. 

  int NumElem = thegrid.numberElements();

  vector<Array2D<double>> du;
  du.resize(NumElem);
  for(int elemnum=0; elemnum<NumElem; elemnum++){
    int indL = 0; //index of leftmost node of that element
    int indR = uh.pointsDim()-1; //index of rightmost node of that element
    double nL = -1.0; //normal to the leftmost node
    double nR = 1.0; //normal to the rightmost node

    //Dimension of the components of the differential equation with 
    //spatial derivatives (dimension of the trimmed A matrices)
    int DdimL = AleftBoundaries[elemnum].getDdim();
    int DdimR = ArightBoundaries[elemnum].getDdim();
    
    //vmin and vmax are min and max indices in vector dimension (psi, rho, pi)
    int vmaxL = AleftBoundaries[elemnum].getAdim() - 1;
    int vminL = vmaxL - DdimL + 1; //neglect zero rows at top of A matrix
    int vmaxR = ArightBoundaries[elemnum].getAdim() - 1;
    int vminR = vmaxR - DdimR + 1; //neglect zero rows at top of A matrix

    Array1D<double> uintL(DdimL); //internal u at left boundary
    Array1D<double> uintR(DdimR); //internal u at right boundary
    Array1D<double> uextL(DdimL); //external u at left boundary
    Array1D<double> uextR(DdimR); //external u at right boundary
    
    uintL = uh.getVectorAsArray1D(modenum,elemnum, indL, vminL, vmaxL, 0); 
    uintR = uh.getVectorAsArray1D(modenum,elemnum, indR, vminR, vmaxR, 0);

    
    if(elemnum > 0) {
      uextL = uh.getVectorAsArray1D(modenum, elemnum - 1, indR, vminL, 
                                    vmaxL, 0); 
      //external u, left boundary
    }else{
      uextL = uh.getVectorAsArray1D(modenum, NumElem - 1, indR, vminL, 
                                    vmaxL, 0); 
      //periodic boundary conditions
    }

    if(elemnum < NumElem - 1) {
      uextR = uh.getVectorAsArray1D(modenum, elemnum + 1, indL, vminR, 
                                    vmaxR, 0); 
      //external u, right boundary
    }else{
      uextR = uh.getVectorAsArray1D(modenum, 0, indL, vminR, vmaxR, 0); 
      //periodic boundary conditions
    }
    
    //Initialize plus and minus components of lambda matrix to zero at both
    //boundaries
    Array2D<double> lambdaminusL(DdimL, DdimL, 0.0);
    Array2D<double> lambdaminusR(DdimR, DdimR, 0.0);
    Array2D<double> lambdaplusL(DdimL, DdimL, 0.0);
    Array2D<double> lambdaplusR(DdimR, DdimR, 0.0);

    Array2D<double> lambdaL= AleftBoundaries[elemnum].getLambda();
    Array2D<double> lambdaR= ArightBoundaries[elemnum].getLambda();

    //lambda minus contains outward moving wave components
    //lambda plus contains inward moving wave components
    //Might be an incorrect summary. Trust the math, not the words
    //See pg 35 of Hesthaven and Warburten
    for(int j = 0; j < DdimL; j++) {
      if(nL * lambdaL[j][j] <= 0) {
        lambdaminusL[j][j] = nL * lambdaL[j][j];
      } else {
        lambdaplusL[j][j] = nL * lambdaL[j][j];
      }
      
      if(nR * lambdaR[j][j] <= 0) {
        lambdaminusR[j][j] = nR * lambdaR[j][j];
      } else {
        lambdaplusR[j][j] = nR * lambdaR[j][j];
      }
    }

    //S and S inverse matrices at both boundaries
    Array2D<double> sinvL = AleftBoundaries[elemnum].getSinv();
    Array2D<double> sinvR = ArightBoundaries[elemnum].getSinv();
    Array2D<double> SL = AleftBoundaries[elemnum].getS();
    Array2D<double> SR = ArightBoundaries[elemnum].getS();
    
    //Numerical fluxes at both boundaries 
    //See Hesthaven and Warburten pg 35 (n*F)
    Array1D<double> nfluxL1 = matmult(lambdaplusL, matmult(sinvL, uintL));
    Array1D<double> nfluxL2 =  matmult(lambdaminusL, matmult(sinvL, uextL));
    Array1D<double> nfluxL = matmult(SL, nfluxL1+nfluxL2);

    Array1D<double> nfluxR1 = matmult(lambdaplusR, matmult(sinvR, uintR));
    Array1D<double> nfluxR2 = matmult(lambdaminusR, matmult(sinvR, uextR));
    Array1D<double> nfluxR = matmult(SR, nfluxR1 + nfluxR2);

    Array2D<double> AtrimmedL= AleftBoundaries[elemnum].getAtrimmed();
    Array2D<double> AtrimmedR= AleftBoundaries[elemnum].getAtrimmed();
    Array2D<double> duelem(AtrimmedR.dim1(), 2, 0.0);
    
    //This gets multiplied by lift matrix to calculate flux
    Array1D<double> duL = nL * matmult(AtrimmedL, uintL) - nfluxL; 
    Array1D<double> duR = nR * matmult(AtrimmedR, uintR) - nfluxR; 

    //create a 2D array with one column corresponding to the left
    //boundary and one column corresponding to the right boundary
    //of du for the the various components of u. 
    insert_1D_into_2D(duelem, duL, 0, false);
    insert_1D_into_2D(duelem, duR, 1, false);

    du[elemnum] = duelem;
  }
  return du;
}

void DiffEq::RHS(int modenum, Grid& thegrid,
                 TwoDVectorGridFunction<double>& uh, 
                 TwoDVectorGridFunction<double>& RHStdvgf, double t, 
                 vector<Array2D<double>>& du )
{

  //We can loop over RHS in an external function independent of mode 
  //because in general the right hand side of the differential equation
  //does not mix spherical harmonic modes. Neither does the numerical flux. 
   

  int NumElem = thegrid.numberElements();

  for(int elemnum = 0; elemnum < NumElem; elemnum++){
    //Maximum index for both A and B matrix
    int vmaxAB = ArightBoundaries[elemnum].getAdim() - 1;
    //Minimum index for use with trimmed A matrix. 
    //Minimum index for B matrix is zero
    int vminA = vmaxAB - ArightBoundaries[elemnum].getDdim() + 1;
    
    //The B matrix component of the RHS. 
    Array2D<double> RHSB(uh.pointsDim(), 
                           ArightBoundaries[elemnum].getAdim());
    for(int nodenum = 0; nodenum < uh.pointsDim(); nodenum++){
      Array1D<double> RHSBpernode;
      //Multiply the B matrix times a "vector" of u for each node
      RHSBpernode = matmult(Bmatrices.get(modenum, elemnum, nodenum),
                            uh.getVectorAsArray1D(modenum,elemnum, nodenum, 
                                                  0, vmaxAB, 0));
      //Insert that result into the rows of a larger matrix
      insert_1D_into_2D(RHSB, RHSBpernode, nodenum, false);
          
    }//This can be sped up by skipping the insert step and reading directly 
    //from per node

    //A contribution:
    Array2D<double> RHSA1(uh.pointsDim(), ArightBoundaries[elemnum].getDdim());
    
    //The A contribution needs to be multiplied one node at a time by the
    //trimmed A matrix in a similar manner to the B contribution. But first,
    //we take the spatial derivative across all nodes. 

    Array2D<double> RHSA1preA = thegrid.jacobian(elemnum) 
      * (matmult(thegrid.refelem.getD(),
                 uh.getVectorNodeArray2D(modenum, elemnum, vminA, vmaxAB, 0)));

    
    //Multiply each row of RHSA1preA by a different a Atrimmed matrix
    for(int nodenum=0; nodenum < uh.pointsDim(); nodenum++)
      {
        int M = trimmedAmatrices.get(elemnum,nodenum).dim1();
        int N = trimmedAmatrices.get(elemnum,nodenum).dim2();
        int K = RHSA1preA.dim1();
        int L = RHSA1preA.dim2();

        Array2D<double> tA(M,N); 
        tA= trimmedAmatrices.get(elemnum, nodenum);

        for (int i=0; i<M; i++){
            double sum = 0;
            for (int k=0; k<N; k++)
              sum += -tA[i][k] * RHSA1preA[nodenum][k];
            RHSA1[nodenum][i] = sum;
          }
        //Copied and pasted from TmatmultT in TNT2 a
        //with modification of variable first matrix
        //TmatmultT was copied and pasted from matmult in tnt itself
      }
    //Negative sign is because of definition of tA
    //A is definied to appear on the left hand side of the differential
    //equation, but this routine calculates the right hand side

    //Needs a multiplication by an A matrix before D 
    //but A is position dependent. 
    
    //This is the contribution due to du, or the numerical flux
    Array2D<double> RHSA2 = thegrid.jacobian(elemnum) 
      * matmult(thegrid.refelem.getLift(), du[elemnum]);

    

    //RHSA and RHSB will have different sizes due to the different 
    //number of diffeq variables stored in each. sum them using a 
    //for loop while assigning values to the RHStdvgf vector grid function

    Array2D<double> RHSA = RHSA1 + RHSA2;

    //Sum the contributions from B, derivative, and flux, 
    //accounting for different matrix dimensions
    for(int vecnum = 0; vecnum < RHStdvgf.vectorDim(); vecnum++){
        for(int nodenum = 0; nodenum < RHStdvgf.pointsDim(); nodenum++){
          if(vecnum<vminA){
            RHStdvgf.set(modenum, vecnum,elemnum,nodenum,-RHSB[nodenum][vecnum]);
          } else {
            RHStdvgf.set(modenum, vecnum, elemnum, nodenum, -RHSB[nodenum][vecnum]
                         + RHSA[nodenum][vecnum - vminA]);
          }//-sign in B because it is on the left hand side of the 
          //equation in the definition supplied in DiffEq.cpp
        }
    }
  }

}

void DiffEq::modeRHS(Grid& thegrid,
                     TwoDVectorGridFunction<double>& uh,
                     TwoDVectorGridFunction<double>& RHStdvgf, 
                     double t)
{

  for(int modenum = 0; modenum < uh.modesDim(); modenum++) {
    vector<Array2D<double>> du;
    du = move(characteristicflux(modenum, thegrid, uh));
    RHS(modenum, thegrid, uh, RHStdvgf, t, du);
  }
}

