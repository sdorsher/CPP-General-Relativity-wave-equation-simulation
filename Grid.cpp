#include "Grid.h"

//du/dt + A du/dx + Bu = 0

Grid::Grid(int elemorder, int numelements, double lowerlim, double upperlim):
  order{elemorder},
  NumElem{numelements},
  nodeLocs{0,elemorder + 1}, 
  Amatrices{numelements, elemorder + 1},
  Bmatrices{numelements, elemorder + 1},
  refelem{elemorder},
  trimmedAmatrices(numelements, elemorder + 1)

{

  //assign evenly spaced element boundaries
  for(int i = 0; i <= numelements; i++) {
    elementBoundaries.push_back(lowerlim + i * (upperlim - lowerlim) 
                                / float(numelements));
  }
  
  //Get physical positions of nodes from the reference element
  Array1D<double> physicalPosition(elemorder + 1);
  for(int elem = 0; elem < numelements; elem++){
    physicalPosition = ((elementBoundaries[elem + 1] 
                         - elementBoundaries[elem]) / 2.0)
      *refelem.getr()
      +((elementBoundaries[elem + 1] + elementBoundaries[elem]) / 2.0);
    nodeLocs.append(physicalPosition);
  }
  
  //Calculate the jacobian associated with the transformation each element
  //from the reference element to physical space
  calcjacobian();
  
  //Setup the A and B matrices in DiffEq.cpp
  Amatrices = setupAmatrix(nodeLocs);
  Bmatrices = setupBmatrix(nodeLocs);

  //Get the A matrix with its zero dimensions removed for each node
  for(int i = 0; i < nodeLocs.gridDim(); i++){
    for(int j = 0; j < nodeLocs.pointsDim(); j++)
      {
        CharacteristicFlux nodechar(Amatrices.get(i,j));
        trimmedAmatrices.set(i, j, (nodechar.getAtrimmed()));
      }
  }

  //Get all characteristic equation information for each boundary node
  for (int i = 0; i < nodeLocs.gridDim(); i++){
    CharacteristicFlux left(Amatrices.get(i, 0));
    CharacteristicFlux right(Amatrices.get(i, nodeLocs.pointsDim() - 1));
    AleftBoundaries.push_back(left);
    ArightBoundaries.push_back(right);
  }
  
  //Allocate memory for the left and right boundary du that contributes
  //to the characteristic flux when multiplied by the lift matrix. 
  //There will be one Array1D of length 2 for each boundary of each element.
  //this is not a GridFunction because although it shares the same 
  //implementation format, it doesn't share the same conceptual format. 
  //It does not embody a function that has values over all nodes of the grid.
  duL.resize(NumElem);
  duR.resize(NumElem);

}

vector<TNT::Array2D<double>> 
Grid::characteristicflux(VectorGridFunction<double>& uh)
{
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
    
    uintL = uh.getVectorAsArray1D(elemnum, indL, vminL, vmaxL); 
    uintR = uh.getVectorAsArray1D(elemnum, indR, vminR, vmaxR);

    
    if(elemnum > 0) {
      uextL = uh.getVectorAsArray1D(elemnum - 1, indR, vminL, vmaxL); 
      //external u, left boundary
    }else{
      uextL = uh.getVectorAsArray1D(NumElem - 1, indR, vminL, vmaxL); 
      //periodic boundary conditions
    }

    if(elemnum < NumElem - 1) {
      uextR = uh.getVectorAsArray1D(elemnum + 1, indL, vminR, vmaxR); 
      //external u, right boundary
    }else{
      uextR = uh.getVectorAsArray1D(0, indL, vminR, vmaxR); 
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
    Array1D<double> nfluxL = matmult(lambdaplusL, matmult(sinvL, uintL));
    nfluxL += matmult(lambdaminusL, matmult(sinvL, uextL));
    nfluxL = matmult(SL, nfluxL);

    Array1D<double> nfluxR = matmult(lambdaplusR, matmult(sinvR, uintR));
    nfluxR += matmult(lambdaminusR, matmult(sinvR, uextR));
    nfluxR = matmult(SR, nfluxR);

    Array2D<double> AtrimmedL= AleftBoundaries[elemnum].getAtrimmed();
    Array2D<double> AtrimmedR= AleftBoundaries[elemnum].getAtrimmed();
    Array2D<double> duelem(AtrimmedR.dim1(), 2, 0.0);


    //This gets multiplied by lift matrix to calculate flux
    Array1D<double> duL = nL * matmult(AtrimmedL, uintL) - nfluxL; 
    Array1D<double> duR = nR * matmult(AtrimmedR, uintR) - nfluxR; 

    insert_1D_into_2D(duelem, duL, 0, false);
    insert_1D_into_2D(duelem, duR, 1, false);

    du[elemnum] = duelem;
  }
  return du;
}

void Grid::RHS(VectorGridFunction<double>& uh, 
               VectorGridFunction<double>& RHSvgf, double t, 
               vector<Array2D<double>>& du )
{

  for(int elemnum = 0; elemnum < NumElem; elemnum++){
    //Maximum index for both A and B matrix
    int vmaxAB = ArightBoundaries[elemnum].getAdim() - 1;
    //Minimum index for use with trimmed A matrix. 
    //Minimum index for B matrix is zero
    int vminA = vmaxAB - ArightBoundaries[elemnum].getDdim() + 1;

    //The B matrix component of the RHS. 
    Array2D<double> RHSB(uh.pointsDim(), ArightBoundaries[elemnum].getAdim());
                         
    for(int nodenum = 0; nodenum < uh.pointsDim(); nodenum++){
      Array1D<double> RHSBpernode;
      //Multiply the B matrix times a "vector" of u for each node
      RHSBpernode = matmult(Bmatrices.get(elemnum, nodenum),
                          uh.getVectorAsArray1D(elemnum, nodenum, 0, vmaxAB));
      //Insert that result into the rows of a larger matrix
      insert_1D_into_2D(RHSB, RHSBpernode, nodenum, false);

    }//This can be sped up by skipping the insert step and reading directly 
    //from per node


    //A contribution:
    Array2D<double> RHSA1(uh.pointsDim(), ArightBoundaries[elemnum].getDdim());
    
    //The A contribution needs to be multiplied one node at a time by the
    //trimmed A matrix in a similar manner to the B contribution. But first,
    //we take the spatial derivative across all nodes. 
    Array2D<double> RHSA1preA = jacobian(elemnum) * (matmult(refelem.getD(),
                            uh.getVectorNodeArray2D(elemnum, vminA, vmaxAB)));
    
    
    //Multiply each row of RHSA1preA by a different a Atrimmed matrix
    for(int nodenum=0; nodenum < uh.pointsDim(); nodenum++)
      {
        int M = trimmedAmatrices.get(elemnum,nodenum).dim1();
        int N = trimmedAmatrices.get(elemnum,nodenum).dim2();
        int K = RHSA1preA.dim1();
        int L = RHSA1preA.dim2();

        Array2D<double> tA = trimmedAmatrices.get(elemnum, nodenum);
            
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
    Array2D<double> RHSA2 = jacobian(elemnum) 
      * matmult(refelem.getLift(), du[elemnum]);

    //RHSA and RHSB will have different sizes due to the different 
    //number of diffeq variables stored in each. sum them using a 
    //for loop while assigning values to the RHSvgf vector grid function

    Array2D<double> RHSA = RHSA1 + RHSA2;

    //Sum the contributions from B, derivative, and flux, 
    //accounting for different matrix dimensions
    for(int vecnum = 0; vecnum < RHSvgf.vectorDim(); vecnum++){
        for(int nodenum = 0; nodenum < RHSvgf.pointsDim(); nodenum++){
          if(vecnum<vminA){
            RHSvgf.set(vecnum,elemnum,nodenum,-RHSB[nodenum][vecnum]);
          } else {
            RHSvgf.set(vecnum, elemnum, nodenum, -RHSB[nodenum][vecnum]
                         + RHSA[nodenum][vecnum - vminA]);
          }//-sign in B because it is on the left hand side of the 
          //equation in the definition supplied in DiffEq.cpp
        }
    }
  }

}

GridFunction<double> Grid::gridNodeLocations()
{
  return nodeLocs;
}

vector<double> Grid::gridBoundaries()
{
  return elementBoundaries;
}

void Grid::calcjacobian()
{
  for(int elem = 0; elem < NumElem; elem++){
    double rx = 2.0 / (elementBoundaries[elem + 1] 
                       - elementBoundaries[elem]);
    drdx.push_back(rx);
  }
}

double Grid::jacobian(int elemnum)
{
  return drdx[elemnum];
}
