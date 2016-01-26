#include "DiffEq.h"

/*
Fla spacetime wave equation:
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

DiffEq::DiffEq(Grid& thegrid, Modes& lmmodes, int nmodetotal):
  Amatrices{params.grid.numelems, params.grid.elemorder + 1},
  Bmatrices{nmodetotal, params.grid.numelems, 
      params.grid.elemorder + 1},
  trimmedAmatrices{params.grid.numelems, params.grid.elemorder + 1},
  source{nmodetotal, params.grid.numelems, params.grid.elemorder+1,{0.0,0.0}},
  window{params.grid.numelems, params.grid.elemorder+1},
  dwindow{params.grid.numelems, params.grid.elemorder+1},
  d2window{params.grid.numelems, params.grid.elemorder+1}

  {
    
  //set up the A and B matrices


  setupABmatrices(thegrid, lmmodes);
  

  if(params.metric.schwarschild){
  //Print out the A and B coefficients in the same format as the Fortran file
   ofstream fs;
  fs.open("ABcoeffs.txt");
  for(int i= 0; i<params.grid.numelems; i++){
    for(int j =0; j <= params.grid.elemorder; j++){
      Array2D<double> A = Amatrices.get(i,j);
      Array2D<double> tA = trimmedAmatrices.get(i,j);
      Array2D<double> B0 = Bmatrices.get(0,i,j);
      Array2D<double> B1 = Bmatrices.get(1,i,j);
      Array2D<double> B2 = Bmatrices.get(2,i,j);
      Array2D<double> B3 = Bmatrices.get(3,i,j);
      fs << thegrid.gridNodeLocations().get(i,j) << " " 
         << -A[2][1] << " " << -A[2][2] << " " << -B0[2][1] << " " 
         << -B0[2][2] <<" " << -B0[2][0] << " "
         << -B1[2][0] << " " << -B2[2][0] << " " << -B3[2][0] << endl;
    }
  }
  fs.close();
  cout << "ABcoeffs output to file" << endl <<endl;
  }
 
  //Get all characteristic equation information for each boundary node
  for (int i = 0; i < thegrid.numberElements(); i++){
    CharacteristicFlux left(Amatrices.get(i, 0), trimmedAmatrices.get(i,0));
    CharacteristicFlux right(Amatrices.get(i, thegrid.nodeOrder()),
                             trimmedAmatrices.get(i, thegrid.nodeOrder()));
    AleftBoundaries.push_back(left);
    ArightBoundaries.push_back(right);
  }
}

void DiffEq::setupABmatrices(Grid& thegrid, Modes& lmmodes)
{
  double Omega, Omegap, H, Hp, eL, eLp, fT, fTp, fTpp,rm2M;
  for(int i = 0; i < thegrid.gridNodeLocations().GFvecDim(); i++){
    for(int j = 0; j < thegrid.gridNodeLocations().GFarrDim(); j++){
      //regular wave equation
      if(params.metric.flatspacetime){
        Array2D<double> A(3, 3, 0.0);
        A[1][2] = -pow(params.waveeq.speed, 2.0);
        A[2][1] = -1.0;
        Amatrices.set(i, j, A);
        Array2D<double> tA(2,2,0.0);
        tA[0][1] = -pow(params.waveeq.speed, 2.0);
        tA[1][0] = -1.0;
        trimmedAmatrices.set(i, j, tA);
        //vector dimension of B is actually number of modes
        for(int k = 0; k < Bmatrices.VGFdim(); k++) {
          Array2D<double> B(3, 3, 0.0);
          B[0][1] = -1.0;
          Bmatrices.set(k, i, j, B);
        }
	
       } else if (params.metric.schwarschild) {
        
        int region;
        if (thegrid.gridNodeLocations().get(i,j)==Sminus) {region = 0;}
        else if ((thegrid.gridNodeLocations().get(i,j)>Sminus)&&(thegrid.gridNodeLocations().get(i,j)<Rminus)) {region=1;}
        else if ((thegrid.gridNodeLocations().get(i,j)>=Rminus)&&(thegrid.gridNodeLocations().get(i,j)<=Rplus)) {region=2;}
        else if ((thegrid.gridNodeLocations().get(i,j)>Rplus)&&(thegrid.gridNodeLocations().get(i,j)<Splus)){region=3;}
        else if (thegrid.gridNodeLocations().get(i,j)==Splus) {region=4;}
        
        
        double Omega, Omegap, eL, eLp, H, Hp, term1, term2;
        //thegrid.gridNodeLocations() = rho
        //horizon
        switch (region){
        case 0:
          {
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
	    Array2D<double> tA(2, 2, 0.0);
	    tA[0][1] = -1.0;
	    tA[1][0] = 0.0;
	    tA[1][1] = -1.0;
	    trimmedAmatrices.set(i, j, tA);
	    for(int k = 0; k < lmmodes.ntotal; k++) {
            Array2D<double> B(3, 3, 0.0);
            B[0][2] = -1.0;
            Bmatrices.set(k, i, j, B);
	    }
	    break;
          }
	case 1:
	  {
	    //inner hyperboloidal layer
	    transition(thegrid.gridNodeLocations().get(i, j), Rminus, 
		       Sminus, fT, fTp, fTpp);
	    Omega = 1.0 - thegrid.gridNodeLocations().get(i, j) / Sminus * fT;
	    Omegap = -(fT + thegrid.gridNodeLocations().get(i, j) * fTp) / Sminus;
	    eL = 1.0 + pow(thegrid.gridNodeLocations().get(i, j), 2.0) * fTp 
	      / Sminus;
	    eLp = thegrid.gridNodeLocations().get(i, j) 
	      * (2.0 * fTp + thegrid.gridNodeLocations().get(i, j) * fTpp)
            / Sminus;
	    H = -1.0 + pow(Omega, 2.0) / eL;
	    Hp = (2.0 * Omega * Omegap * eL - pow(Omega, 2.0) * eLp) 
	      / pow(eL, 2.0);
	    thegrid.rstar.set(i, j, thegrid.gridNodeLocations().get(i,j) / Omega);
	    rm2M = invert_tortoise(thegrid.rstar.get(i, j), params.schw.mass);
	    thegrid.rschw.set(i, j, 2.0 * params.schw.mass + rm2M);
	    term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i,j),3.0));
	    term2 = 2.0 * params.schw.mass / thegrid.rschw.get(i,j);
	    Array2D<double> A(3, 3, 0.0);
	    A[1][2] = -1.0;
	    A[2][1] = -(1.0 + H) / (1.0 - H);
	    A[2][2] = 2.0 * H / (1.0 - H);
	    Amatrices.set(i,j,A);
	    Array2D<double> tA(2, 2, 0.0);
	    tA[0][1] = -1.0;
	    tA[1][0] = -(1.0 + H) / (1.0 - H);
	    tA[1][1] = 2.0 * H / (1.0 - H);
	    trimmedAmatrices.set(i,j,tA);
	    for(int k = 0; k < lmmodes.ntotal; k++) {
	      Array2D<double> B(3, 3, 0.0);
	      B[0][2] = -1.0;
	      B[2][1] = -Hp / (1.0 - H);
	      B[2][2] = Hp / (1.0 - H);
	      B[2][0] = 1.0 / (1.0 - pow(H,2.0)) * pow(Omega, 2.0) * term1
               *( lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
	      Bmatrices.set(k, i, j, B);
	    }
	    
	    break;
	  }
         case 2:
           {
	     
           //central tortoise region
	     Omega = 1.0;
	     Omegap = 0.0;
	     eL = 1.0;
	     eLp = 0.0;
	     H = 0.0;
	     Hp = 0.0;
	     thegrid.rstar.set(i, j, thegrid.gridNodeLocations().get(i, j));
	     rm2M = invert_tortoise(thegrid.rstar.get(i, j), params.schw.mass);
	     thegrid.rschw.set(i, j, 2.0 * params.schw.mass + rm2M);
	     term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i, j), 3.0));
	     term2 = 2.0 * params.schw.mass / thegrid.rschw.get(i, j);
	     Array2D<double> A(3, 3, 0.0);
	     A[1][2] = -1.0;
	     A[2][1] = -1.0;
	     Amatrices.set(i, j, A);
	     Array2D<double> tA(2, 2, 0.0);
	     tA[0][1] = -1.0;
	     tA[1][0] = -1.0;
	     trimmedAmatrices.set(i, j, tA);
	     for(int k = 0; k < lmmodes.ntotal; k++) {
	       Array2D<double> B(3, 3, 0.0);
	       B[0][2] = -1.0;
	       B[2][0] = 1.0 / (1.0 - pow(H, 2.0)) * pow(Omega, 2.0) * term1
		 * (lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
	       Bmatrices.set(k, i, j, B);
	     }
	     break;
           }
	case 3:
	  {
	    
	    //outer hyperboloidal region
	    transition(thegrid.gridNodeLocations().get(i,j), Rplus, Splus, fT, fTp, fTpp);
	    Omega = 1.0 - thegrid.gridNodeLocations().get(i, j) / Splus * fT;
	    Omegap = -(fT + thegrid.gridNodeLocations().get(i, j) * fTp) / Splus;
	    eL = 1.0 + pow(thegrid.gridNodeLocations().get(i, j), 2.0) * fTp / Splus; 
	    eLp = thegrid.gridNodeLocations().get(i, j) 
            * (2.0 * fTp + thegrid.gridNodeLocations().get(i, j) * fTpp) / Splus;
	    H = 1.0 - pow(Omega, 2.0) / eL;
	    Hp = -(2.0 * Omega * Omegap * eL - pow(Omega, 2.0) * eLp) 
	      / pow(eL, 2.0);
	    thegrid.rstar.set(i, j, thegrid.gridNodeLocations().get(i, j) / Omega);
	    rm2M = invert_tortoise(thegrid.rstar.get(i, j), params.schw.mass);
	    thegrid.rschw.set(i, j, 2.0 * params.schw.mass + rm2M);
	    term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i, j),3.0));
	    term2 = 2.0 * params.schw.mass / thegrid.rschw.get(i, j);
	    Array2D<double> A(3, 3, 0.0);
	    A[1][2] = -1.0;
	    A[2][1] = -(1.0 - H) / (1.0 + H);
	    A[2][2] = 2.0 * H / (1.0 + H);
	    Amatrices.set(i, j, A);
	    Array2D<double> tA(2, 2, 0.0);
	    tA[0][1] = -1.0;
	    tA[1][0] = -(1.0 - H) / (1.0 + H);
	    tA[1][1] = 2.0 * H / (1.0 + H);
	    trimmedAmatrices.set(i, j, tA);
	    for(int k = 0; k < lmmodes.ntotal; k++) {
	      Array2D<double> B(3, 3, 0.0);
	      B[0][2] = -1.0;
	      B[2][1] = Hp / (1.0 + H);
	      B[2][2] = Hp / (1.0 + H);
	      B[2][0] = 1.0 / (1.0 - pow(H, 2.0)) * pow(Omega, 2.0) * term1
		* (lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
	      Bmatrices.set(k, i, j, B);
	    }
	    break;
             }
        case 4:
          {
	    
	    //scri-plus
	    Omega = 0.0;
	    Omegap = 0.0;
	    eL = 1.0;
	    eLp = 0.0;
	    H = 1.0;
	    Hp = 0.0;
	    thegrid.rstar.set(i, j, -DBL_MAX);
	    thegrid.rschw.set(i,j, -DBL_MAX);
	    term1 = 1.0/ pow(thegrid.gridNodeLocations().get(i,j),2.0);
	    term2 = 0.0;
	    Array2D<double> A(3,3,0.0);
	    A[1][2]=-1.0;
	    A[2][2]=1.0;
	    Amatrices.set(i,j,A);
	    Array2D<double> tA(2,2,0.0);
	    tA[0][1]=-1.0;
	    tA[1][1]=1.0;
	    trimmedAmatrices.set(i,j,tA);
	    for(int k= 0; k < lmmodes.ntotal; k++) {
            Array2D<double> B(3,3,0.0);
            B[0][2]=-1.0;
            B[2][0]=lmmodes.ll[k]*(lmmodes.ll[k]+1.0)/(2.0*pow(Splus,2.0));
            Bmatrices.set(k,i,j,B);
	    }
            break;
          }
	default:
	  {
              throw logic_error("AB matrix region not defined");
	      break;
	  }
        }//end switch case
      }//end inner for            
    }//end outer for
  }//end if schw
}//end function setab


vector<TNT::Array2D<complex<double>>> 
DiffEq::characteristicflux(double t, int modenum, 
                           Grid& thegrid,
                           TwoDVectorGridFunction<complex<double>>& uh, 
                           bool output)
{
  
  //We can loop over characteristicFlux in an external function because
  //in general, neither the RHS of the differential equation nor du mixes
  //spherical harmonic modes. 

  int NumElem = thegrid.numberElements();
  
  vector<Array2D<complex<double>>> du;
  du.resize(NumElem);
#pragma omp parallel for shared(du,NumElem,thegrid,output,modenum,uh) if(uh.TDVGFdim()<=thegrid.numberElements())
  for(int elemnum=0; elemnum<NumElem; elemnum++){
    int indL = 0; //index of leftmost node of that element
    int indR = uh.GFarrDim()-1; //index of rightmost node of that element
    double nL = -1.0; //normal to the leftmost node
    double nR = 1.0; //normal to the rightmost node
    
    //Dimension of the components of the differential equation with 
    //spatial derivatives (dimension of the trimmed A matrices)
    //two for schwarszchild
    int DdimL = AleftBoundaries[elemnum].getDdim();
    int DdimR = ArightBoundaries[elemnum].getDdim();
    
    //vmin and vmax are min and max indices in vector dimension (psi, rho, pi)
    //one and two respectively for schwarzschild
    int vmaxL = AleftBoundaries[elemnum].getAdim() - 1;
    int vminL = vmaxL - DdimL + 1; //neglect zero rows at top of A matrix
    int vmaxR = ArightBoundaries[elemnum].getAdim() - 1;
    int vminR = vmaxR - DdimR + 1; //neglect zero rows at top of A matrix
    
    Array1D<complex<double>> uintL(DdimL); //internal u at left boundary
    Array1D<complex<double>> uintR(DdimR); //internal u at right boundary
    Array1D<complex<double>> uextL(DdimL); //external u at left boundary
    Array1D<complex<double>> uextR(DdimR); //external u at right boundary
    
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
    Array1D<complex<double>> nfluxL1 = matmult(lambdaplusL,
					       matmult(sinvL, uintL));
    Array1D<complex<double>> nfluxL2 =  matmult(lambdaminusL,
						matmult(sinvL, uextL));
    
    Array1D<complex<double>> nfluxL = matmult(SL, nfluxL1+nfluxL2);

 
    Array1D<complex<double>> nfluxR1 = matmult(lambdaplusR,
					       matmult(sinvR, uintR));
    Array1D<complex<double>> nfluxR2 = matmult(lambdaminusR,
					       matmult(sinvR, uextR));
    Array1D<complex<double>> nfluxR = matmult(SR, nfluxR1 + nfluxR2);

    Array2D<double> AtrimmedL= AleftBoundaries[elemnum].getAtrimmed();
    Array2D<double> AtrimmedR= ArightBoundaries[elemnum].getAtrimmed();
    Array2D<complex<double>> duelem(AtrimmedR.dim1(), 2, {0.0,0.0});

    
    //This gets multiplied by lift matrix to calculate flux
    Array1D<complex<double>> temp10 = nL* matmult(AtrimmedL,uintL);
    Array1D<complex<double>> duL = nL * matmult(AtrimmedL, uintL) - nfluxL; 
    Array1D<complex<double>> duR = nR * matmult(AtrimmedR, uintR) - nfluxR; 
    //DU SHOULD BE ZERO AFTER FIRST CALL TO RHS

    
    //    cout << duL[0].real() << " " <<duL[1].real() << " " << duR[0].real() << " " << duR[1].real() << endl;
    
    //    Array1D<complex<double>> fL = matmult(AtrimmedL,uintL);
    //  Array1D<complex<double>> fR = matmult(AtrimmedR,uintR);

    /*    
    if((output)&&(modenum==1)){
      ofstream fs3, fs5, fs6;
      fs3.open("du.txt",ios::app);
      fs3 << setprecision(16);
      fs3 << modenum <<  " " <<  thegrid.gridNodeLocations().get(elemnum, indL)<< " " << duL[0].real() << " " << duL[1].real() << " " <<
        thegrid.gridNodeLocations().get(elemnum, indR)<<" " << duR[0].real() << " " << duR[1].real() << endl;
      fs3.close();
      fs6.open("uL.txt", ios::app);
      fs6 << setprecision(16);
      fs6 << t << " " << modenum << " " << thegrid.gridNodeLocations().get(elemnum, indL) <<  " " << uintL[0].real() << " " << uextL[0].real() << " " << uintL[1].real() << " " << uextL[1].real()<<endl;
      fs6.close();
      fs5.open("uR.txt", ios::app);
      fs5 << setprecision(16);
      fs5 << uintR[0].real() << " " << uextR[0].real() << " " << uintR[1].real() << " " << uextR[1].real()<<endl;
      fs5.close();

    }
    */
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
                 TwoDVectorGridFunction<complex<double>>& uh, 
                 TwoDVectorGridFunction<complex<double>>& RHStdvgf, double t, 
                 vector<Array2D<complex<double>>>& du , bool output)
{

  //We can loop over RHS in an external function independent of mode 
  //because in general the right hand side of the differential equation
  //does not mix spherical harmonic modes. Neither does the numerical flux. 
   

  int NumElem = thegrid.numberElements();
  int vmaxAB = ArightBoundaries[0].getAdim()-1;
  int vminA = vmaxAB - ArightBoundaries[0].getDdim() +1;
  Array2D<double> D;
  D = thegrid.refelem.getD();
  Array2D<double> Lift;
  Lift = thegrid.refelem.getLift();
  
  //pragma omp parallel if(NumElem>lmmodes.ntotal) for

#pragma omp parallel for shared(uh,modenum,thegrid,SOURCE_VECNUM,RHStdvgf,params) if(uh.TDVGFdim()<=thegrid.numberElements())
  for(int elemnum = 0; elemnum < NumElem; elemnum++){
    //Maximum index for both A and B matrix
    //    int vmaxAB = ArightBoundaries[elemnum].getAdim() - 1;
    //Minimum index for use with trimmed A matrix. 
    //Minimum index for B matrix is zero
    //    int vminA = vmaxAB - ArightBoundaries[elemnum].getDdim() + 1;
    
    //The B matrix component of the RHS. 
    Array2D<complex<double>> RHSB(uh.GFarrDim(), 
                                  ArightBoundaries[elemnum].getAdim());
    for(int nodenum = 0; nodenum < uh.GFarrDim(); nodenum++){
      Array1D<complex<double>> RHSBpernode;
      
      //FIX THIS. THIS CAN BE SPED UP BY NOT COPYING IN GETVECTORASARRAY1D
      //Multiply the B matrix times a "vector" of u for each node
      RHSBpernode = matmult(Bmatrices.get(modenum, elemnum, nodenum),
                            uh.getVectorAsArray1D(modenum,elemnum, nodenum, 
                                                  0, vmaxAB, 0));

      //Insert that result into the rows of a larger matrix
      insert_1D_into_2D(RHSB, RHSBpernode, nodenum, false);
      

    }//This can be sped up by skipping the insert step and reading directly 
    //from per node


    //NOT CHECKED
    /*    if((output)&&(modenum==1)){
      cout << RHSB[nodenum][0] << " " << RHSB[nodenum][1] << " " <<  RHSB[nodenum][2] << endl;
      }*/
    
      

    //A contribution:
    Array2D<complex<double>> RHSA1(uh.GFarrDim(), ArightBoundaries[elemnum].getDdim());
    
    //The A contribution needs to be multiplied one node at a time by the
    //trimmed A matrix in a similar manner to the B contribution. But first,
    //we take the spatial derivative across all thegrid.gridNodeLocations(). 


    //FIX THIS. THIS CAN BE SPED UP BY NOT COPYING IN GETVECTORNODEARRAY
    //    Array2D<complex<double>> RHSA1preA = move(thegrid.jacobian(elemnum) 
    //					      * move((matmult(thegrid.refelem.getD(),
    //							      uh.getVectorNodeArray2D(modenum, elemnum, vminA, vmaxAB, 0)))));

    Array2D<complex<double>> RHSA1preA = thegrid.jacobian(elemnum) 
					      * matmult(D,
						 uh.getVectorNodeArray2D(modenum, elemnum, vminA, vmaxAB, 0));

    //RHSA1preA is du2drho and du2dpi.
    /*    if(modenum==1){
      for(int nodenum=0; nodenum < uh.GFarrDim(); nodenum++){ 
	cout << setprecision(15);
	cout << t << " " << thegrid.gridNodeLocations().get(elemnum, nodenum) << " " << RHSA1preA[nodenum][0].real() << " " << RHSA1preA[nodenum][1].real() << endl;
      }
      }*/
    //Multiply each row of RHSA1preA by a different a Atrimmed matrix
    for(int nodenum=0; nodenum < uh.GFarrDim(); nodenum++)
      {
        int M = trimmedAmatrices.get(elemnum,nodenum).dim1();
        int N = trimmedAmatrices.get(elemnum,nodenum).dim2();
        int K = RHSA1preA.dim1();
        int L = RHSA1preA.dim2();
        
	//        Array2D<double> tA(M,N); 
	// tA= 
	
        for (int i=0; i<M; i++){
          complex<double> sum = 0;
	  for (int k=0; k<N; k++)
	    sum += -trimmedAmatrices.get(elemnum, nodenum)[i][k]
	      * RHSA1preA[nodenum][k];
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

    //get rid of moves?? PARALLEL PROBLEM?
    //This is the contribution due to du, or the numerical flux
    //    Array2D<complex<double>> RHSA2 = move(thegrid.jacobian(elemnum) 
    //					  *move(matmult(thegrid.refelem.getLift(), du[elemnum])));
    Array2D<complex<double>> RHSA2 = thegrid.jacobian(elemnum) 
      *matmult(Lift, du[elemnum]);

      
    //RHSA and RHSB will have different sizes due to the different 
    //number of diffeq variables stored in each. sum them using a 
    //for loop while assigning values to the RHStdvgf vector grid function
    
    //    Array2D<complex<double>> RHSA = RHSA1 + RHSA2;
    
      //Sum the contributions from B, derivative, and flux, 
      //accounting for different matrix dimensions
    for(int vecnum = 0; vecnum < RHStdvgf.VGFdim(); vecnum++){
      for(int nodenum = 0; nodenum < RHStdvgf.GFarrDim(); nodenum++){
	complex<double> tot;
	if(vecnum<vminA){
	  tot = -RHSB[nodenum][vecnum];
	} else {
	  /*if((output)&&(modenum==1)&&(vecnum==2)){
	    cout << setprecision(15);
	    cout << nodenum << " " << RHSA[nodenum][vecnum-vminA].real() << endl;
	    }*/
	  tot=-RHSB[nodenum][vecnum]+RHSA1[nodenum][vecnum-vminA]
	    + RHSA2[nodenum][vecnum-vminA];
	}//-sign in B because it is on the left hand side of the 
	if((params.opts.useSource)&&(vecnum==SOURCE_VECNUM)){
	  tot += source.get(modenum, elemnum, nodenum);
                       
	  //equation in the definition supplied in DiffEq.cpp
	}
	RHStdvgf.set(modenum, vecnum, elemnum, nodenum, tot);
      }
    }
  }
}

void DiffEq::modeRHS(Grid& thegrid,
                     TwoDVectorGridFunction<complex<double>>& uh,
                     TwoDVectorGridFunction<complex<double>>& RHStdvgf,
                     double t, bool output)
{
  if(params.opts.useSource) {
    fill_source_all(thegrid, t, uh.TDVGFdim(), source, window,
		    dwindow, d2window);
  }
  if (output&&params.opts.useSource){
    for (int k = 0; k<source.VGFdim(); k++){
      ofstream fs;
      fs.precision(16);
      ostringstream oss;
      oss << "source" << "." << k << ".txt";
      fs.open(oss.str(),ios::app);
      fs << endl << endl;
      fs << " #time = " << t << endl;
      
      
      for (int i = 0; i < source.GFvecDim(); i++){
	for(int j = 0; j < source.GFarrDim(); j++){
	  //Print out at select time steps
	  fs << thegrid.gridNodeLocations().get(i, j) << " "
	     << source.get(k, i, j).real() << " " 
	     << source.get(k, i, j).imag() <<endl;
	}
      }
    }
  }

  

#pragma omp parallel for if(uh.TDVGFdim()>thegrid.numberElements())
  for(int modenum = 0; modenum < uh.TDVGFdim(); modenum++) {
    //    cout << modenum << endl;
    double max_speed = 1.0;
    vector<Array2D<complex<double>>> du;
    du = characteristicflux(t, modenum, thegrid, uh, output);
    RHS(modenum, thegrid, uh, RHStdvgf, t, du, output);
  }
}

