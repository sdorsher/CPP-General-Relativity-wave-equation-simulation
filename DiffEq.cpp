#include "DiffEq.h"

/*
Fla spacetime wave equation:
drho/dt = c^2 dpi/dx
dpi/dt=drho/dx
drho/dt=psi


also dpi/dx=psi, but this is not needed for evolution
*/


vector<double> DiffEq::getA(int gridindex, int pointsindex)
{
  return Amatrices.get(gridindex, pointsindex);
}

vector<double> DiffEq::getB(int modesindex, int gridindex, int pointsindex)
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
vector<double> DiffEq::getAtrimmed(int gridindex, int pointsindex)
{
  return trimmedAmatrices.getVector(gridindex, pointsindex);
}

// The A matrix must be formatted such that zero rows are at the top

DiffEq::DiffEq(Grid& thegrid, Modes& lmmodes, int nmodetotal, Coordinates& coordobj):
  Amatrices{params.grid.Adim*params.grid.Adim,params.grid.numelems, params.grid.elemorder + 1,0.},
  Bmatrices{nmodetotal,params.grid.Adim*params.grid.Adim, params.grid.numelems, 
	    params.grid.elemorder + 1,0.},
  trimmedAmatrices{params.grid.Ddim*params.grid.Ddim,params.grid.numelems, params.grid.elemorder + 1,0.},
  source{nmodetotal, params.grid.numelems, params.grid.elemorder+1,{0.0,0.0}},
  uextL{params.grid.Ddim},uextR{params.grid.Ddim},uintL{params.grid.Ddim},
  uintR{params.grid.Ddim}
  {
    //set up the A and B matrices


    printf("setting up A and B matrices\n");
    setupABmatrices(thegrid, lmmodes,coordobj);
    printf("A and B matrices established\n");
    
    //Get all characteristic equation information for each boundary node
    for (int i = 0; i < thegrid.numberElements(); i++){

      CharacteristicFlux left(Amatrices.getVector(i, 0), trimmedAmatrices.getVector(i,0));

      CharacteristicFlux right(Amatrices.getVector(i, thegrid.nodeOrder()),
			       trimmedAmatrices.getVector(i, thegrid.nodeOrder()));
      AleftBoundaries.push_back(left);
      ArightBoundaries.push_back(right);
    }
  }

void DiffEq::setupABmatrices(Grid& thegrid, Modes& lmmodes, Coordinates& coordobj)
{
  double Omega, Omegap, H, Hp, eL, eLp, fT, fTp, fTpp,rm2M;
  for(int i = 0; i < thegrid.gridNodeLocations().GFvecDim(); i++){
    for(int j = 0; j < thegrid.gridNodeLocations().GFarrDim(); j++){
      //regular wave equation
      if(params.metric.flatspacetime){
        //Array2D<double> A(3, 3, 0.0);
        //A[1][2] = -pow(params.waveeq.speed, 2.0);
        //A[2][1] = -1.0;
        //Amatrices.set(i, j, A);
	Amatrices.set(1*params.grid.Adim+2,i,j,-pow(params.waveeq.speed, 2.0));
	Amatrices.set(2*params.grid.Adim+1,i,j,-1.0);
        //Array2D<double> tA(2,2,0.0);
        //tA[0][1] = -pow(params.waveeq.speed, 2.0);
        //tA[1][0] = -1.0;
        //trimmedAmatrices.set(i, j, tA);

	trimmedAmatrices.set(0*params.grid.Ddim+1,i,j,-pow(params.waveeq.speed,2.0));
	trimmedAmatrices.set(1*params.grid.Ddim+0,i,j,-1.0);
	
        //TDVGF dimension of B is actually number of modes
        for(int k = 0; k < Bmatrices.TDVGFdim(); k++) {
          //Array2D<double> B(3, 3, 0.0);
          //B[0][1] = -1.0;
          //Bmatrices.set(k, i, j, B);
	  Bmatrices.set(k,0*params.grid.Adim+1,i,j,-1.0);
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
	    thegrid.rstar.set(i,j,-DBL_MAX);
	    thegrid.rschw.set(i,j,2.0*params.schw.mass); 
	    term1 = 0.0;
	    term2 = 1.0;
	    //Array2D<double> A(3, 3, 0.0);
	    //A[1][2] = -1.0;
	    //A[2][1] = 0.0;
	    //A[2][2] = -1.0;
	    //Amatrices.set(i, j, A);
	    Amatrices.set(1*params.grid.Adim+2,i,j,-1.0);
	    Amatrices.set(2*params.grid.Adim+2,i,j,-1.0);
	    //Array2D<double> tA(2, 2, 0.0);
	    //tA[0][1] = -1.0;
	    //tA[1][0] = 0.0;
	    //tA[1][1] = -1.0;
	    //trimmedAmatrices.set(i, j, tA);
	    trimmedAmatrices.set(0*params.grid.Ddim+1,i,j,-1.0);
	    trimmedAmatrices.set(1*params.grid.Ddim+1,i,j,-1.0);
	    for(int k = 0; k < lmmodes.ntotal; k++) {
	      //Array2D<double> B(3, 3, 0.0);
	      //B[0][2] = -1.0;
	      //Bmatrices.set(k, i, j, B);
	      Bmatrices.set(k,0*params.grid.Adim+2,i,j,-1.0);
	    }
	    break;
          }
	case 1:
	  {
	    //inner hyperboloidal layer
	    coordobj.transition(thegrid.gridNodeLocations().get(i, j), Rminus, 
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
	    rm2M = coordobj.invert_tortoise(thegrid.rstar.get(i, j), params.schw.mass);
	    thegrid.rschw.set(i, j, 2.0 * params.schw.mass + rm2M);
	    term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i,j),3.0));
	    term2 = 2.0 * params.schw.mass / thegrid.rschw.get(i,j);
	    //Array2D<double> A(3, 3, 0.0);
	    //A[1][2] = -1.0;
	    //A[2][1] = -(1.0 + H) / (1.0 - H);
	    //A[2][2] = 2.0 * H / (1.0 - H);
	    Amatrices.set(1*params.grid.Adim+2,i,j,-1.0);
	    Amatrices.set(2*params.grid.Adim+1,i,j,-(1.0 + H) / (1.0 - H));
	    Amatrices.set(2*params.grid.Adim+2,i,j,2.0 * H / (1.0 - H));
	    //Amatrices.set(i,j,A);
	    //Array2D<double> tA(2, 2, 0.0);
	    //tA[0][1] = -1.0;
	    //tA[1][0] = -(1.0 + H) / (1.0 - H);
	    //tA[1][1] = 2.0 * H / (1.0 - H);
	    //trimmedAmatrices.set(i,j,tA);
	    trimmedAmatrices.set(0*params.grid.Ddim+1,i,j,-1.0);
	    trimmedAmatrices.set(1*params.grid.Ddim+0,i,j, -(1.0 + H) / (1.0 - H));
	    trimmedAmatrices.set(1*params.grid.Ddim+1,i,j,2.0 * H / (1.0 - H));
	    for(int k = 0; k < lmmodes.ntotal; k++) {
	      //Array2D<double> B(3, 3, 0.0);
	      //B[0][2] = -1.0;
	      //B[2][1] = -Hp / (1.0 - H);
	      //B[2][2] = Hp / (1.0 - H);
	      //B[2][0] = 1.0 / (1.0 - pow(H,2.0)) * pow(Omega, 2.0) * term1
              // *( lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
	      //Bmatrices.set(k, i, j, B);
	      Bmatrices.set(k,0*params.grid.Adim+2,i,j,-1.0);
	      Bmatrices.set(k,2*params.grid.Adim+0,i,j, term1/(1.0-pow(H,2.0))*(term2+lmmodes.ll[k]*(lmmodes.ll[k]+1)));
	      Bmatrices.set(k,2*params.grid.Adim+1,i,j,-Hp / (1.0 - H));
	      Bmatrices.set(k,2*params.grid.Adim+2,i,j,Hp / (1.0 - H));	      double temp=1.0 / (1.0 - pow(H,2.0)) * pow(Omega, 2.0) * term1*( lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);

	      //	      cout << i << " " << j << " " << k << " " << temp << endl;
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
	     rm2M = coordobj.invert_tortoise(thegrid.rstar.get(i, j), params.schw.mass);
	     thegrid.rschw.set(i, j, 2.0 * params.schw.mass + rm2M);
	     term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i, j), 3.0));
	     term2 = 2.0 * params.schw.mass / thegrid.rschw.get(i, j);
	     //Array2D<double> A(3, 3, 0.0);
	     //A[1][2] = -1.0;
	     //A[2][1] = -1.0;
	     //Amatrices.set(i, j, A);
	     Amatrices.set(1*params.grid.Adim+2,i,j,-1.0);
	     Amatrices.set(2*params.grid.Adim+1,i,j,-1.0);
	     //Array2D<double> tA(2, 2, 0.0);
	     //tA[0][1] = -1.0;
	     //tA[1][0] = -1.0;
	     //trimmedAmatrices.set(i, j, tA);
	     trimmedAmatrices.set(0*params.grid.Ddim+1,i,j,-1.0);
	     trimmedAmatrices.set(1*params.grid.Ddim+0,i,j,-1.0);
	     for(int k = 0; k < lmmodes.ntotal; k++) {
	       //Array2D<double> B(3, 3, 0.0);
	       //B[0][2] = -1.0;
	       //B[2][0] = 1.0 / (1.0 - pow(H, 2.0)) * pow(Omega, 2.0) * term1
	       // * (lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
	       //Bmatrices.set(k, i, j, B);
	       Bmatrices.set(k,0*params.grid.Adim+2,i,j,-1.0);
	       double temp =term1 * (lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
	       Bmatrices.set(k,2*params.grid.Adim+0,i,j,temp);
	       //cout << i << " " << j << " " << k <<" " << temp << endl;
	     }
	     break;
           }
	case 3:
	  {
	    
	    //outer hyperboloidal region
	    coordobj.transition(thegrid.gridNodeLocations().get(i,j), Rplus, Splus, fT, fTp, fTpp);
	    Omega = 1.0 - thegrid.gridNodeLocations().get(i, j) / Splus * fT;
	    Omegap = -(fT + thegrid.gridNodeLocations().get(i, j) * fTp) / Splus;
	    eL = 1.0 + pow(thegrid.gridNodeLocations().get(i, j), 2.0) * fTp / Splus; 
	    eLp = thegrid.gridNodeLocations().get(i, j) 
            * (2.0 * fTp + thegrid.gridNodeLocations().get(i, j) * fTpp) / Splus;
	    H = 1.0 - pow(Omega, 2.0) / eL;
	    Hp = -(2.0 * Omega * Omegap * eL - pow(Omega, 2.0) * eLp) 
	      / pow(eL, 2.0);
	    thegrid.rstar.set(i, j, thegrid.gridNodeLocations().get(i, j) / Omega);
	    rm2M = coordobj.invert_tortoise(thegrid.rstar.get(i, j), params.schw.mass);
	    thegrid.rschw.set(i, j, 2.0 * params.schw.mass + rm2M);
	    term1 = rm2M / (pow(Omega, 2.0) * pow(thegrid.rschw.get(i, j),3.0));
	    term2 = 2.0 * params.schw.mass / thegrid.rschw.get(i, j);
	    //Array2D<double> A(3, 3, 0.0);
	    //A[1][2] = -1.0;
	    //A[2][1] = -(1.0 - H) / (1.0 + H);
	    //A[2][2] = 2.0 * H / (1.0 + H);
	    //Amatrices.set(i, j, A);
	    Amatrices.set(1*params.grid.Adim+2,i,j,-1.0);
	    Amatrices.set(2*params.grid.Adim+1,i,j,-(1.0 - H) / (1.0 + H));
	    Amatrices.set(2*params.grid.Adim+2,i,j,2.0 * H / (1.0 + H));
	    
	    //Array2D<double> tA(2, 2, 0.0);
	    //tA[0][1] = -1.0;
	    //tA[1][0] = -(1.0 - H) / (1.0 + H);
	    //tA[1][1] = 2.0 * H / (1.0 + H);
	    //trimmedAmatrices.set(i, j, tA);

	    trimmedAmatrices.set(0*params.grid.Ddim+1,i,j,-1.0);
	    trimmedAmatrices.set(1*params.grid.Ddim+0,i,j,-(1.0 - H) / (1.0 + H));
	    trimmedAmatrices.set(1*params.grid.Ddim+1,i,j,2.0 * H / (1.0 + H));
	    for(int k = 0; k < lmmodes.ntotal; k++) {
	      //Array2D<double> B(3, 3, 0.0);
	      //B[0][2] = -1.0;
	      //B[2][1] = Hp / (1.0 + H);
	      //B[2][2] = Hp / (1.0 + H);
	      //B[2][0] = 1.0 / (1.0 - pow(H, 2.0)) * pow(Omega, 2.0) * term1
		//* (lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
	      //Bmatrices.set(k, i, j, B);
	      Bmatrices.set(k,0*params.grid.Adim+2,i,j,-1.0);
	      Bmatrices.set(k,2*params.grid.Adim+0,i,j,term1/(1.0-pow(H,2.0))*(lmmodes.ll[k]*(lmmodes.ll[k]+1)+term2));
	      Bmatrices.set(k,2*params.grid.Adim+1,i,j,Hp / (1.0 + H));
	      Bmatrices.set(k,2*params.grid.Adim+2,i,j,Hp / (1.0 + H));
	      double temp = 1.0 / (1.0 - pow(H, 2.0)) * pow(Omega, 2.0) * term1* (lmmodes.ll[k] * (lmmodes.ll[k] + 1.0) + term2);
	      //Bmatrices.set(k,2*params.grid.Adim+0,i,j,temp);
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
	    thegrid.rstar.set(i, j, DBL_MAX);
	    thegrid.rschw.set(i,j, DBL_MAX);
	    term1 = 1.0/ pow(thegrid.gridNodeLocations().get(i,j),2.0);
	    term2 = 0.0;
	    //Array2D<double> A(3,3,0.0);
	    //A[1][2]=-1.0;
	    //A[2][2]=1.0;
	    //Amatrices.set(i,j,A);
	    Amatrices.set(1*params.grid.Adim+2,i,j,-1.0);
	    Amatrices.set(2*params.grid.Adim+2,i,j,1.0);
	    //Array2D<double> tA(2,2,0.0);
	    //tA[0][1]=-1.0;
	    //tA[1][1]=1.0;
	    //trimmedAmatrices.set(i,j,tA);
	    trimmedAmatrices.set(0*params.grid.Ddim+1,i,j,-1.0);
	    trimmedAmatrices.set(1*params.grid.Ddim+1,i,j,1.0);
	    for(int k= 0; k < lmmodes.ntotal; k++) {
	      //Array2D<double> B(3,3,0.0);
            //B[0][2]=-1.0;
            //B[2][0]=lmmodes.ll[k]*(lmmodes.ll[k]+1.0)/(2.0*pow(Splus,2.0));
            //Bmatrices.set(k,i,j,B);
	      Bmatrices.set(k,0*params.grid.Adim+2,i,j,-1.0);

	      //FIX ME: THIS DOESN"T LOOK RIGHT
	      double temp =lmmodes.ll[k]*(lmmodes.ll[k]+1.0)/(2.0*pow(Splus,2.0));
	      Bmatrices.set(k,2*params.grid.Adim+0,i,j,temp);
			  
	      
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

void DiffEq::set_coefficients(Grid &thegrid, EllipticalOrbit* orb, Coordinates &coords, double& maxspeed, int elemnum, Modes& lmmodes, double xp, double xip, double dxpdt, double d2xpdt2, vector<double>& dxdt, vector<double>& dxdxi)


{
  double ne = params.grid.elemorder;

  vector<double> rm2m(ne+1), x(ne+1), d2xdt2(ne+1), d2xdxi2(ne+1), d2xdtdxi(ne+1);

  coords.coord_trans(Rminus, Rplus, thegrid, xp, xip, dxpdt, d2xpdt2, x, dxdt, dxdxi, d2xdt2, d2xdxi2, d2xdtdxi, elemnum);


  for(int i = 0; i<=ne; i++){
    double dxdxiinv = 1.0/dxdxi.at(i);
    double dxdxiinv2 = dxdxiinv*dxdxiinv;
    double dxdxiinv3 = dxdxiinv2*dxdxiinv;
    double dxdt2 = dxdt.at(i)*dxdt.at(i);
    
    double coeff1 = -(1.-dxdt2)*dxdxiinv2;
    double coeff2 = -2.*dxdt.at(i)*dxdxiinv;
    double coeff3 = -(d2xdxi2.at(i)*(dxdt2-1.)
		      +dxdxi.at(i)*(-2.*dxdt.at(i)*d2xdtdxi.at(i)
				    +dxdxi.at(i)*d2xdt2.at(i)))*dxdxiinv3;
    double coeff4 = 0.0;


    cout <<  elemnum << " " << coeff1 << " " << coeff2 << " " << coeff3 << endl;
    
    Amatrices.set(1*params.grid.Adim+2,elemnum,i,coeff1);
    Amatrices.set(2*params.grid.Adim+2,elemnum,i,coeff2);
    for(int modenum=0; modenum< lmmodes.ntotal; modenum++){
      Bmatrices.set(modenum,1*params.grid.Adim+2,elemnum,i,coeff3);
      Bmatrices.set(modenum,2*params.grid.Adim+2,elemnum,i,coeff4);
    }
    trimmedAmatrices.set(0*params.grid.Ddim+1,elemnum,i,coeff1);
    trimmedAmatrices.set(1*params.grid.Ddim+1,elemnum,i,coeff2);

    int boundary;
    
    if(i==0){
      AleftBoundaries[elemnum].LambdaV.at(0*params.grid.Ddim+0)=-(1.0+dxdt.at(i))*dxdxiinv;
      AleftBoundaries[elemnum].LambdaV.at(1*params.grid.Ddim+1)=(1.-dxdt.at(i))*dxdxiinv;
      maxspeed=max(max(abs(AleftBoundaries[elemnum].LambdaV.at(0*params.grid.Ddim+0)),abs(AleftBoundaries[elemnum].LambdaV.at(1*params.grid.Ddim+1))),maxspeed);
      AleftBoundaries[elemnum].SmatrixV.at(0*params.grid.Ddim+0)=(1.+dxdt.at(i))*dxdxiinv;
      AleftBoundaries[elemnum].SmatrixV.at(0*params.grid.Ddim+1)=(1.0+dxdt.at(i))*dxdxiinv;
      AleftBoundaries[elemnum].SmatrixV.at(1*params.grid.Ddim+0)=(-1.0+dxdt.at(i))*dxdxiinv;
      AleftBoundaries[elemnum].SmatrixV.at(0*params.grid.Ddim+1)=1.0;
      AleftBoundaries[elemnum].SmatrixV.at(1*params.grid.Ddim+1)=1.0;
      AleftBoundaries[elemnum].SinvV.at(0*params.grid.Ddim+0)=0.5*dxdxi.at(i);
      AleftBoundaries[elemnum].SinvV.at(0*params.grid.Ddim+1)=0.5*(1.0-dxdt.at(i));
      AleftBoundaries[elemnum].SinvV.at(1*params.grid.Ddim+0)=-0.5*dxdxi.at(i);
      AleftBoundaries[elemnum].SinvV.at(1*params.grid.Ddim+1)=0.5*(1.+dxdt.at(i));
    } else if(i==params.grid.elemorder){
      ArightBoundaries[elemnum].LambdaV.at(0*params.grid.Ddim+0)=-(1.0+dxdt.at(i))*dxdxiinv;
      ArightBoundaries[elemnum].LambdaV.at(1*params.grid.Ddim+1)=(1.-dxdt.at(i))*dxdxiinv;
      maxspeed=max(max(abs(ArightBoundaries[elemnum].LambdaV.at(0*params.grid.Ddim+0)),abs(ArightBoundaries[elemnum].LambdaV.at(1*params.grid.Ddim+1))),maxspeed);
      ArightBoundaries[elemnum].SmatrixV.at(0*params.grid.Ddim+0)=(1.+dxdt.at(i))*dxdxiinv;
      ArightBoundaries[elemnum].SmatrixV.at(0*params.grid.Ddim+1)=(1.0+dxdt.at(i))*dxdxiinv;
      ArightBoundaries[elemnum].SmatrixV.at(1*params.grid.Ddim+0)=(-1.0+dxdt.at(i))*dxdxiinv;
      ArightBoundaries[elemnum].SmatrixV.at(0*params.grid.Ddim+1)=1.0;
      ArightBoundaries[elemnum].SmatrixV.at(1*params.grid.Ddim+1)=1.0;
      ArightBoundaries[elemnum].SinvV.at(0*params.grid.Ddim+0)=0.5*dxdxi.at(i);
      ArightBoundaries[elemnum].SinvV.at(0*params.grid.Ddim+1)=0.5*(1.0-dxdt.at(i));
      ArightBoundaries[elemnum].SinvV.at(1*params.grid.Ddim+0)=-0.5*dxdxi.at(i);
      ArightBoundaries[elemnum].SinvV.at(1*params.grid.Ddim+1)=0.5*(1.+dxdt.at(i));

    }
  }

  double mass = params.schw.mass;
  
  for(int i=0; i<=params.grid.elemorder; i++){

    rm2m.at(i)= coords.invert_tortoise(thegrid.rstar.get(elemnum,i),mass);
  }

  for(int i=0; i<=params.grid.elemorder; i++){
    thegrid.rschw.set(elemnum,i,2.0*mass+rm2m.at(i));
    double rfac=rm2m.at(i)/pow(thegrid.rschw.get(elemnum,i),4.);
    for(int k=0; k<lmmodes.ntotal; k++){
      double coeffl = (lmmodes.ll.at(k)*(lmmodes.ll.at(k)+1)*thegrid.rschw.get(elemnum,i)+2.*mass)*rfac;
      Bmatrices.set(k,1*params.grid.Adim+2,elemnum,i,coeffl);
    }
  }

  //finds right side of element boundary in time dependent region. for dxdt and dxdxi
  coords.dxdxibL0.at(elemnum)=dxdxi.at(0);
  coords.dxdxibL1.at(elemnum)=dxdt.at(0);
  coords.dxdxibR0.at(elemnum)=dxdxi.at(params.grid.elemorder);
  coords.dxdxibR1.at(elemnum)=dxdt.at(params.grid.elemorder);
  //cout << "coordsdx : " << elemnum << " "<< coords.dxdxibL0.at(elemnum) <<" " << coords.dxdxibL1.at(elemnum) <<" " << coords.dxdxibR0.at(elemnum) <<" " << coords.dxdxibL1.at(elemnum) << endl; 

  // cout << setprecision(15);
  //cout << elemnum <<  " " << dxdt.at(0) << " " << dxdxi.at(0) <<  " " << dxdt.at(params.grid.elemorder) << " " <<  dxdxi.at(params.grid.elemorder) << endl;

  //GOOD THROUGH HERE
  if(abs(thegrid.gridNodeLocations().get(elemnum,0)-coords.xip)<1e-10){
    orb->drdlambda_particle=dxdt.at(0);
    orb->drdxi_particle=dxdxi.at(0);
  }
}

  
void DiffEq::RHS(int modenum, Grid& thegrid,
                 TwoDVectorGridFunction<complex<double>>& uh, 
                 TwoDVectorGridFunction<complex<double>>& RHStdvgf, double t, bool output, Coordinates& coords, WorldTube* wt)
{


  
  int NumElem = thegrid.numberElements();
  int vmaxAB = params.grid.Adim-1;
  int vminA = vmaxAB - params.grid.Ddim +1;
  vector<double> Dmatrix; //was Array2D<double>
  Dmatrix = thegrid.refelem.getD();
  vector<double> Lift; //was Array2D<double>
  Lift = thegrid.refelem.getLift();

  //We can loop over RHS in an external function independent of mode 
  //because in general the right hand side of the differential equation
  //does not mix spherical harmonic modes. Neither does the numerical flux. 
  
  //int NumElem = thegrid.numberElements();
  
  //du.resize(NumElem);
  //#pragma omp parallel for shared(du,NumElem,thegrid,output,modenum,uh) if(uh.TDVGFdim()<=thegrid.numberElements())

  for(int elemnum=0; elemnum<NumElem; elemnum++){
     vector<complex<double>> du; //inner was Array2D
    du.resize(2*params.grid.Ddim);
    int indL = 0; //index of leftmost node of that element
    int indR = uh.GFarrDim()-1; //index of rightmost node of that element
    double nL = -1.0; //normal to the leftmost node
    double nR = 1.0; //normal to the rightmost node
    
    
    
    //Dimension of the components of the differential equation with 
    //spatial derivatives (dimension of the trimmed A matrices)
    //two for schwarszchild
        
    //vmin and vmax are min and max indices in vector dimension (psi, rho, pi)
    //one and two respectively for schwarzschild
    int vmaxL = params.grid.Adim - 1;
    int vminL = vmaxL - params.grid.Ddim + 1; //neglect zero rows at top of A matrix
    int vmaxR = params.grid.Adim - 1;
    int vminR = vmaxR - params.grid.Ddim + 1; //neglect zero rows at top of A matrix

    double sstre, sstim, ssrre, ssrim, alpha;
 
    int nodenumFound;
    bool bleft,bright,add, sub, found;
    //although subtraction and adding from left and right do not occur at same element number, always want to consider them at current element number to make sure to iterate over all of them

    
    if(params.opts.use_world_tube){
      if(elemnum==0){
	//cout << "zero" << endl;// FIX ME
	bright=false;
	bleft = wt->addSingFieldToLeftElemExt.at(elemnum)||wt->subSingFieldFromLeftElemExt.at(elemnum);
	add = wt->addSingFieldToLeftElemExt.at(elemnum);
	sub= wt->subSingFieldFromLeftElemExt.at(elemnum);
      } else if (elemnum==NumElem-1){
	//cout << "max" << endl;//FIX ME
        bleft = false;
	bright = false;
	add =false; //wt->addSingFieldToLeftElemExt.at(elemnum);
	sub = false;//wt->subSingFieldFromLeftElemExt.at(elemnum);
      }else {
	//right or left side of boundary (j or j+1 in WorldTube.cpp)
	bright = wt->subSingFieldFromRightElemExt.at(elemnum)||wt->addSingFieldToRightElemExt.at(elemnum);
	bleft = wt->addSingFieldToLeftElemExt.at(elemnum)||wt->subSingFieldFromLeftElemExt.at(elemnum);
	
	add = wt->addSingFieldToLeftElemExt.at(elemnum)||wt->addSingFieldToRightElemExt.at(elemnum); // add ext sing field to this eleement
	sub=wt->subSingFieldFromRightElemExt.at(elemnum)||wt->subSingFieldFromLeftElemExt.at(elemnum); //subtract external singular field from this element
	//cout << elemnum << " " << (bleft ? 1:0) << " " << (bright ? 1 : 0) << " " << (add ? 1:0) << " "<< (sub ? 1:0) << endl;
      }
      if(bleft){
	nodenumFound=params.grid.elemorder;
      }else if(bright){
	nodenumFound=0;
      }
      
      double ssign;
      if(add){
	ssign=1.;
      }else if(sub){
	ssign = -1.;
      }


      found = bleft || bright;

      if(found){

	double node = thegrid.rschw.get(elemnum, nodenumFound);
	dPhi_dt(&modenum, &node, &sstre, &sstim);
	dPhi_dr(&modenum, &node, &ssrre, &ssrim);
	alpha = (thegrid.rschw.get(elemnum,nodenumFound)-2.*params.schw.mass)
	  /thegrid.rschw.get(elemnum,nodenumFound);

	sstre = ssign*sstre;
	sstim = ssign*sstim;
	ssrre = alpha *ssign*ssrre;
	ssrim = alpha *ssign*ssrim;
	//	if((modenum==1&found){
	//cout << setprecision(15);
	//cout << modenum  << " " << elemnum << " " << nodenumFound<< " " << " "<< ssign << " " << alpha << " " << t << " " << thegrid.gridNodeLocations().get(elemnum,nodenumFound) << " " << sstre <<  " " << sstim << " " << ssrre << " " << ssrim <<  " " << (add ? 1 : 0 ) << " " << ( sub ? 1 : 0) << " " << (bleft ? 1:0) << " " << (right ? 1:0)<< endl;
	//}
      }
    }

    uintL = uh.getVectorRange(modenum,elemnum, indL, vminL, vmaxL, 0);
    uintR = uh.getVectorRange(modenum,elemnum, indR, vminR, vmaxR, 0);


    if(elemnum > 0) {
      uextL = uh.getVectorRange(modenum, elemnum - 1, indR, vminL,
				vmaxL, 0);
      //external u, left boundary
    }else{
      uextL = uh.getVectorRange(modenum, NumElem - 1, indR, vminL,
				vmaxL, 0);
      //periodic boundary conditions
    }


    if(elemnum < NumElem - 1) {
      uextR = uh.getVectorRange(modenum, elemnum + 1, indL, vminR,
				vmaxR, 0);
      //external u, right boundary
    }else{
      uextR = uh.getVectorRange(modenum, 0, indL, vminR, vmaxR, 0);
      //periodic boundary conditions
    }

    //     if((elemnum==3)&&(output))
    //{
    //cout << "LHS0 ext=" << uextL[0] << " uint=" << uintL[0] << endl;
    //cout << "RHS0 uext=" << uextR[0] << " uint=" << uintR[0] << endl;
    //cout << "LHS1 ext=" << uextL[1] << " uint=" << uintL[1] << endl;
    //cout << "RHS1 uext=" << uextR[1] << " uint=" << uintR[1] << endl;
    //}
    if(params.metric.schwarschild){
      //specializing to schwazschild case
      if(params.opts.use_generic_orbit){
	if(bleft&&add){//position 
	  uextL.at(0) = uextL.at(0)/coords.dxdxibR0.at(elemnum+1);
	  uextL.at(1)= uextL.at(1)-coords.dxdxibR1.at(elemnum+1)*uextL.at(0);
	}else if (bright&&sub){//position 4

	  uextR.at(0)=uextR.at(0)/coords.dxdxibL0.at(elemnum-1);
	  uextR.at(1)=uextR.at(1)-coords.dxdxibL1.at(elemnum-1)*uextR.at(0);
	}
      }
    }
    //inside the world tube, first handle the singular field.

    if(params.metric.schwarschild){
      if(params.opts.use_world_tube){
	if(found){
	  if(params.opts.useSource){
	    //cout << "DiffEq:use_world_tube:ss " << sstre << " " << sstim << " " << ssrre << " " << ssrim << endl;
	    complex<double> sst(sstre, sstim);
	    complex<double> ssr(ssrre,ssrim);
	    if(bright){
	      uextL.at(1)=uextL.at(1)+sst;
	      uextL.at(0)=uextL.at(0)+ssr;
	    }
	    if(bleft){
	      uextR.at(1)=uextR.at(1)+sst;
	      uextR.at(0)=uextR.at(0)+ssr;
	    }
	  }
	}
      }
    }



    if(params.metric.schwarschild){
      //specializing to schwazschild case
      if(params.opts.use_generic_orbit){
	if(bleft&&sub){
	  uextL.at(1)=uextL.at(1)+uextL.at(0)*coords.dxdxibR1.at(elemnum);
	  uextL.at(0)=uextL.at(0)*coords.dxdxibR0.at(elemnum);
	}else if(bright&&add){
	  uextR.at(1)=uextR.at(1)+uextR.at(0)*coords.dxdxibL1.at(elemnum);
	  uextR.at(0)=uextR.at(0)*coords.dxdxibR0.at(elemnum);
	}
      }
    }



    //below were Array2Ds
    //Initialize plus and minus components of lambda matrix to zero at both
    //boundaries
    vector<double> lambdaminusL(params.grid.Ddim*params.grid.Ddim, 0.0);
    vector<double> lambdaminusR(params.grid.Ddim*params.grid.Ddim, 0.0);
    vector<double> lambdaplusL(params.grid.Ddim*params.grid.Ddim, 0.0);
    vector<double> lambdaplusR(params.grid.Ddim*params.grid.Ddim, 0.0);

    vector<double> lambdaL= AleftBoundaries[elemnum].getLambda();
    vector<double> lambdaR= ArightBoundaries[elemnum].getLambda();


    //lambda minus contains outward moving wave components
    //lambda plus contains inward moving wave components
    //Might be an incorrect summary. Trust the math, not the words
    //See pg 35 of Hesthaven and Warburten


    for(int j = 0; j < params.grid.Ddim; j++) {
      if(nL * lambdaL[j*params.grid.Ddim+j] <= 0) {
	lambdaminusL[j*params.grid.Ddim+j] = nL * lambdaL[j*params.grid.Ddim+j];
      } else {
	lambdaplusL[j*params.grid.Ddim+j] = nL * lambdaL[j*params.grid.Ddim+j];
      }

      if(nR * lambdaR[j*params.grid.Ddim+j] <= 0) {
	lambdaminusR[j*params.grid.Ddim+j] = nR * lambdaR[j*params.grid.Ddim+j];
      } else {
	lambdaplusR[j*params.grid.Ddim+j] = nR * lambdaR[j*params.grid.Ddim+j];
      }
    }


    //S and S inverse matrices at both boundaries
    vector<double> sinvL = AleftBoundaries[elemnum].getSinv();
    vector<double> sinvR = ArightBoundaries[elemnum].getSinv();
    vector<double> SL = AleftBoundaries[elemnum].getS();
    vector<double> SR = ArightBoundaries[elemnum].getS();


    //Numerical fluxes at both boundaries
    //See Hesthaven and Warburten pg 35 (n*F)
    //were Array1D's
    vector<complex<double>> temp = matmul(sinvL, uintL,params.grid.Ddim,params.grid.Ddim,1);
    vector<complex<double>> nfluxL1 = matmul(lambdaplusL, temp, params.grid.Ddim,params.grid.Ddim,1);
    temp = matmul(sinvL, uextL,params.grid.Ddim,params.grid.Ddim,1);
    vector<complex<double>> nfluxL2 =  matmul(lambdaminusL,temp,params.grid.Ddim,params.grid.Ddim,1);

    //    cout << modenum << " " <<elemnum << " " << nfluxL1.at(0) << " " << nfluxL1.at(1) << " " << nfluxL2.at(0) << " " << nfluxL2.at(1) << endl;
    
    temp = vecsum(nfluxL1,nfluxL2);
    vector<complex<double>> nfluxL = matmul(SL,temp,params.grid.Ddim,params.grid.Ddim,1);
    temp = matmul(sinvR, uintR,params.grid.Ddim,params.grid.Ddim,1);
    vector<complex<double>> nfluxR1 = matmul(lambdaplusR,temp,params.grid.Ddim,params.grid.Ddim,1);


    temp = matmul(sinvR, uextR,params.grid.Ddim,params.grid.Ddim,1);
    vector<complex<double>> nfluxR2 = matmul(lambdaminusR,temp,params.grid.Ddim,params.grid.Ddim,1);
    temp = vecsum(nfluxR1,nfluxR2);
    vector<complex<double>> nfluxR = matmul(SR,temp,params.grid.Ddim,params.grid.Ddim,1);


    if(found&&(modenum==1)){
      //cout << modenum << " " << elemnum << " " <<  nfluxL.at(0) << " " << nfluxL.at(1) << " " << nfluxR.at(0) << " " << nfluxR.at(1) << endl;
    }

    //were Array2D's
    vector<double> AtrimmedL= AleftBoundaries[elemnum].getAtrimmed();
    vector<double> AtrimmedR= ArightBoundaries[elemnum].getAtrimmed();
    vector<complex<double>> duelem(params.grid.Adim*2, {0.0,0.0});

    //    if(output){
    //  cout << elemnum << " " << SR[0] << " " << SR[1] << " " << SR[2] << " "<< SR[3] << endl;
    // }
    //were Array1Ds
    //This gets multiplied by lift matrix to calculate flux
    vector<complex<double>> temp1 = matmul(AtrimmedL,uintL,params.grid.Ddim,params.grid.Ddim,1);
    vector<complex<double>> temp2 = matmul(AtrimmedR, uintR, params.grid.Ddim,params.grid.Ddim,1);
    vector<complex<double>> temp3= scalarmult(nL,temp1);
    vector<complex<double>> temp4 = scalarmult(nR,temp2);
    vector<complex<double>> duL = vecdiff(temp3,nfluxL);
    vector<complex<double>> duR = vecdiff(temp4,nfluxR);

    //DU SHOULD BE ZERO AFTER FIRST CALL TO RHS BECAUSE SAME ON EITHER SIDE OF BOUNDARY

    //  if((elemnum==3)&&(output)){
    // for(int i=0; i<AtrimmedR.size(); i++){
    //cout << AtrimmedL[i] << endl;
    // }
    // }

    insert_1D_into_2D_vec(du,duL,2,params.grid.Ddim,0,false);
    insert_1D_into_2D_vec(du,duR,2,params.grid.Ddim,1,false);

    //if(output){
    //cout << modenum << " " << t << " " << thegrid.gridNodeLocations().get(elemnum,nodenumFound) << " " << duL[0].real() << " " << duL[1].real() << " " << duR[0].real() << " "  << duR[1].real() << endl;
    //}


    vector<complex<double>> uint2D(2*params.grid.Ddim);
    insert_1D_into_2D_vec(uint2D,uintL,2,params.grid.Ddim,0,false);
    insert_1D_into_2D_vec(uint2D,uintR,2,params.grid.Ddim,1,false);

    //    if((modenum==0)&&(output)){
    // for(int i=0; i<du.size(); i++){
    //cout << elemnum << "\t" << du[i] << "\t" << uint2D[i] << endl;
    // }
    // }

    //du correct-- roundoff

    /*cout << elemnum << "\t";
    for(int i = 0; i<du.size(); i++){
    cout << du[i].real() << "\t";
    }
    cout << endl;
    */
    //create a 2D array with one column corresponding to the left
    //boundary and one column corresponding to the right boundary
    //of du for the the various components of u.

    //END characteristicFlux
    //BEGIN RHS

    //pragma omp parallel if(NumElem>lmmodes.ntotal) for

    //#pragma omp parallel for shared(uh,modenum,thegrid,SOURCE_VECNUM,RHStdvgf,params) if(uh.TDVGFdim()<=thegrid.numberElements())
    //Maximum index for both A and B matrix
    //    int vmaxAB = ArightBoundaries[elemnum].getparams.grid.Adim() - 1;
    //Minimum index for use with trimmed A matrix.
    //Minimum index for B matrix is zero
    //    int vminA = vmaxAB - ArightBoundaries[elemnum].getDdim() + 1;

    //The B matrix component of the RHS.

    //was Array2D
    vector<complex<double>> RHSB(uh.GFarrDim()
				 *params.grid.Adim);


    for(int nodenum = 0; nodenum < uh.GFarrDim(); nodenum++){
      vector<complex<double>> RHSBpernode; //was Array1D


      //FIX THIS. THIS CAN BE SPED UP BY NOT COPYING IN GETVECTORASARRAY1D
      //Multiply the B matrix times a "vector" of u for each node
      vector<complex<double>> temp = uh.getVectorRange(modenum,elemnum, nodenum,0, vmaxAB, 0);

      vector<double> tempd = Bmatrices.getVector(modenum,elemnum,nodenum);


      RHSBpernode = matmul(tempd,temp,params.grid.Adim,params.grid.Adim,1);

      //if((output)&&(modenum==0)){
      //cout << setprecision(12);
      //cout << RHSBpernode[0] << " " << RHSBpernode[1] << " " << RHSBpernode[2] << endl;
      //}



      //Insert that result into the rows of a larger matrix

      //if((elemnum==3)&&(output)){
      //for(int i =0; i<RHSBpernode.size(); i++){
      //  cout <<i << "\t" <<  RHSBpernode.at(i) << endl;
      //}
      //}

      insert_1D_into_2D_vec(RHSB, RHSBpernode,uh.GFarrDim(),params.grid.Adim,nodenum, false);


    }//end node loop

    //    if (output){
    //for(int i=0; i<RHSB.size(); i++){
    //cout << setprecision(15);
    //cout << RHSB[i].real() << endl;
    // }
    //}
    //A contribution:
    vector<complex<double>> RHSA1(uh.GFarrDim()*params.grid.Ddim);//was Array2D

    //The A contribution needs to be multiplied one node at a time by the
    //trimmed A matrix in a similar manner to the B contribution. But first,
    //we take the spatial derivative across all thegrid.gridNodeLocations().

    int vecNodeDim1=params.grid.Ddim;
    int vecNodeDim2=uh.GFarrDim();
    temp = uh.getVectorNode2D(modenum, elemnum, vminA, vmaxAB, 0, vecNodeDim1,vecNodeDim2);
    temp2 = matmul(Dmatrix,temp,uh.GFarrDim(),uh.GFarrDim(),params.grid.Ddim);
    vector<complex<double>> RHSA1preA = scalarmult(thegrid.jacobian(elemnum), temp2);
    //was Array2D



    //if((modenum==0)&&(output)){
    //for(int i=0; i<RHSA1preA.size(); i+=2){
    //cout << setprecision(15);
    //cout << thegrid.gridNodeLocations().get(elemnum,i/2) << " " << RHSA1preA.at(i).real() << " " << RHSA1preA.at(i+1).real() << endl;
    //}
    //}

    //Multiply each row of RHSA1preA by a different a Atrimmed matrix
    for(int nodenum=0; nodenum < uh.GFarrDim(); nodenum++)
      {
	int M = params.grid.Ddim;
	int N = params.grid.Ddim;
	int K = uh.GFarrDim(); //RHSApreAdim1
	int L = params.grid.Ddim; //RHSApreAdim2

	//M through L are the same

	//cout << M << " " << N << " " << K << " " << L << endl;

	vector<double> tA = trimmedAmatrices.getVector(elemnum,nodenum);

	//cout << tA.size() << endl;

	//tA is same size

	//if((output)&&(modenum==0)){
	//cout << setprecision(15);

	//cout << thegrid.gridNodeLocations().get(elemnum,nodenum) << " " << tA[0] << " " << tA[1] << " " << tA[2] << " " <<tA[3] << endl;
	//}

	for (int i=0; i<M; i++){
	  complex<double> sum = 0;
	  for (int k=0; k<N; k++)
	        sum += -tA[i*N+k]
		  * RHSA1preA.at(nodenum*L+k);
	  RHSA1.at(nodenum*L+i) = sum;
	}
	//Copied and pasted from TmatmultT in TNT2 a
	//with modification of variable first matrix
	//TmatmultT was copied and pasted from matmult in tnt itself
      }//end nodenum loop
    //Negative sign is because of definition of tA
    //A is definied to appear on the left hand side of the differential
    //equation, but this routine calculates the right hand side

    //Needs a multiplication by an A matrix before D
    //but A is position dependent.

    //get rid of moves?? PARALLEL PROBLEM?
    //This is the contribution due to du, or the numerical flux
    //    Array2D<complex<double>> RHSA2 = move(thegrid.jacobian(elemnum)
    //  *move(matmult(thegrid.refelem.getLift(), du[elemnum])));

    //RHSA1.txt
    //if(output){
    //for(int i=0; i<RHSA1.size(); i++){
    //  cout << setprecision(15);
    //cout << elemnum*RHSA1.size()+i << "\t" << RHSA1.at(i).real() << endl;
    //}
    //}

    temp =matmul(Lift,du,uh.GFarrDim(),2,2);
    vector<complex<double>> RHSA2 = scalarmult(thegrid.jacobian(elemnum),temp);

    //if(output){
    //for(int i=0; i< RHSA2.size(); i++){
    //cout << setprecision(15);
    // cout << elemnum*RHSA2.size()+i << "\t" << RHSA2.at(i).real() << endl;
    //}
    //}


    //RHSA2 matches. Problem was transposed conversion from vectors to TNT in CharacteristicFlux.cpp. trimmedA may still not match to 1e-7 not sure

    //  if(output){
    // cout << elemnum << "\t" << du[0].real() << "\t" << du[1].real() << "\t" << du[2].real()<< "\t" << du[3].real() << endl;
    //}

    //if(output){
    //for(int i=0; i<temp.size(); i++){
    //cout << setprecision(15);
    //cout << elemnum*temp.size()+i << "\t" << temp.at(i).real() << endl;
    //}
    //}

    //LIFT IS FINE
    //if(output){
    // for(int i=0; i<Lift.size(); i++){
    //cout << setprecision(15);
    //cout << elemnum*Lift.size()+i << "\t" << Lift.at(i) << endl;

    //}
    //}


    //if((elemnum==3)&&(output)){
    // for(int i=0; i<RHSA2.size(); i++){
    //cout << i << "\t" << RHSA2.at(i) << endl;
    // }
    //}

    //RHSA and RHSB will have different sizes due to the different
    //number of diffeq variables stored in each. sum them using a
    //for loop while assigning values to the RHStdvgf vector grid function

    //    Array2D<complex<double>> RHSA = RHSA1 + RHSA2;

    //Sum the contributions from B, derivative, and flux,
    //accounting for different matrix dimensions
    for(int vecnum = 0; vecnum < RHStdvgf.VGFdim(); vecnum++){
      for(int nodenum = 0; nodenum < RHStdvgf.GFarrDim(); nodenum++){
	complex<double> tot=-RHSB.at(nodenum*params.grid.Adim+vecnum);
	if(vecnum>=vminA){
	  /*if((output)&&(modenum==1)&&(vecnum==2)){
	    cout << setprecision(15);
	    cout << nodenum << " " << RHSA[nodenum][vecnum-vminA].real() << endl;
	    }*/
	  tot+=RHSA1.at(nodenum*params.grid.Ddim+(vecnum-vminA))
	    + RHSA2.at(nodenum*params.grid.Ddim+(vecnum-vminA));
	}//-sign in B because it is on the left hand side of the
	if(vecnum==0){
	  //  cout << modenum << " " << nodenum << " " << tot << endl;
	}
	// I AM HERE
	if((params.opts.useSource)&&(vecnum==SOURCE_VECNUM)){
	  complex<double> temp = source.get(modenum,elemnum,nodenum);
	  //if((modenum==1)&&(elemnum==15)&&(nodenum==16)){
	  //  cout << t << " " << temp.real() << " " <<  temp.imag() << endl;
	  //}
	  tot += temp;

	  //equation in the definition supplied in DiffEq.cpp
	}
	RHStdvgf.set(modenum, vecnum, elemnum, nodenum, tot);
	//cout << setprecision(15);
	//cout << modenum << " " <<  vecnum << " " << elemnum << " " <<  nodenum << " " <<  tot.real() << " " << tot.imag() << endl;
      }
    }


  }

}

void DiffEq::modeRHS(Grid& thegrid,
                     TwoDVectorGridFunction<complex<double>>& uh,
                     TwoDVectorGridFunction<complex<double>>& RHStdvgf,
                     double t, bool output, Orbit* orb, WorldTube* wt, Coordinates& coords, double& max_speed, Modes& lmmodes)
{
  double maxspeed=1.0; 
  if(orb->orbType()==elliptical){
    double rp, drpdt, d2rpdt2;
    EllipticalOrbit * eorb = dynamic_cast<EllipticalOrbit *>(orb);
    eorb->orb_of_t(coords, rp, drpdt, d2rpdt2);
    coords.timedep_to_rstar(rp, drpdt, d2rpdt2);
    vector<double> dxdt(params.grid.elemorder+1),dxdxi(params.grid.elemorder+1);
      
    for(int elemnum=1; elemnum<params.grid.numelems; elemnum++){
      if(coords.timeDepTrans.at(elemnum-1)){
	set_coefficients(thegrid, eorb, coords, maxspeed, elemnum, lmmodes, rp,rstar_orb, drpdt, d2rpdt2, dxdt, dxdxi);
	//cout << elemnum << " " << dxdxi.at(0) << " " << dxdt.at(0) << " " << dxdxi.at(params.grid.elemorder) << " " << dxdt.at(params.grid.elemorder) << endl;
	coords.dxdxibL0.at(elemnum)=dxdxi.at(0);
	coords.dxdxibL1.at(elemnum)=dxdt.at(0);
	coords.dxdxibR0.at(elemnum)=dxdxi.at(params.grid.elemorder);
	coords.dxdxibR1.at(elemnum)=dxdt.at(params.grid.elemorder);
	//	cout << elemnum << " " << coords.dxdxibL0.at(elemnum) << " " << coords.dxdxibL1.at(elemnum) << " "<< coords.dxdxibR0.at(elemnum) <<" " << coords.dxdxibR1.at(elemnum) << endl;
	
      }
      
      if(!params.opts.use_world_tube){
	vector<double> rschwv = thegrid.rschw.get(elemnum);
	double * rarr = &rschwv[0];
      
	
	vector<double> windowv = thegrid.window.get(elemnum);
	vector<double> dwindowv = thegrid.dwindow.get(elemnum);
	vector<double> d2windowv = thegrid.d2window.get(elemnum);
	
	
	double * winarr = &windowv[0];
	double * dwinarr = &dwindowv[0];
	double * d2winarr = &d2windowv[0];
	
	  
	calc_window(params.grid.elemorder+1, rarr, winarr, dwinarr, d2winarr);
	thegrid.window.set(elemnum,windowv);
	thegrid.dwindow.set(elemnum,dwindowv);
	thegrid.d2window.set(elemnum,d2windowv);
	thegrid.rschw.set(elemnum,rschwv);
	//put inside if time dependent loop?
	setupABmatrices(thegrid,lmmodes,coords);
      }
      max_speed=max(maxspeed,max_speed);
      
    }
  }
  

  if(params.opts.useSource) {
    
    fill_source_all(thegrid, t, uh.TDVGFdim(), source, thegrid.window,
		    thegrid.dwindow, thegrid.d2window, orb, lmmodes);
  }
  

  typedef numeric_limits<double> dbl;
  /*for(int k = 0; k<lmmodes.ntotal; k++){
    ofstream fss;
    fss.precision(dbl::max_digits10);
    fss.precision(16);
    ostringstream oss;
    oss << "source" << "." << k << ".txt";
    fss.open(oss.str(), ios::app);
    fss << endl << endl;
    fss << " #time = " << t << endl;
    
  
  
    for(int i=0; i<source.GFvecDim(); i++){
      for(int j=0; j<source.GFarrDim(); j++){
	fss << thegrid.gridNodeLocations().get(i,j) << " "<< source.get(k,i,j).real() << " " << source.get(k,i,j).imag() << endl;
      }
    }
    
    fss.close();
    }*/
  /*if (params.opts.useSource){
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
    }*/
  
  
  //#pragma omp parallel for if(uh.TDVGFdim()>thegrid.numberElements())


  
  for(int modenum = 0; modenum < uh.TDVGFdim(); modenum++) {
       //  double max_speed = 1.0;

    RHS(modenum, thegrid, uh, RHStdvgf, t, output, coords,wt);
  }
}
