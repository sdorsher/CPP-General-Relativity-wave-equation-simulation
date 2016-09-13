#include "Coordinates.h"


Coordinates::Coordinates(),dxdxib{params.grid.numelems,0.0},dxidbL{params.grid.elemorder+1,0.},dxidbR{params.grid.elemorder+1,0.}{}

double Coordinates::rstar_of_r(double r, double mass)
{
  return r + 2.0*mass * log(r/(2.0*mass)-1.0);
}

/* Function to calculate Lambert's W-function
Algorithm found at
http://en.citizendium.org/wiki/Lambert_W_function#Numerical_calculation
*/

double Coordinates::Lambert(double z)
{
  double wcurrent, wnew, expw, diff;
  const double eps = 1.0e-12;
  wcurrent =1.0;
  bool finished=false;

  while(!finished){
    expw = exp(wcurrent);
    diff = wcurrent * expw - z;
    wnew = wcurrent - diff / (expw * (wcurrent + 1.) 
                              - ((wcurrent + 2.) * diff 
                                 / (2. * wcurrent + 2.)));
    if(fabs(wnew-wcurrent) < eps) finished=true;
    wcurrent = wnew;
  }
  return wnew;
}

/*Function to invert the tortoise radius as a function of Schwarzchild radius.
Returns rschw-2M.
Simple Newton root finding algorithm. Fortran version has problems 
converging for negative 
rstar as it overshoots and rschw-2M can become negative. 
*/

double Coordinates::rschw(double z, double mass)
{
  double xcurrent, xnew;
  const double eps=1.0e-12;
  
  xnew=0.0;
  xcurrent = 1.0;

  bool finished =false;

  while(!finished){
    xnew = (z - 2.0 * mass * log(0.5 * xcurrent / mass))
      /(1. + 2. * mass / xcurrent);
    if(fabs(xnew - xcurrent) < eps) 
      {
        finished=true;
      }
    xcurrent = xnew;
  }
  return xnew;
}



/* Function to invert the tortoise radius as a function of Schwarzchild radius.
Returns rschw-2M. 
Uses rschw for rstar>=0 and Lambert for rstar<0. The Lambert method has
problems for large rstar due to numerical overflow of it's arguments. 
*/
//invert_tortoise
double Coordinates::invert_tortoise(double rstar, double mass)
{
  double invtort;
  if (rstar >= 0.0) {
    invtort = rschw(rstar, mass);
  } else {
    invtort = 2.0 * mass * Lambert( exp(0.5 * rstar / mass - 1.0));
  }
  return invtort;
}

void Coordinates::transition(double rho, double R, double S, double& fT, double& fTp,
                  double& fTpp) {
  const double es = 1.5;
  const double q2 = 1.3;
  const double half = 0.5;
  double x, tanx, cotx, fac, tanhfac, cscx2, secx2, sechfac2;

  x = half * PI * (rho - R) / (S - R);
  tanx = tan(x);
  cotx = 1.0 / tanx;
  fac = es / PI * (tanx - q2 * cotx);
  tanhfac = tanh(fac);
  fT = half + half * tanhfac; //transition function
  cscx2 = 1.0 / pow(sin(x),2.0);
  secx2 = 1.0 / pow(cos(x),2.0);
  sechfac2 = 1.0 / pow(cosh(fac),2.0);

  //derivative of transition function wrt coordinate rho
  fTp = 0.25 * es * (q2 * cscx2 + secx2) * sechfac2 / (S - R);

  //second derivative of transition function wrt rho
  fTpp = -0.25 * es * sechfac2 * 
    (PI * (q2 * cotx * cscx2 - secx2 * tanx) 
     + es * pow((q2 * cscx2 + secx2), 2.0) * tanhfac) / pow((S - R),2.);
}


void Coordinates::timedep_to_rstar(Orbit& orb){
  orb.d2rpdt2 = ( -2.0 * params.schw.mass * pow(orb.drpdt,2.) 
	      + orb.rp * ( orb.rp - 2.0 * params.schw.mass ) * orb.d2rpdt2 )
    /pow(( orb.rp - 2.0 * params.schw.mass ),2.);
  orb.drpdt = orb.rp / ( orb.rp - 2.0 * params.schw.mass ) * orb.drpdt;
  orb.rp = rstar_of_r ( orb.rp, params.schw.mass );
}


void Coordinates::coord_trans(Coordinates &coords, Grid& thegrid, vector<double> & dxdxi, vector<double> & d2dxdt2, vector<double> & d2dxdxi2,vector<double> & d2xdtdxi, int elemnum){
  //time dep coord transf
  double xpma, xipma, xima, bmxp, bmxip, bmxi, bma, ximxip, xipmxp, xipmainv, xipamulbmxipinv, dtfac;
  for(int i=0; i<xi.size(); i++){
    
    //UNFINISHED HERE. Size is one params.grid.elemorder+1. pass in one element at a time. look at reference for arrays peter gave me. 

    j= elemnum;
    double a = coords.a;
    double b = coords.b;
    double xp = coords.xp;
    xi=thegrid.rho.get(j);
    x = thegrid.rstar.get(j);
    
    xpma=xp-a;
    xipma=xp-a;
    xima=xi.at(i)-a;
    bmxp=b-xp;
    bmxi=b-xi.at(i);
    bmxip=b-xip;
    bma=b-a;
    ximxip=xi.at(i)-xip;
    xipmainv=1./xipma;
    xipmamulbmxipinv=xipmainv/bmxip;
    xipmxp=xip-xp;
    xi.at(i)=a+xpma*xipmainv*xima
      +(bmxp*xipma-xpma*bmxip)*xipmamulbmxipinv/bma*xima*ximxip;
    dtfac=xima*bmxi*xipmamulbmxipinv;
    dxdt.at(i)=dtfac*dxpdt;
    dxdxi.at(i)=((2.*xi.at(i)-xip-a)*xipmxp+xpma*bmxip)
      *xipmamulbmxipinv;
    d2xdt2.at(i)=dtfac*d2xpdt2;
    d2xdxi2.at(i)=2.*xipmxp*xipmamulbmxipinv;
    d2dtdxi(i)=(a+b+2.*xi.at(i))*xipmamulbmxipinv*dxpdt;
  }
}
