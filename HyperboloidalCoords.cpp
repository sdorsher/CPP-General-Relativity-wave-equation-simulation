#include "HyperboloidalCoords.h"


/* Function to calculate Lambert's W-function
Algorithm found at
http://en.citizendium.org/wiki/Lambert_W_function#Numerical_calculation
*/

double Lambert(double z)
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

double rschw(double z, double mass)
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
double invert_tortoise(double rstar, double mass)
{
  double invtort;
  if (rstar >= 0.0) {
    invtort = rschw(rstar, mass);
  } else {
    invtort = 2.0 * mass * Lambert( exp(0.5 * rstar / mass - 1.0));
  }
  return invtort;
}

void transition(double rho, double R, double S, double& fT, double& fTp,
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
  fT = half + half * tanhfac;
  cscx2 = 1.0 / pow(sin(x),2.0);
  secx2 = 1.0 / pow(cos(x),2.0);
  sechfac2 = 1.0 / pow(cosh(fac),2.0);
  fTp = 0.25 * es * (q2 * cscx2 + secx2) * sechfac2 / (S - R);
  fTpp = -0.25 * es * sechfac2 * 
    (PI * (q2 * cotx * cscx2 - secx2 * tanx) 
     + es * pow((q2 * cscx2 + secx2), 2.0) * tanhfac) / pow((S - R),2.);
}
