#include "numerics.h"

void time_window(const double  time, const double tsigma, const int norder, double& tfac, double & dtfac_dt, double & d2tfac_dt2) {
  //tested and confirmed
  double tfactor, expfactor;
  tfactor = pow((time/tsigma), double(norder));
  expfactor = exp(-tfactor);
  tfac = 1.0 - expfactor;
  dtfac_dt = norder * pow(time, (norder-1.0)) / pow(tsigma, double(norder)) 
    * expfactor;
  d2tfac_dt2 = -norder * (1.0 + norder * (tfactor - 1.0))
    * pow(time, (norder-2.0))/pow(tsigma, double(norder)) * expfactor;
}
