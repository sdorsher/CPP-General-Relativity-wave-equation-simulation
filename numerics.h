#ifndef NUMERICS_H
#define NUMERICS_H

#include <cmath>

using namespace std;

void time_window(const double time,const  double tsigma,const  int norder, double& tfac, double & dtfac_dt, double & d2tfac_dt2);  

#endif
