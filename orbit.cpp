#include "orbit.h"

void initialize_orbit() {
  p = params.schw.p_orb/params.schw.mass;
  e = params.schw.ecc;
  chi = acos(-1.0);
  phi = 0.0;
}

double phi_of_t(double t)
{
  double omega = sqrt(params.schw.mass/pow(params.schw.p_orb,3.0));
  return omega*t;
}

double chi_of_t(double t)
{
  double chiomega = params.schw.mass * sqrt(params.schw.p_orb/params.schw.mass
                                     -6.0)/pow(params.schw.p_orb, 2.0);
  return chiomega*t;
}

