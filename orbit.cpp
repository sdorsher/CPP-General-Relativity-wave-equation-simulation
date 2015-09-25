#include "orbit.h"

//Set initial orbit parameters
void initialize_orbit() {
  p = params.schw.p_orb/params.schw.mass; //semi-latus rectum
  e = params.schw.ecc; //eccentricity
  chi = acos(-1.0); //2pi in one full radial oscillation
  phi = 0.0; //2pi in one full angular oscillation
}

// phi as a function of time
double phi_of_t(double t)
{
  double omega = sqrt(params.schw.mass/pow(params.schw.p_orb,3.0));
  return omega*t;
}

// chi as a function of time
double chi_of_t(double t)
{
  double chiomega = params.schw.mass * sqrt(params.schw.p_orb/params.schw.mass
                                     -6.0)/pow(params.schw.p_orb, 2.0);
  return chiomega*t;
}

