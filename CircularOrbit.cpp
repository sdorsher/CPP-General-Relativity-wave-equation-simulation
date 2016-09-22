#include "CircularOrbit.h"

//Set initial orbit parameters
CircularOrbit() {
  otype = circular;
  p = params.schw.p_orb/params.schw.mass; //semi-latus rectum
  e = params.schw.ecc; //eccentricity
  chi = acos(-1.0); //2pi in one full radial oscillation
  phi = 0.0; //2pi in one full angular oscillation
}

// phi as a function of time
double CircularOrbit::phi_of_t(double t)
{
  //  cout << setprecision(15);
  //  cout << params.schw.mass << " " << params.schw.p_orb << endl;

  //  cout << setprecision(15);
  //  cout << t << endl;
  

  double omega = sqrt(params.schw.mass/pow(params.schw.p_orb,3.0));

  //  cout << setprecision(15);
  //  cout << "before: " << t << " " << params.schw.mass << " " << params.schw.p_orb << " " << omega*t << endl;

  return omega*t;
  
  
}

// chi as a function of time
double CircularOrbit::chi_of_t(double t)
{
  double chiomega = params.schw.mass * sqrt(params.schw.p_orb/params.schw.mass
                                     -6.0)/pow(params.schw.p_orb, 2.0);
  return chiomega*t;
}

OrbitType CircularOrbit::orbType(){
  return otype;
}
