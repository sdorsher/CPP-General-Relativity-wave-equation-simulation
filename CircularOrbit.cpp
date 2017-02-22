#include "CircularOrbit.h"

//Set initial orbit parameters
CircularOrbit::CircularOrbit() {
  otype = circular;
  p = params.schw.p_orb; //semi-latus rectum
  e = params.schw.ecc; //eccentricity
  chi = acos(-1.0); //2pi in one full radial oscillation
  phi = 0.0; //2pi in one full angular oscillation
  E=circ_E();
  L=circ_L();
}

// phi as a function of time
double CircularOrbit::phi_of_t(double t)
{
  //  cout << setprecision(15);
  //  cout << params.schw.mass << " " << params.schw.p_orb << endl;

  //  cout << setprecision(15);
  //  cout << t << endl;
  
  double omega = sqrt(params.schw.mass/pow(p*params.schw.mass,3.0));
  //  cout << setprecision(15);
  //  cout << "before: " << t << " " << params.schw.mass << " " << params.schw.p_orb << " " << omega*t << endl;
  
  
  phi=omega*t;
  if(phi>2.0*PI) phi-=2.0*PI;
  return phi;
  
  
}

// chi as a function of time
double CircularOrbit::chi_of_t(double t)
{
  double chiomega = params.schw.mass * sqrt(p-6.0)/pow(p*params.schw.mass, 2.0);
  chi=chiomega*t;
  if (chi>2.0*PI) chi-=2.0*PI;
  return chi;
}



double CircularOrbit::circ_E(){
  double r = p * params.schw.mass;
  double m = params.schw.mass;
  E=(r-2*m)/sqrt(r)/sqrt(r-3*m);
  return E;

}

double CircularOrbit::circ_L(){
  double r=p * params.schw.mass;
  double m = params.schw.mass;
  if (p>6){
    L = r*sqrt(m)/sqrt(r-3*m);
  }else {
    L=0;
    //fix me!
  }
  return L;
}
    
