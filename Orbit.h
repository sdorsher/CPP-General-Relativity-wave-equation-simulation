

#ifndef ORBIT_H
#define ORBIT_H

#include "globals.h"
#include "ConfigParams.h"


enum OrbitType{circular, elliptical};

class Orbit{
 protected:
  OrbitType otype;
 public:
  double p; //simulatus rectum
  double e; //eccentricity -- zero for circular orbits
  double chi; //angular parameter for radial oscillations
  double phi; //angular parameter for angular oscillations
  double drdlambda_particle, drdxi_particle;
  double E; //particle energy, neglecting background spacetime reactions or
  // scalar field energy density
  double L; //angular momentum, same deal. See Wald chapter on Schwarzschild geodesics or Pound and Poisson. 
  virtual ~Orbit(){};
  Orbit();
  OrbitType orbType();
  
  
};


#endif
