

#ifndef ORBIT_H
#define ORBIT_H

enum OrbitType{circular, elliptical};

class Orbit{
 protected:
  OrbitType otype;
 public:
  double rp, drpdt, d2rpdt2;
  double p; //simulatus rectum
  double e; //eccentricity -- zero for circular orbits
  double chi; //angular parameter for radial oscillations
  double phi; //angular parameter for angular oscillations
  double drdlambda_particle, drdxi_particle;
  virtual ~Orbit(){};
  Orbit();
  virtual OrbitType orbType();
  
};


#endif
