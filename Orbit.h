

#ifndef ORBIT_H
#define ORBIT_H

class Orbit{
 public:
  double rp, drpdt, d2rpdt2;
  double p; //simulatus rectum
  double e; //eccentricity -- zero for circular orbits
  double chi; //angular parameter for radial oscillations
  double phi; //angular parameter for angular oscillations
  double drdlambda_particle, drdxi_particle;
  virtual ~Orbit(){};
  Orbit();

};


#endif
