#ifndef ELLIPTICAL_ORBIT_H
#define ELLIPTICAL_ORBIT_H
#include "namespaces.h"
#include "ConfigParams.h"
#include "Orbit.h"

using namespace std;
//using namespace orbit;

class EllipticalOrbit:public Orbit{
 public:
  double dtdchi, dphidchi, dchidt, dphidt;
  double dchidt, d2chidt2, drpdt, d2rpdt2, dphidt, reschi, resphi;
  void EllipticalOrbit();
  void dorbdchi();
  void dorbdt();
  void orb_of_t();
};
#endif

  
