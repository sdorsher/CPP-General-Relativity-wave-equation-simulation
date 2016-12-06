#ifndef ELLIPTICAL_ORBIT_H
#define ELLIPTICAL_ORBIT_H
#include "namespaces.h"
#include "ConfigParams.h"
#include "Orbit.h"
#include "Coordinates.h"

using namespace std;
//using namespace orbit;

class EllipticalOrbit:public Orbit{
 public:
  double dtdchi, dphidchi;
  double dchidt, d2chidt2, drpdt, d2rpdt2, dphidt, reschi, resphi;
  EllipticalOrbit();
  void dorbdchi();
  void dorbdt();
  void orb_of_t(Coordinates & coords);
};
#endif

  
