#ifndef ELLIPTICAL_ORBIT_H
#define ELLIPTICAL_ORBIT_H
#include "namespaces.h"


using namespace std;
using namespace orbit;

class EllipticalOrbit{
  void EllipticalOrbit();
  void dorbdchi(double &dtdchi, double &dphidchi);
  void dorbdt();
  void orb_of_t(double &rp, double &drpdt, double &d2rpdt2);
}
#endif

  
