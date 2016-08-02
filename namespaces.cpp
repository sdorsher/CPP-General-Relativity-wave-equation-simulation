#include "namespaces.h"

namespace orbit{
  double p; //simulatus rectum
  double e; //eccentricity -- zero for circular orbits
  double chi; //angular parameter for radial oscillations
  double phi; //angular parameter for angular oscillations
  double dchidt, dphidt, reschi, resphi;
  double drdlambda_particle, drdxi_particle;
}

namespace layers{
  double Splus; //position of Scri plus in hyperboloidal coordinates
  double Sminus; //position of the horizon in hyperbolidal coordinates
  double Rplus; // position of the transition in either coordinate
  double Rminus; // position of the transition in either coordinate
  double Wplus; //boundaries of the window function in tortoise coordinates
  double Wminus;
}
namespace window
{//Window parameters in Schwarzschild coordinates
  double R1, R2, w1, w2, s1, s2, q1, q2, nmodes;
}
