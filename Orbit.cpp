#include "Orbit.h"


Orbit::Orbit():
 phi{0.},chi{PI}
{  //initialize orbit
  p=params.schw.p_orb/params.schw.mass;
  e=params.schw.ecc;
  //chi apastron
}
