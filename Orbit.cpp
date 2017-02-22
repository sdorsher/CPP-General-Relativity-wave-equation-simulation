#include "Orbit.h"


Orbit::Orbit()
{  //initialize orbit

  phi=0.;
  chi=PI;
  p=params.schw.p_orb/params.schw.mass;
  e=params.schw.ecc;
  //chi apastron
}

OrbitType Orbit::orbType(){
  return otype;
}

