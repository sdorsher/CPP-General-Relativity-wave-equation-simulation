#ifndef CIRCULARORBIT_H
#define CIRCULARORBIT_H

#include "namespaces.h"
#include <cmath>
#include "ConfigParams.h"
#include <iomanip>
#include "Orbit.h"

//using namespace orbit;

class CircularOrbit:public Orbit{
 public:
  //set parameters in orbit namespace
  CircularOrbit();
  
  //calculate phi as a function of time
  double phi_of_t(double t);
  
  //calculate chi as a function of time
  double chi_of_t(double t);


  double circ_E();
  double circ_L();

};
 
#endif
