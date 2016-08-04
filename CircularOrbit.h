#ifndef CIRCULARORBIT_H
#define CIRCULARORBIT_H

#include "namespaces.h"
#include <cmath>
#include "ConfigParams.h"
#include <iomanip>


using namespace orbit;

class CircularOrbit{
  //set parameters in orbit namespace
  void CircularOrbit();

  //calculate phi as a function of time
  double CircularOrbit::phi_of_t(double t);
  
  //calculate chi as a function of time
  double CircularOrbit::chi_of_t(double t);

}
 
#endif
