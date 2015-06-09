#ifndef CONFIGPARAMS_H
#define CONFIGPARAMS_H

#include <fstream>
#include <libconfig.h++>
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace libconfig;

struct WaveEqParams{
  int pdenum; //number of components in coupled differential equations
  double speed; //speed of wave
  bool isgaussian; //true if wave is gaussian and zero time derivative
  bool issinusoid; //true if wave is sinusoidal and travelling
};

struct SinusoidParams{
  double phase; 
  double amp;
  double wavelength; //most sensible if this is a standing wave
};

struct GaussParams{
  double amp;
  double mu;
  double sigma; //sigma greater than or equal to element size works best
};

struct GridParams{
  double lowerlim; //lower boundary of grid
  double upperlim; //upper boundary of grid
  int numelems; //number of elements in grid
  int elemorder; //order of each element (all are the same)
  bool readfromfile; //broken right now, true or false
};

struct TimeParams{
  double dt; //time step if time step is fixed
  double courantfac; //ratio between smallest spatial step and time step 
                     //if time step is not fixed
  double t0; //initial time
  double tmax; //end time
  double outputinterval; //output with this time step
  int comparisoncount; //use the wave data to evaluate the L2 norm after 
                       //this integer number of output intervals
  bool usefixedtimestep; 
};

struct FileParams
{//file names
  string pdesolution; //stores the pde solution at each output time
  string oneperioderror; //the error as a function of position after 
                         //one oscillation
  string L2error; //order, timestep, number of elements, L2error
  string initialconditions; //initial conditions as a function of position
};

//main structure to be used as global variable
struct ConfigParams {
  WaveEqParams waveeq;
  SinusoidParams sine;
  GaussParams gauss;
  GridParams grid;
  TimeParams time;
  FileParams file;

  ConfigParams(const std::string& configFileName);

  template <typename T>
  T getConfigFromFile( const std::string& configFileName, 
                       const std::string& structkey, const std::string& key);
};

extern const ConfigParams params; //global params variable

#endif
