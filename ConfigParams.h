#ifndef CONFIGPARAMS_H
#define CONFIGPARAMS_H

#include <fstream>
#include <libconfig.h++>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>

using namespace std;
using namespace libconfig;

struct OptionsParams{
  bool useSource; //use the effective source
  bool turn_on_source_smoothly; //use a window function in time to turn on
  //the effective source smoothly
};

struct MetricParams{
  bool flatspacetime; //use flat space-time
  bool schwarschild; //use a Schwarzschild background
};

struct WaveEqParams{ //set only one of issinusoid or isgaussian to true
  double speed; //speed of wave
  bool isgaussian; //true if wave is gaussian and zero time derivative
  bool issinusoid; //true if wave is sinusoidal and travelling
};

struct SinusoidParams{ //used if issinusoid is set to true
  double phase; 
  double amp;
  double wavelength; //most sensible if this is a standing wave
};

struct GaussParams{ //used if isgaussian is set to true
  double amp;
  double mu;
  double sigma; //sigma greater than or equal to element size works best
};

struct SchwParams{ //used if schwarzschild is set to true under metric
  double mass;
  double sigma;
  double p_orb; //location of orbit
  double ecc; //eccentricity of orbit
};

struct WindowParams{
  int noffset;
};

struct TimeWindowParams{
  double tsigma;
  int torder;
};

struct ModeParams{
  int lmax;
};

struct GridParams{
  int pdenum; //number of components in coupled differential equations
  double lowerlim; //lower boundary of grid
  double upperlim; //upper boundary of grid
  int numelems; //number of elements in grid
  int elemorder; //order of each element (all are the same)
  bool readfromfile; //broken right now, true or false
  double outputradius;
};

struct HyperbParams{ //hyperboloidal coordinate parameters
  double Splus; //scri plus
  double Sminus; //scri minus
  double Rplus; //boundary of outer hyperboloidal region
  double Rminus; //boundary of inner hyperboloidal region
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
  int outputevery; // output every N iterations

};

struct FileParams
{//file names
  string pdesolution; //stores the pde solution at each output time
  string oneperioderror; //the error as a function of position after 
                         //one oscillation
  string L2error; //order, timestep, number of elements, L2error
  string initialconditions; //initial conditions as a function of position
  bool outputtimefixed; 
  bool outputradiusfixed; 
  string fixedradiusfilename;

};

//main structure to be used as global variable
struct ConfigParams {
  OptionsParams opts;
  MetricParams metric;
  WaveEqParams waveeq;
  SinusoidParams sine;
  GaussParams gauss;
  ModeParams modes;
  GridParams grid;
  HyperbParams hyperb;
  TimeParams time;
  FileParams file;
  SchwParams schw;
  WindowParams window;
  TimeWindowParams timewindow;

  ConfigParams(const std::string& configFileName);

  template <typename T>
  T getConfigFromFile( const std::string& configFileName, 
                       const char *structkey, const char *key);
};

extern const ConfigParams params; //global params variable

#endif
