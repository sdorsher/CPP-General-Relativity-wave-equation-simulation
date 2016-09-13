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
  int L2outputcount;
  bool useSource; //use the effective source
  bool turn_on_source_smoothly; //use a window function in time to turn on
  //the effective source smoothly
  bool use_generic_orbit;
  bool use_world_tube;
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
  double phase; //phase of sinusoid
  double amp; //amplitude of sinusoid
  double wavelength; //most sensible if this is a standing wave
};

struct GaussParams{ //used if isgaussian is set to true
  double amp; //amplitude of gaussian
  double mu; //mean of gaussian
  double sigma; //sigma greater than or equal to element size works best
};

struct SchwParams{ //used if schwarzschild is set to true under metric
  double mass; //mass of black hole
  double sigma; //width of initial scalar field perturbation if effective source is not used
  double p_orb; //location of orbit
  double ecc; //eccentricity of orbit

};

struct WindowParams{
  int noffset; //Offset of end of window around particle from end of tortoise region.
  //The window ends inside the tortoise region. This is given in units of grid points.
};

struct TimeWindowParams{
  double tsigma; //Width of time window.
  int torder; //Order of time window.
};

struct ModeParams{
  int lmax; //Maximum l mode to be computed.
};


struct GridParams{
  int Ddim;
  int Adim;
  int pdenum; //number of components in coupled differential equations
  double lowerlim; //lower boundary of grid
  double upperlim; //upper boundary of grid
  int numelems; //number of elements in grid
  int elemorder; //order of each element (all are the same)
  bool readfromfile; //broken right now, true or false
  double outputradius; //Radius at which to output data
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
  int comparisoncount; //use the wave data to evaluate the L2 norm after 
                       //this integer number of output intervals
  int outputevery; // output every N iterations

};

struct FileParams
{//file names
  string pdesolution; //stores the pde solution at each output time
  string oneperioderror; //the error as a function of position after 
                         //one oscillation
  string L2error; //order, timestep, number of elements, L2error
  string initialconditions; //initial conditions as a function of position
  bool outputtimefixed; //output at a fixed time: true or false
  bool outputradiusfixed; //output at a fixed radius: true or false
  string fixedradiusfilename; //name of file for outputing at a fixed radius

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
  // use these as, for example, params.schw.mass
  
  ConfigParams(const std::string& configFileName);

  //this function reads in parameters of a specified template type
  // from a file using libconfig
  template <typename T>
  T getConfigFromFile( const std::string& configFileName, 
                       const char *structkey, const char *key);
};

//the global variable params, available if "ConfigParams.h" is included
extern const ConfigParams params; //global params variable

#endif
