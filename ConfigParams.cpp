#include "ConfigParams.h"
#include <string>

const ConfigParams params("params.cfg");

ConfigParams::ConfigParams(const std::string& configFileName)
{ //Use template function that uses libconfig commands to read in global
  //parameters into a params structure, for example, params.waveq.pdenum

  //Read parameters associated with the metric
  metric.flatspacetime = getConfigFromFile<bool>(configFileName, "metric", 
                                                 "flatspacetime");
  metric.schwarschild =  getConfigFromFile<bool>(configFileName, "metric", 
                                                 "schwarschild");
  //if schwartschild, also hyperboloidal

  if(metric.flatspacetime) {
    //Read parameters associated with the wave equation
    waveeq.pdenum=getConfigFromFile<int>(configFileName, "waveeq", "pdenum");
    waveeq.speed=getConfigFromFile<double>(configFileName, "waveeq", "speed");
    waveeq.isgaussian=getConfigFromFile<bool>(configFileName,"waveeq",
                                              "isgaussian");
    waveeq.issinusoid=getConfigFromFile<bool>(configFileName, "waveeq",
                                              "issinusoid");

    //If it is gaussian, read parameters associated with gaussian 
    //initial conditions
    if(waveeq.isgaussian){
      gauss.amp=getConfigFromFile<double>(configFileName, "gauss", 
                                          "amplitude");
      gauss.mu= getConfigFromFile<double>(configFileName, "gauss", "mu");
      gauss.sigma= getConfigFromFile<double>(configFileName, "gauss", "sigma");
      
    } else if(waveeq.issinusoid) { 
      //If sinusoidal initial conditions, read in parameters 
    //associated with those
      sine.amp=getConfigFromFile<double>(configFileName, "sine", "amplitude");
      sine.phase=getConfigFromFile<double>(configFileName, "sine", "phase");
      sine.wavelength=getConfigFromFile<double>(configFileName, "sine",
                                                "wavelength");
    }else{
      //Currently no other options
      throw invalid_argument("No initial conditions set in parameter file.");
    }
  } else {
    schw.mass = getConfigFromFile<double>(configFileName, "schw", "mass");
    schw.sigma = getConfigFromFile<double>(configFileName, "schw", "sigma");
  }

  //Read in parameters associated with the modes
  modes.lmax = getConfigFromFile<int>(configFileName,"modes", "lmax");

  //Read in parameters associated with the grid
  grid.lowerlim=getConfigFromFile<double>(configFileName, "grid", "lowerlim");
  grid.upperlim=getConfigFromFile<double>(configFileName, "grid", "upperlim");
  grid.numelems=getConfigFromFile<int>(configFileName, "grid", "numelems");
  grid.elemorder=getConfigFromFile<int>(configFileName, "grid", "elemorder");
  grid.readfromfile=getConfigFromFile<bool>(configFileName, "grid",
                                            "readfromfile");
  //Read in parameters associated with hyperboloidal coordinates
  hyperb.Sminus = getConfigFromFile<double>(configFileName, "hyperb", "Sminus");
  hyperb.Splus = getConfigFromFile<double>(configFileName, "hyperb", "Splus");
  hyperb.Rplus = getConfigFromFile<double>(configFileName, "hyperb", "Rplus");
  hyperb.Rminus = getConfigFromFile<double>(configFileName, "hyperb", "Rminus");
  

  
  //Read in parameters associated with time evolution
  time.dt=getConfigFromFile<double>(configFileName, "time", "dt");
  time.courantfac=getConfigFromFile<double>(configFileName, "time",
                                            "courantfac");
  time.t0=getConfigFromFile<double>(configFileName, "time", "t0");
  time.tmax=getConfigFromFile<double>(configFileName, "time", "tmax");
  time.outputinterval=getConfigFromFile<double>(configFileName, "time",
                                                "outputinterval");
  time.comparisoncount=getConfigFromFile<int>(configFileName, "time", 
                                              "comparisoncount");
  time.usefixedtimestep=getConfigFromFile<bool>(configFileName, "time",
                                                "usefixedtimestep");

  //Read in filename parameters
  file.pdesolution=getConfigFromFile<string>(configFileName, "file",
                                             "pdesolution");
  file.oneperioderror=getConfigFromFile<string>(configFileName, "file",
                                                "oneperioderror");
  file.L2error=getConfigFromFile<string>(configFileName, "file",
                                         "L2error");
  file.initialconditions=getConfigFromFile<string>(configFileName, "file",
                                                   "initialconditions");
}

//Template function that handles libconfig input from file for different
//types of data
template <typename T>
T ConfigParams::getConfigFromFile( const std::string& configFileName, 
                                   const std::string& structkey, 
                                   const std::string& key)
{
  T value;
  Config cfg;
  cfg.readFile("params.cfg");
  const Setting& root=cfg.getRoot();
  const Setting &structparams=root[structkey];
  structparams.lookupValue(key,value);
  return value;

  /*

This function returns value for a param.cfg file with a segment in the
following format.  value can be any data type but must match the
template parameter as well as the declaration of the params
sub-structure in ConfigParams.h.



structkey:
{
key=value;
key2=value2;
}

  */

}
