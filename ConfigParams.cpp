#include "ConfigParams.h"
#include <string>

const ConfigParams params("params.cfg");

ConfigParams::ConfigParams(const std::string& configFileName)
{ 
  waveeq.pdenum=getConfigFromFile<int>(configFileName,"waveeq","pdenum");
  waveeq.speed=getConfigFromFile<double>(configFileName,"waveeq","speed");
  waveeq.isgaussian=getConfigFromFile<bool>(configFileName,"waveeq","isgaussian");
  waveeq.issinusoid=getConfigFromFile<bool>(configFileName,"waveeq","issinusoid");
  
  if(waveeq.isgaussian)
    {
      gauss.amp=getConfigFromFile<double>(configFileName,"gauss","amplitude");
      gauss.mu= getConfigFromFile<double>(configFileName,"gauss","mu");
      gauss.sigma= getConfigFromFile<double>(configFileName,"gauss","sigma");
    }
  else if(waveeq.issinusoid)
    {
      sine.amp=getConfigFromFile<double>(configFileName,"sine","amplitude");
      sine.phase=getConfigFromFile<double>(configFileName,"sine","phase");
      sine.wavelength=getConfigFromFile<double>(configFileName,"sine","wavelength");
    }else
    {
      throw invalid_argument("No initial conditions set in parameter file.");
    }

  grid.lowerlim=getConfigFromFile<double>(configFileName,"grid","lowerlim");
  grid.upperlim=getConfigFromFile<double>(configFileName,"grid","upperlim");
  grid.numelems=getConfigFromFile<int>(configFileName,"grid","numelems");
  grid.elemorder=getConfigFromFile<int>(configFileName,"grid","elemorder");
  grid.readfromfile=getConfigFromFile<bool>(configFileName,"grid","readfromfile");
  
  time.dt=getConfigFromFile<double>(configFileName,"time","dt");
  time.courantfac=getConfigFromFile<double>(configFileName,"time","courantfac");
  time.t0=getConfigFromFile<double>(configFileName,"time","t0");
  time.tmax=getConfigFromFile<double>(configFileName,"time","tmax");
  time.outputinterval=getConfigFromFile<double>(configFileName,"time","outputinterval");
  time.comparisoncount=getConfigFromFile<int>(configFileName,"time","comparisoncount");
  time.usefixedtimestep=getConfigFromFile<int>(configFileName,"time","usefixedtimestep");


  file.pdesolution=getConfigFromFile<string>(configFileName,"file","pdesolution");
  file.oneperioderror=getConfigFromFile<string>(configFileName,"file","oneperioderror");
  file.L2error=getConfigFromFile<string>(configFileName,"file","L2error");
  file.initialconditions=getConfigFromFile<string>(configFileName,"file","initialconditions");
}

template <typename T>
T ConfigParams::getConfigFromFile( const std::string& configFileName, const std::string& structkey, const std::string& key)
{
  T value;
  Config cfg;
  cfg.readFile("params.cfg");
  const Setting& root=cfg.getRoot();
  const Setting &structparams=root[structkey];
  structparams.lookupValue(key,value);
  return value;
}
