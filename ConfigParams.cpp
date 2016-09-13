#include "ConfigParams.h"
#include <string>


//see .h file for definitions of variables.


const ConfigParams params("params.cfg");

ConfigParams::ConfigParams(const std::string& configFileName)
{ //Use template function that uses libconfig commands to read in global
  //parameters into a params structure, for example, params.waveq.pdenum
  
  //read in options
  opts.useSource = getConfigFromFile<bool>(configFileName, "options", 
					   "useSource");
  opts.L2outputcount = getConfigFromFile<int>(configFileName, "options", 
					     "Ltwooutputcount");
  cout << opts.L2outputcount << endl;
  opts.turn_on_source_smoothly = getConfigFromFile<bool>(configFileName,
							 "options",
							 "turn_on_source_smoothly");
  opts.use_generic_orbit=getConfigFromFile<bool>(configFileName, "options",
						 "use_generic_orbit");
  opts.use_world_tube=getConfigFromFile<bool>(configFileName, "options",
						 "use_world_tube");							

  //Read parameters associated with the metric
  metric.flatspacetime = getConfigFromFile<bool>(configFileName, "metric", 
                                                 "flatspacetime");
  metric.schwarschild =  getConfigFromFile<bool>(configFileName, "metric", 
                                                 "schwarschild");
  //if schwartschild, also hyperboloidal

  if(metric.flatspacetime) {
    //Read parameters associated with the wave equation
    waveeq.speed=getConfigFromFile<double>(configFileName, "waveeq", "speed");
    waveeq.isgaussian=getConfigFromFile<bool>(configFileName,"waveeq",
                                              "isgaussian");
    waveeq.issinusoid=getConfigFromFile<bool>(configFileName, "waveeq",
                                              "issinusoid");
    grid.lowerlim=getConfigFromFile<double>(configFileName, "grid", "lowerlim");
    grid.upperlim=getConfigFromFile<double>(configFileName, "grid", "upperlim");

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
  } else { //if metric is not flat (Schwarzschild)
    schw.mass = getConfigFromFile<double>(configFileName, "schw", "mass");
    schw.sigma = getConfigFromFile<double>(configFileName, "schw", "sigma");
    schw.p_orb = getConfigFromFile<double>(configFileName, "schw", "p_orb");
    schw.ecc = getConfigFromFile<double>(configFileName, "schw", "ecc");
    
    window.noffset = getConfigFromFile<int>(configFileName, "window", "noffset");
    timewindow.tsigma = getConfigFromFile<double>(configFileName, "timewindow",
						  "tsigma");
    timewindow.torder = getConfigFromFile<int>(configFileName, "timewindow",
					       "torder");

  //Read in parameters associated with hyperboloidal coordinates
    hyperb.Sminus = getConfigFromFile<double>(configFileName, "hyperb",
					      "Sminus");
    hyperb.Splus = getConfigFromFile<double>(configFileName, "hyperb", "Splus");
    hyperb.Rplus = getConfigFromFile<double>(configFileName, "hyperb", "Rplus");
    hyperb.Rminus = getConfigFromFile<double>(configFileName, "hyperb",
					      "Rminus");

    grid.lowerlim=hyperb.Sminus; 
    grid.upperlim=hyperb.Splus;
    //in this case, the boundaries on the grid should correspond to the horizon
    //and Scri-plus

  }


  
  //Read in parameters associated with the modes
  modes.lmax = getConfigFromFile<int>(configFileName,"modes", "lmax");

  //Read in parameters associated with the grid
  grid.Adim = getConfigFromFile<int>(configFileName,"grid","Adim");
  grid.Ddim = getConfigFromFile<int>(configFileName,"grid","Ddim");
  
  grid.pdenum=getConfigFromFile<int>(configFileName, "grid", "pdenum");

    grid.numelems=getConfigFromFile<int>(configFileName, "grid", "numelems");
    grid.elemorder=getConfigFromFile<int>(configFileName, "grid", "elemorder");
    grid.readfromfile=getConfigFromFile<bool>(configFileName, "grid",
                                             "readfromfile");
    grid.outputradius = getConfigFromFile<double>(configFileName, "grid", 
                                                  "outputradius");
   
  
  //Read in parameters associated with time evolution
    time.dt=getConfigFromFile<double>(configFileName, "time", "dt");
    time.courantfac=getConfigFromFile<double>(configFileName, "time",
                                              "courantfac");
    time.t0=getConfigFromFile<double>(configFileName, "time", "t0");
    time.tmax=getConfigFromFile<double>(configFileName, "time", "tmax");
    time.comparisoncount=getConfigFromFile<int>(configFileName, "time", 
                                                "comparisoncount");
    time.outputevery=getConfigFromFile<int>(configFileName, "time",
                                                  "outputevery");

    //Read in filename parameters

    file.outputtimefixed =getConfigFromFile<bool>(configFileName, "file",
                                                    "outputtimefixed");

    file.outputradiusfixed = getConfigFromFile<bool>(configFileName, "file",
                                                       "outputradiusfixed");

    file.fixedradiusfilename = getConfigFromFile<string>(configFileName, "file",
                                                         "fixedradiusfilename");


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
				   const char *structkey, 
				   const char *key)
{
  T value;
  Config cfg;
  cfg.readFile("params.cfg");
  const Setting& root=cfg.getRoot();
  const Setting& structparams=root[structkey];
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
