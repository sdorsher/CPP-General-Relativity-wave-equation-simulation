#ifndef CONFIGPARAMS_H
#define CONFIGPARAMS_H

#include <fstream>
#include <libconfig.h++>
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace libconfig;

struct WaveEqParams{
  int pdenum;
  double speed;
  bool isgaussian;
  bool issinusoid;
};

struct SinusoidParams{
  double phase;
  double amp;
  double wavelength;
};

struct GaussParams{
  double amp;
  double mu;
  double sigma;
};

struct GridParams{
  double lowerlim;
  double upperlim;
  int numelems;
  int elemorder;
  bool readfromfile;
};

struct TimeParams{
  double dt;
  double courantfac;
  double t0;
  double tmax;
  double outputinterval;
  int comparisoncount;
  bool usefixedtimestep;
};

struct FileParams
{
  string pdesolution;
  string oneperioderror;
  string L2error;
  string initialconditions;
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

extern const ConfigParams params;

#endif
