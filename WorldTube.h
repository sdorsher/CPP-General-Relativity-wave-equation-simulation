#ifndef WORLDTUBE_H
#define WORLDTUBE_H

#include "Grid.h"
#include "Coordinates.h"
#include "globals.h"
#include "ConfigParams.h"
#include "GridFunction.h"
#include "namespaces.h"


using namespace layers;

class WorldTube{

 public:

  WorldTube(Grid &thegrid, Coordinates &coords);

  vector<bool> addSingFieldToLeftElemExt;
  vector<bool> addSingFieldToRightElemExt;
  vector<bool> subSingFieldFromLeftElemExt;
  vector<bool> subSingFieldFromRightElemExt;
  vector<bool> inWorldTube;
  //vector<int> elem_index;
  

  void init_world_tube(Grid &thegrid, Coordinates & coords);
  
  void set_world_tube_window(Grid &thegrid, Coordinates &coords);
};


#endif
