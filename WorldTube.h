#ifndef WORLDTUBE_H
#define WORLDTUBE_H

#include "Grid.h"
#include "Coordinates.h"
#include "globals.h"
#include "ConfigParams.h"
#include "GridFunction.h"



class WorldTube{

  vector<bool> timeDepTrans;
  vector<bool> addSingFieldtoLeftElemExt;
  vector<bool> addSingFiledtoRightElemExt;
  vector<bool> subSingFieldFromLeftElemExt;
  vector<bool> subSingFiledFromRightElemExt;
  vector<bool> inWorldTube;
  //vector<int> elem_index;
  
  
  WorldTube(Grid thegrid, Coordinates coords);
  
  void set_world_tube_window(Grid thegrid, Coordinates coords);
}


#endif
