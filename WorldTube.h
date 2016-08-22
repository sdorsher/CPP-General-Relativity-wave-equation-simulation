#ifndef WORLDTUBE_H
#define WORLDTUBE_H

#include "Grid.h"
#include "Coordinates.h"
#include "globals.h"
#include "ConfigParams.h"
#include "GridFunction.h"

enum PMzero{-1,0,1};


class WorldTube{

  vector<bool> timeDepTrans;
  vector<PMzero> addSingFieldtoLeft;
  vector<bool> inWorldTube;
  
  WorldTube(Grid thegrid, Coordinates coords);
  
  void set_world_tube_window(Grid thegrid, Coordinates coords);
}


#endif
