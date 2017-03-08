#include "WorldTube.h"

WorldTube::WorldTube(Grid& thegrid, Coordinates &coords){}

void WorldTube::init_world_tube(Grid & thegrid, Coordinates &coords)
{
  addSingFieldToLeftElemExt.resize(params.grid.numelems-1,false);
  addSingFieldToRightElemExt.resize(params.grid.numelems-1,false);
  subSingFieldFromLeftElemExt.resize(params.grid.numelems-1,false);
  subSingFieldFromRightElemExt.resize(params.grid.numelems-1,false);
  inWorldTube.resize(params.grid.numelems-1,false);
  //boolean flags represent boundaries excluding scri+ and scri-

  
  for (int j=0; j<params.grid.numelems-1; j++){
    double rho = thegrid.gridNodeLocations().get(j+1,0);
    if((rho>=Rminus)&&(rho<=Rplus))
      {
	coords.timeDepTrans.at(j)=true;
      }
    if((rho>=Wminus)&&(rho<Wplus)){
      inWorldTube.at(j)=true;
    }
    if(abs(rho-Wminus)<1.e-10){
      cout << "Wminus " << j << endl;
      addSingFieldToLeftElemExt.at(j+1)=true;//position 2, in world tube at left
      subSingFieldFromRightElemExt.at(j)=true;//position 1, out of world tube at left
    }else if(abs(rho-Wplus)<1.e-10){
      cout << "Wplus " << j << endl;
      addSingFieldToRightElemExt.at(j)=true;//position 4, out of world tube at right
      subSingFieldFromLeftElemExt.at(j+1)=true;//position 3, in world tube at right
    }
  }
  for (int j=0; j<params.grid.numelems-1; j++){
  if((!coords.timeDepTrans.at(j))&&(coords.timeDepTrans.at(j+1))){
     coords.timeDepTrans.at(j)=true;
   }
  }
  
}

void WorldTube::set_world_tube_window(Grid &thegrid, Coordinates &coords){

  for(int j=0; j<params.grid.numelems; j++){
    for(int i=0; i<=params.grid.elemorder; i++){
    
      if(j==0){
	thegrid.window.set(j,i,0.);
	thegrid.dwindow.set(j,i,0.);
	thegrid.d2window.set(j,i,0.);
      }
      else if(j==params.grid.numelems-1){
 	thegrid.window.set(j,i,0.);
	thegrid.dwindow.set(j,i,0.);
	thegrid.d2window.set(j,i,0.);
      } else if (inWorldTube.at(j-1)||(inWorldTube.at(j))){
	thegrid.window.set(j,i,1.);
	thegrid.dwindow.set(j,i,0.);
	thegrid.d2window.set(j,i,0.);
      } else {
	thegrid.window.set(j,i,0.);
	thegrid.dwindow.set(j,i,0.);
	thegrid.d2window.set(j,i,0.);
      }
    }
  }
}
 
