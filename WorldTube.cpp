#include "WorldTube.h"

WorldTube(Grid thegrid, Coordinates coords),timeDepTrans{params.grid.numelems,false},addSingFieldToLeft{params.grid.numelems,0}
{
  for (int j=1; j<params.grid.numelem; j++){
    double rho = thegrid.gridNodeLocations().get(j,0);
    if((abs(rho-coords.Rminus)<1.e-10)||(abs(rho-coords.Rplus)<1.e-10))
      {
	timeDepTrans.at(j)=true;
      }
    if(abs(rho-coords.Wminus)<1.e-10){
      addSingFieldToLeft.at(j)=1;
    }else if(abs(rho-coords.Wplus)<1.e-10){
      addSingFieldToLeft.at(j)=-1;
    }
  }
      
}


void WorldTube::set_world_tube_window(Grid thegrid, Coordinates coords, DiffEq theeq){

  double wind, dwind, d2wind;
  
  for(int i=0; i<=params.grid.elemorder; i++){
    for(int j=0; j<params.grid.numelems; j++){
      if((i==0)&&(j==0)){
	wind=0.;
	dwind=0.;
	d2wind=0.;
      } else if((i==params.grid.elemorder)&&(j==params.grid.numelems-1)){
	wind=0.;
	dwind=0.;
	d2wind=0.;
      } else if((timeDepTrans.at(j))&&(timeDepTrans.at(j+1))){
	wind=1.;
	dwind=1.;
	d2wind=1.;
      } else {
	wind=0.;
	dwind=0.;
	d2wind=0.;
      }
      theeq.window.set(j,i,wind);
      theeq.dwindow.set(j,i,dwind);
      theeq.d2window.set(j,i,d2wind);
    }
  }
}
 
