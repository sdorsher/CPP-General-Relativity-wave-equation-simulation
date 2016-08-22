#include "WorldTube.h"

WorldTube(Grid thegrid, Coordinates coords),timeDepTrans{params.grid.numelems,false},addSingFieldToLeft{params.grid.numelems,0}
{
  for (int j=1; j<params.grid.numelem; j++){
    double rho = thegrid.gridNodeLocations().get(j,0);
    if((abs(rho-coords.Rminus)<1.e-10)||(abs(rho-coords.Rplus)<1.e-10))
      {
	timeDepTrans.at(j)=true;
      }
    
    if((rho>coords.Wminus)&&(rho<coords.Wplus)){
      inWorldTube.at(j)=true;
    }else if(abs(rho-coords.Wminus)<1.e-10){
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
	theeq.window.set(j,i,0.);
	theeq.dwindow.set(j,i,0.);
	theeq.d2window.set(j,i,0.);
      }
      else if((i==params.grid.elemorder)&&(j=params.grid.numelems-1)){
	theeq.window.set(j,i,1.);
	theeq.dwindow.set(j,i,1.);
	theeq.d2window.set(j,i,1.);
      }else if((!inWorldTube.at(j))&(inWorldTube.at(j+1))&&(i==params.grid.elemorder)){
	theeq.window.set(j,i,1.);
	theeq.dwindow.set(j,i,1.);
	theeq.d2window.set(j,i,1.);
      } else if ((!inWorldTube.at(j))&&(inWorldTube.at(j-1))&&(i==0)){
  	theeq.window.set(j,i,1.);
	theeq.dwindow.set(j,i,1.);
	theeq.d2window.set(j,i,1.);
      } else if (inWorldTube.at(j)){
	theeq.window.set(j,i,1.);
	theeq.dwindow.set(j,i,1.);
	theeq.d2window.set(j,i,1.);
      } else {
	theeq.window.set(j,i,0.);
	theeq.dwindow.set(j,i,0.);
	theeq.d2window.set(j,i,0.);
      }
    }
  }
}
 
