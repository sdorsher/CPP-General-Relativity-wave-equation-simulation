#include "TNT2.h"
#include "Grid.h"
#include "ReferenceElement.h"
#include "GridFunction.h"
#include "VectorGridFunction.h"
#include "Evolution.h"
#include "globals.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include "ConfigParams.h"
#include "DiffEq.h"
#include "TwoDVectorGridFunction.h"
#include "Modes.h"
#include "namespaces.h"
#include "Orbit.h"
#include <complex>
#include "WriteFile.h"
#include "vecMatrixTools.h"
#include "EllipticalOrbit.h"
#include "Coordinates.h"
#include "WorldTube.h"
#include "CircularOrbit.h"
#include "OutputIndices.h"

using namespace std;
using namespace layers;
//using namespace orbit;
using namespace window;
using namespace orbit;


//Initial condition options
void initialGaussian(TwoDVectorGridFunction<complex<double>>& uh, Grid& grd);
void initialSinusoid(TwoDVectorGridFunction<complex<double>>& uh, Grid& grd);
void initialSchwarzchild(TwoDVectorGridFunction<complex<double>>& uh, Grid& grd, DiffEq& eqn, vector<EffectiveSource*>& effsource);


//Characterization of convergence, error using L2 norm
double LTwoError(Grid thegrid, TwoDVectorGridFunction<complex<double>>& uh0, 
                 TwoDVectorGridFunction<complex<double>>& uhend);

int main()
{

  cout << "Start" << endl;
   //setup the modes
  Modes lmmodes(params.modes.lmax);

  Coordinates coords;

  double rmin = params.schw.p_orb / (1.0 + params.schw.ecc);
  double rmax = params.schw.p_orb / (1.0 - params.schw.ecc);
  xip = 0.5*(rmin+rmax);
  Sminus = params.hyperb.Sminus;
  rstar_orb = coords.rstar_of_r(xip, params.schw.mass);
  


  double deltar = (rstar_orb - Sminus)* 2.0/params.grid.numelems; 
  Splus = rstar_orb +round(0.5* params.grid.numelems) *deltar;
  Rminus = rstar_orb 
    - round(0.175 * params.grid.numelems) * deltar;
  Rplus = rstar_orb 
    + round(0.125 * params.grid.numelems) * deltar;
  if(params.opts.use_world_tube){
    Wminus=Rminus;
    Wplus=Rplus;
  }else{
    Wminus = Rminus + params.window.noffset * deltar;
    Wplus = Rplus - params.window.noffset * deltar;
  }
  cout << "R_star orbit" << endl;
  cout << rstar_orb << endl << endl;
  
  cout << "Sminus Rminus Rplus Splus Wminus Wplus" << endl;
  cout << Sminus << " " << Rminus << " " 
       << Rplus <<" " << Splus << " " 
       << Wminus << " " << Wplus << endl;
  cout << endl;

  
  Grid thegrid(params.grid.elemorder, params.grid.numelems, lmmodes.ntotal);
  cout << "grid established " << endl;

  EllipticalOrbit * eorb;
  CircularOrbit * corb;

  if(params.opts.use_generic_orbit){
    eorb = new EllipticalOrbit();
  }else{
    corb= new CircularOrbit();
    cout << "circular initiated" << endl;
  }


  if(params.opts.use_generic_orbit){
    
    cout << "p =" << eorb->p << endl;
    cout << "e =" << eorb->e << endl;
    cout << "chi =" << eorb->chi << endl;
    cout << "phi = " << eorb->phi << endl;
    cout << endl;
  } else{
    cout << "p =" << corb->p << endl;
    cout << "e =" << corb->e << endl;
    cout << "chi =" << corb->chi << endl;
    cout << "phi = " << corb->phi << endl;
    cout << endl;
  }
   

  
  vector<EffectiveSource*> effsource;
  if(params.opts.useSource){
    for(int i=0; i<lmmodes.ntotal; i++){
      effsource.push_back(new EffectiveSource(lmmodes.ll[i], lmmodes.mm[i], params.schw.mass));
    }
  }

  if(params.metric.schwarschild){
    double deltar = thegrid.gridNodeLocations().get(0,1)-thegrid.gridNodeLocations().get(0,0);
    R1= coords.invert_tortoise(Rminus, params.schw.mass) + 2.0*params.schw.mass;
    R2 = coords.invert_tortoise(Rplus, params.schw.mass) + 2.0* params.schw.mass;
    w1 = params.schw.p_orb-(coords.invert_tortoise(2.0*deltar, params.schw.mass)
			    +2.0*params.schw.mass)-R1;
    w2 = R2 - (params.schw.p_orb + coords.invert_tortoise(2.0*deltar, params.schw.mass)+2.0*params.schw.mass);
    nmodes = lmmodes.ntotal;


    cout << "R1 R2 w1 w2" << endl;
    cout << R1 << " " << R2 << " " << w1 <<  " " << w2 << endl << endl;
    if(params.opts.useSource){
      for(auto& x: effsource){
	x->set_window(R1,w1,1.0,1.5,R2,w2,1.0,1.5);
      }
    }
  }


  //Declaration of calculation variables and 
  //Initialization to either zero or value read from file
  //Solution to PDE, possibly a vector
  TwoDVectorGridFunction<complex<double>> uh(lmmodes.ntotal,
					     params.grid.pdenum,
                                             params.grid.numelems,
                                             params.grid.elemorder+1,
					     {0.0,0.0});

  //Initialization of initial value of calculation variables.
  TwoDVectorGridFunction<complex<double>> uh0(lmmodes.ntotal,
                                              params.grid.pdenum,
                                              params.grid.numelems,
                                              params.grid.elemorder + 1,
                                              {0.0,0.0});
  //Right hand side of PDE
  TwoDVectorGridFunction<complex<double>> RHStdvgf(lmmodes.ntotal,
						   params.grid.pdenum,
						   params.grid.numelems,
						   params.grid.elemorder + 1,
						   {0.0,0.0}); 



  OutputIndices ijoutput;

  double dx = thegrid.gridNodeLocations().get(0, 1) - thegrid.gridNodeLocations().get(0, 0);

  if(params.metric.schwarschild){
    //  int ifinite, iSplus, jfinite, jSplus;
    thegrid.find_extract_radii(rstar_orb, Splus, 
			       ijoutput, dx);
    cout << "Oribital radius and output radius in Schwarzschild coords" << endl;
    cout << params.schw.p_orb<< " " << params.grid.outputradius << endl << endl;
    cout << "Output indices for finite and scri-plus radii" << endl;
    cout << ijoutput.ifinite << " " << ijoutput.jfinite << " " 
	 << ijoutput.iSplus << " " << ijoutput.jSplus << endl << endl;
  }
  
  
  if(params.opts.use_generic_orbit){
    cout << "generic orbit set particle" << endl;
    double rp, drpdt, d2rpdt2;
    eorb->orb_of_t(coords, rp, drpdt, d2rpdt2);
    
    if(params.opts.useSource){
      //replace set_particle here! fix me! need different method for eliptical orbits
      for(auto& x: effsource){
	double rp = (eorb->p)* params.schw.mass;
	double phip = eorb->phi;
	x->set_particle(rp,phip,0.0, eorb->elip_E(),eorb->elip_L(),0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      }   //      set_particle(eorb->p,eorb->e,eorb->chi,eorb->phi,lmmodes.ntotal);
      
    }
  }
  
  
  
  
  cout << "defining the differential equation" << endl;
  //setup the differential equation
  DiffEq theequation(thegrid, lmmodes, lmmodes.ntotal, coords);

  cout << "diff eq established" << endl;

  
  WorldTube * wt;

  
  //Setup initial conditions and initialize window
  if(params.waveeq.issinusoid){
    initialSinusoid(uh, thegrid);
  } else if(params.waveeq.isgaussian) {
    initialGaussian(uh, thegrid);
  } else if(params.metric.schwarschild) {
    initialSchwarzchild(uh, thegrid, theequation, effsource);
        
    if(params.opts.use_world_tube){
      wt = new WorldTube(thegrid, coords);
      wt->init_world_tube(thegrid,coords);
      wt->set_world_tube_window(thegrid,coords); 
      cout << "use world tube " << (coords.timeDepTrans).size() << endl;
     
    }
  }


  uh0 = uh;

  double max_speed=1.0;

  if(params.opts.use_generic_orbit){
    write_fixed_time(0,params.time.t0,uh,RHStdvgf,thegrid,
		     theequation,lmmodes,true,
		     "coords",5, eorb);
     write_fixed_time(0,params.time.t0,uh,RHStdvgf,thegrid,
		      theequation,lmmodes,true,
		     "window",6, eorb);
  }else{
    cout << "circular orbit write initial " << endl;
    write_fixed_time(0,params.time.t0,uh,RHStdvgf,thegrid,
		     theequation,lmmodes,true,
		     "coords",5, corb);
    write_fixed_time(0,params.time.t0,uh,RHStdvgf,thegrid,
		     theequation,lmmodes,true,
		     "window",6, corb);
  }
  //output coords
  //  write_fixed_time(0,params.time.t0,uh,RHStdvgf,thegrid,
  //		   theequation,lmmodes,true,"coords",5);



  if(params.opts.use_generic_orbit){
    theequation.modeRHS(thegrid, uh, RHStdvgf, 0.0, true, eorb, wt, coords, max_speed, lmmodes, effsource);

  }else{
    cout << "circ write rhs initial" << endl;
    theequation.modeRHS(thegrid, uh, RHStdvgf, 0.0, true, corb, wt, coords, max_speed, lmmodes, effsource);

  }
  cout << "first call to RHS succeeded" << endl;

  

  
  
  double deltat=params.time.dt;

  
  if(params.metric.schwarschild){
    deltat = params.time.courantfac * dx/max_speed;

    cout << "deltat : " << deltat  << " " << max_speed << endl;
    cout << "dx: " << dx << endl;
    
  }else if(params.metric.flatspacetime){
    //int nt = ceil((params.time.tmax-params.time.t0) / params.time.courantfac / dx);
    //deltat = (params.time.tmax - params.time.t0) / nt;
    deltat = params.time.dt;
  }


  if (params.metric.schwarschild){
    if (params.opts.use_generic_orbit){ 
      lmmodes.sum_m_modes(uh,0.0, ijoutput.ifinite, ijoutput.jfinite, eorb);
    } else {
      lmmodes.sum_m_modes(uh,0.0, ijoutput.ifinite, ijoutput.jfinite, corb);
    }
  }

    if(params.file.outputtimefixed) {
      
      
	for(int k = 0; k < uh.TDVGFdim(); k++) {
  
	  if(params.opts.use_generic_orbit){
	    write_fixed_time(k,params.time.t0,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,
			     params.file.pdesolution,1, eorb);
	  }else{
	     write_fixed_time(k,params.time.t0,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,
			     params.file.pdesolution,1, corb);
	  }
	}//end for k
    }//end if outputfixed time
  
    
    if(params.file.outputradiusfixed){
      for(int k = 0; k < uh.TDVGFdim(); k++) {
	if(params.opts.use_generic_orbit){
	  write_fixed_radius(ijoutput,k,params.time.t0,uh,RHStdvgf,thegrid,
			     theequation,lmmodes, true,
			     params.file.fixedradiusfilename,1, eorb);
	  //write_fixed_radius(ijoutput,k,params.time.t0,uh,RHStdvgf,thegrid,
	  //		     theequation,lmmodes, true,
	  //		     "upt",2, eorb);
	}else{
	  write_fixed_radius(ijoutput,k,params.time.t0,uh,RHStdvgf,thegrid,
			     theequation,lmmodes, true,
			     params.file.fixedradiusfilename,1, corb);
	  // write_fixed_radius(ijoutput,k,params.time.t0,uh,RHStdvgf,thegrid,
	  //		     theequation,lmmodes, true,
	  //		     "upt",2, corb);
	}
      }//end for k

      if(params.opts.use_generic_orbit){
	write_summed_psi(ijoutput, params.modes.lmax,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes, true,
			 "psil",1, eorb);
	write_summed_psi(ijoutput,params.modes.lmax,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes, true,
			 "psitl",2, eorb);
	write_summed_psi(ijoutput,params.modes.lmax,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes, true,
			 "psiphil",3, eorb);
	write_summed_psi(ijoutput,params.modes.lmax,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes,true,
			 "psirl",4, eorb);
      }else{
	write_summed_psi(ijoutput,params.modes.lmax,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes, true,
			 "psil",1, corb);
	write_summed_psi(ijoutput,params.modes.lmax,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes, true,
			 "psitl",2, corb);
	write_summed_psi(ijoutput,params.modes.lmax,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes, true,
			 "psiphil",3, corb);
	write_summed_psi(ijoutput,params.modes.lmax,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes,true,
			 "psirl",4, corb);
	
	}//end params.opts.usegenericorbit
      
    }//end if outputradiusfixed
    
    //Initialize loop variables to determine when output
    //double output = deltat / 2.0;
    int outputcount =0;
    double t= params.time.t0;
  
    
  
    //BEGIN MAIN LOOP
    while(t<params.time.tmax){
    //Increment the count to determine whether or not to output
    //max_speed=1.0;
      outputcount++;
      
      //Increment the time integration
    if(params.opts.use_generic_orbit){
      rk4lowStorage(thegrid, theequation, uh, RHStdvgf, t, deltat, wt, max_speed,eorb, coords, lmmodes, effsource);
    }else{
      rk4lowStorage(thegrid, theequation, uh, RHStdvgf, t, deltat, wt, max_speed,corb,coords, lmmodes, effsource);
    }
    //Initial conditions, numerical fluxes, boundary conditions handled inside 
    //Evolution.cpp, in RHS.
    
    //increment time
    t+=deltat;
    //assert(0);
    
    //might need fill_source_all here  
    if (outputcount%params.time.outputevery == 0){
      //Output in gnuplot format
      
      for(int k = 0; k < uh.TDVGFdim(); k++) {
        if(params.file.outputtimefixed) {
	  
	  if (params.opts.use_generic_orbit){
	    write_fixed_time(k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,
			     params.file.pdesolution,1,eorb);
	    write_fixed_time(k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,"rhs",3,eorb);
	  }else{
	    write_fixed_time(k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,
			     params.file.pdesolution,1,corb);
	    write_fixed_time(k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,"rhs",3,corb);
	  /*	  write_fixed_time(ijoutput,k,t,uh,RHStdvgf,thegrid,
			   theequation,lmmodes,true,"source",2);
	  write_fixed_time(ijoutput,k,t,uh,RHStdvgf,thegrid,
			   theequation,lmmodes,true,"up",4);
	  */
	  }
	 
	}

	if(params.file.outputradiusfixed){
	  if(params.opts.use_generic_orbit){
	    write_fixed_radius(ijoutput,k,t,uh,RHStdvgf,thegrid,
			       theequation,lmmodes, true,
			       params.file.fixedradiusfilename,1,eorb);
	    // write_fixed_radius(ijoutput,k,t,uh,RHStdvgf,thegrid,
	    //		       theequation,lmmodes, true,
	    //			"upt",2,eorb);
	  }else{
	     write_fixed_radius(ijoutput,k,t,uh,RHStdvgf,thegrid,
			       theequation,lmmodes, true,
				params.file.fixedradiusfilename,1,corb);
	     //  write_fixed_radius(ijoutput,k,t,uh,RHStdvgf,thegrid,
	     //		       theequation,lmmodes, true,
	     //			 "upt",2,corb);
	  }
	    if(k==params.modes.lmax){
	    if(params.opts.use_generic_orbit){
	      lmmodes.sum_m_modes(uh, t, ijoutput.ifinite, ijoutput.jfinite,eorb);
	    }else{
	      lmmodes.sum_m_modes(uh, t, ijoutput.ifinite, ijoutput.jfinite,corb);
	    }

	    if(params.opts.use_generic_orbit){
	      
	      write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			       theequation,lmmodes,true,
			       "psil",1,eorb);
	      write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			       theequation,lmmodes,true,
			       "psitl",2,eorb);
	      write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			       theequation,lmmodes,true,
			      "psiphil",3,eorb);
	      write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,
			       "psirl",4,eorb);
	    }else{
	      	      write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			       theequation,lmmodes,true,
			       "psil",1,corb);
	      write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			       theequation,lmmodes,true,
			       "psitl",2,corb);
	      write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			       theequation,lmmodes,true,
			       "psiphil",3,corb);
	      write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,
			       "psirl",4,corb);
	    }
	  }//end if k==lmax
	}//end if outputradiusfixed
      }//end for k
     ofstream fsL2;
     fsL2.open("L2error.txt", ios::app);
     fsL2.precision(15);
     if (outputcount==params.time.comparisoncount){
       fsL2 << params.grid.elemorder << " " << params.grid.numelems << " " << deltat << " " << LTwoError(thegrid, uh0, uh) << " " << t << " " << outputcount << endl;
     }
     fsL2.close();
    }
  //Set time based on smallest grid spacing. This assumes all elements
  // are equal in size

  dx = thegrid.gridNodeLocations().get(0, 1) - thegrid.gridNodeLocations().get(0, 0);
  




  if(params.metric.schwarschild){
    deltat = params.time.courantfac * dx/max_speed;
    //deltat = params.time.dt;
  }else if(params.metric.flatspacetime){
    deltat = params.time.dt;
  }

  ofstream fstimes;
  fstimes.open("times.out",ios::app);
  fstimes.precision(15);
  fstimes << t << "\t" << deltat << "\t" << dx <<"\t" << max_speed << endl;
  fstimes.close();
  
  }//end while loop

  cout.flush();

  vector<complex<EffectiveSource*>>::size_type s{effsource.size()};
  for (vector<complex<EffectiveSource*>>::size_type i=0; i<s; i++) {
    assert(effsource.at(s-1-i)!=nullptr);
    delete effsource.at(s-1-i);
    effsource.pop_back(); 
  }
  assert(effsource.size()==0);



 if(params.opts.use_generic_orbit){
   delete eorb;
 }else{
   delete corb;
 }

 if(params.opts.use_world_tube){
   delete wt;
 }
 
}//END MAIN
                         

void initialSchwarzchild(TwoDVectorGridFunction<complex<double>>& uh, Grid& grd, DiffEq& eqn, vector<EffectiveSource*>& effsource) {
  GridFunction<double> rho(uh.GFvecDim(), uh.GFarrDim(), false);
  rho=grd.gridNodeLocations();
  ofstream fs;
  GridFunction<double> nodes = grd.gridNodeLocations();
  if(!params.opts.useSource){
    fs.open("initialdata.txt");
  }
  for(int i = 0; i < uh.GFvecDim(); i++) {
    for (int j = 0; j < uh.GFarrDim(); j++) {
      for (int n = 0; n < uh.TDVGFdim(); n++) {
	complex<double> modeval = exp(-0.5 * pow((rho.get(i,j) / params.schw.sigma), 2.0));

	
        //uh.set(n,0,i,j,0.0);
        if(!params.opts.useSource){
          uh.set(n,2,i,j,modeval); //time derivative of field (rho) (definitive answer here)

	
	  fs << setprecision(16);
	  if(n==1){
	    fs << rho.get(i,j) << " " << uh.get(0,2,i,j).real() << endl;
	  }
	}//end if
      }//end n loop modes
    }//end j loop

   
    if(params.opts.useSource){
      vector<double> rschw = grd.rschw.get(i);
      vector<double> window = grd.window.get(i);
      vector<double> dwindow = grd.dwindow.get(i);
      vector<double> d2window = grd.d2window.get(i);
      double * r = &rschw[0];
      double * win = &window[0];
      double * dwin = &dwindow[0];
      double * d2win = &d2window[0];
      auto & temp = effsource.at(0);
      temp->calc_window(params.grid.elemorder+1,r,win, dwin, d2win);
      vector<double> win3(win, win+params.grid.elemorder+1);
      grd.window.set(i,win3);
      vector<double> dwin3(dwin, dwin+params.grid.elemorder+1);
      grd.dwindow.set(i,dwin3);
      vector<double> d2win3(d2win, d2win+params.grid.elemorder+1);
      grd.d2window.set(i,d2win3);
    }
    
    //double dxmin = fabs(nodes.get(0,0)
    //			-nodes.get(0,2));
    
  }//end i loop

  if(!params.opts.useSource){
    fs.close();
  }
  
}
  

void initialSinusoid(TwoDVectorGridFunction<complex<double>>& uh, Grid& grd){
  double omega = 2.0 * PI / params.sine.wavelength;
  
  for(int k = 0; k < uh.TDVGFdim(); k++) {
    for(int i = 0; i < uh.GFvecDim(); i++){
      for (int j = 0; j < uh.GFarrDim(); j++){
        double psi = params.sine.amp * sin(omega * grd.gridNodeLocations().get(i, j)
                                           + params.sine.phase);
        double rhovar = omega * params.sine.amp * cos(omega * grd.gridNodeLocations().get(i, j)
                                                     + params.sine.phase);
        double pivar = -params.waveeq.speed * rhovar;
        //travelling waves
        //are rho and pi backward?
        uh.set(k, 0, i, j, psi);
        uh.set(k, 1, i, j, pivar);
        uh.set(k, 2, i, j, rhovar);
      }
    }
  }
}
void initialGaussian(TwoDVectorGridFunction<complex<double>>& uh, Grid& grd){
  for(int k = 0; k < uh.TDVGFdim(); k++) {
    for(int i = 0; i < uh.GFvecDim(); i++){
      for(int j = 0; j < uh.GFarrDim(); j++){
        double gaussian = params.gauss.amp * exp(-pow((grd.gridNodeLocations().get(i, j)
                                                       - params.gauss.mu), 2.0)
                                                 / 2.0 
                                                 / pow(params.gauss.sigma, 2.0));
        double dgauss = -(grd.gridNodeLocations().get(i, j) - params.gauss.mu)
          / pow(params.gauss.sigma, 2.0) * gaussian;
        uh.set(k, 0, i, j, gaussian);
        uh.set(k, 1, i, j, 0.0); //time derivative is zero
        //starts at center and splits
        uh.set(k, 2, i, j, dgauss);
      }
    }
  }
}

  double LTwoError(Grid thegrid, TwoDVectorGridFunction<complex<double>>& uh0, 
                   TwoDVectorGridFunction<complex<double>>& uhend)
{
  //FIX THIS SO IT DEALS WITH SUM OF MODES and nonrepeating waveforms
  //The square root of the integral of the squared difference
  double L2;
  L2 = 0.0;
  //was Array1D
  vector<double> weights;
  weights = thegrid.refelem.getw();
  for(int i = 0; i < uh0.GFvecDim(); i++){
    for(int j = 0; j < uh0.GFarrDim(); j++){
      double added = weights[j] 
        * pow(abs(uh0.get(0, 0, i, j)
                  - uhend.get(0, 0, i, j)), 2.0)
        / thegrid.jacobian(i);
      L2 += added;
    }
  }
  L2 = sqrt(L2); 
  return L2;
}

