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
#include "orbit.h"
#include <complex>
#include "source_interface.h"
#include "WriteFile.h"
#include "vecMatrixTools.h"
#include "EllipticalOrbit.h"
#include "Coordinates.h"
#include "WorldTube.h"
#include "CircularOrbit.h"


using namespace std;
using namespace layers;
using namespace orbit;
using namespace window;
using namespace source_interface;


//Initial condition options
void initialGaussian(TwoDVectorGridFunction<complex<double>>& uh, Grid& grd);
void initialSinusoid(TwoDVectorGridFunction<complex<double>>& uh, Grid& grd);
void initialSchwarzchild(TwoDVectorGridFunction<complex<double>>& uh, Grid& grd, DiffEq& eqn);


//Characterization of convergence, error using L2 norm
double LTwoError(Grid thegrid, TwoDVectorGridFunction<complex<double>>& uh0, 
                 TwoDVectorGridFunction<complex<double>>& uhend);

int main()
{

  cout << "Start" << endl;
  Grid thegrid(params.grid.elemorder, params.grid.numelems, lmmodes.ntotal,
               lowlim, uplim);
  cout << "grid established " << endl;
   //setup the modes
  Modes lmmodes(params.modes.lmax);

  Orbit orb;
  
  if(params.opts.use_generic_orbit){
    orb = new EllipticalOrbit();
  }else{
    orb= new CircularOrbit();
  }
  
  if(params.opts.useSource) {
    init_source( lmmodes, params.schw.mass);
  }

  Coordinates coords();
  
  
  if(params.opts.use_generic_orbit){
    orb.orb_of_t(rp,drpdt,d2rpdt2);
    cout << "rp = " << rp << " drpdt = " << drpt << " d2rpdt2 = " d2rpdt2 << endl;
    if(params.opts.useSource){
      set_particle(p,e,chi,phi,nmodes);
    }
  }
  
    
  //find the indices associated with the radii to extract the solution at
  OutputIndices ijoutput;
  if(params.metric.schwarschild){
    //  int ifinite, iSplus, jfinite, jSplus;
    thegrid.find_extract_radii(coordobj.rstar_of_r(params.grid.outputradius,
					  params.schw.mass), Splus, 
			       ijoutput);
    cout << "Oribital radius and output radius in Schwarzschild coords" << endl;
    cout << params.schw.p_orb<< " " << params.grid.outputradius << endl << endl;
    cout << "Output indices for finite and scri-plus radii" << endl;
    cout << ijoutput.ifinite << " " << ijoutput.jfinite << " " 
	 << ijoutput.iSplus << " " << ijoutput.jSplus << endl << endl;
  }
  
  
  cout << "defining the differential equation" << endl;
  //setup the differential equation
  DiffEq theequation(thegrid, lmmodes, lmmodes.ntotal, coordobj);

  cout << "diff eq established" << endl;

  
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

  
  //Setup initial conditions and initialize window
  if(params.waveeq.issinusoid){
    initialSinusoid(uh, thegrid);
  } else if(params.waveeq.isgaussian) {
    initialGaussian(uh, thegrid);
  } else if(params.metric.schwarschild) {
    initialSchwarzchild(uh, thegrid, theequation);

    WorldTube worldtb;
    
    if(params.opts.use_world_tube){
      worldtb = new WorldTube(thegrid, coords);
      worldtb.set_world_tube_window(thegrid,coords);
    }
    
  }

  //output window function
  ofstream fs;
  fs.open("window.txt");
  for(int i=0; i<params.grid.numelems; i++) {
    for(int j = 0; j<params.grid.elemorder+1; j++) {
      fs << thegrid.gridNodeLocations().get(i,j) << " " 
         << theequation.window.get(i,j) << " " << theequation.dwindow.get(i,j) << " "
         << theequation.d2window.get(i,j) << endl;
    }
  }
  cout << "window output" << endl;

  uh0 = uh;



  //output coords
  //  write_fixed_time(0,params.time.t0,uh,RHStdvgf,thegrid,
  //		   theequation,lmmodes,true,"coords",5);
  
  theequation.modeRHS(thegrid, uh, RHStdvgf, 0.0, true);

  cout << "first call to RHS succeeded" << endl;
  
  double deltat, max_speed;

  double dx = thegrid.gridNodeLocations().get(0, 1) - thegrid.gridNodeLocations().get(0, 0);
  

  if(params.metric.schwarschild){
    deltat = params.time.courantfac * dx/max_speed;
    //deltat = params.time.dt;
    cout << "set and actual time step, based on courant factor" << endl;
    
    //temporary
    cout << dx << " " << deltat << endl << endl;
  }else if(params.metric.flatspacetime){
    //int nt = ceil((params.time.tmax-params.time.t0) / params.time.courantfac / dx);
    //deltat = (params.time.tmax - params.time.t0) / nt;
    deltat = params.time.dt;
    cout << "deltat set to dt" << endl;
    cout << dx << " " << deltat << endl;
  }

  
  if (params.metric.schwarschild){
    lmmodes.sum_m_modes(uh,0.0, ijoutput.ifinite, ijoutput.jfinite);
  }

  
  for(int k = 0; k < uh.TDVGFdim(); k++) {
    if(params.file.outputtimefixed) {


      write_fixed_time(k,params.time.t0,uh,RHStdvgf,thegrid,
		       theequation,lmmodes,true,
		       params.file.pdesolution,1);
      write_fixed_time(k,params.time.t0,uh,RHStdvgf,thegrid,
		       theequation,lmmodes,true,"rhs",3);
	     /*write_fixed_time(k,params.time.t0,uh,RHStdvgf,thegrid,
		       theequation,lmmodes,true,"source",2);
     write_fixed_time(k,params.time.t0,uh,RHStdvgf,thegrid,
		       theequation,lmmodes,true,"up",4);
	     */
    }

    if(params.file.outputradiusfixed){
      write_fixed_radius(ijoutput,k,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes, true,
			 params.file.fixedradiusfilename,1);

      if(k==params.modes.lmax){
	write_summed_psi(ijoutput,k,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes,true,
			 "psil",1);
	write_summed_psi(ijoutput,k,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes,true,
			 "psitl",2);
	write_summed_psi(ijoutput,k,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes,true,
			 "psiphil",3);
	write_summed_psi(ijoutput,k,params.time.t0,uh,RHStdvgf,thegrid,
			 theequation,lmmodes,true,
			 "psirl",4);
      }//end if k==lmax
    }//end if outputradiusfixed
  }//end for k
  
  //Initialize loop variables to determine when output
  //double output = deltat / 2.0;
  int outputcount =0;
  double t= params.time.t0;

  
  //  for(double t = params.time.t0; t < params.time.tmax + deltat; t += deltat) {
  while(t<params.time.tmax){
    //Increment the count to determine whether or not to output
    outputcount++;
    
    //Increment the time integration
    rk4lowStorage(thegrid, theequation, uh, RHStdvgf, t, deltat, max_speed);
    //Initial conditions, numerical fluxes, boundary conditions handled inside 
    //Evolution.cpp, in RHS.

    //increment time
    t+=deltat;
    
    //might need fill_source_all here  
    if (outputcount%params.time.outputevery == 0){
      //Output in gnuplot format
    
     for(int k = 0; k < uh.TDVGFdim(); k++) {
        if(params.file.outputtimefixed) {

	  write_fixed_time(k,t,uh,RHStdvgf,thegrid,
			   theequation,lmmodes,true,
			   params.file.pdesolution,1);
	  write_fixed_time(k,t,uh,RHStdvgf,thegrid,
			   theequation,lmmodes,true,"rhs",3);
	  /*	  write_fixed_time(ijoutput,k,t,uh,RHStdvgf,thegrid,
			   theequation,lmmodes,true,"source",2);
	  write_fixed_time(ijoutput,k,t,uh,RHStdvgf,thegrid,
			   theequation,lmmodes,true,"up",4);
	  */
	  
	 
	}
      
	if(params.file.outputradiusfixed){
	  write_fixed_radius(ijoutput,k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes, true,
			   params.file.fixedradiusfilename,1);
	  if(k==params.modes.lmax){
	    lmmodes.sum_m_modes(uh, t, ijoutput.ifinite, ijoutput.jfinite);
	    write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,
			     "psil",1);
	    write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,
			     "psitl",2);
	    write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,
			     "psiphil",3);
	    write_summed_psi(ijoutput,k,t,uh,RHStdvgf,thegrid,
			     theequation,lmmodes,true,
			     "psirl",4);
	  }//end if k==lmax
	}//end if outputradiusfixed
      }//end for k
     ofstream fsL2;
     fsL2.open("L2error.txt", ios::app);
     fsL2.precision(15);
     if (outputcount==params.opts.L2outputcount){
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
    cout << "set and actual time step, based on courant factor" << endl;
    
    //temporary
    cout << dx << " " << deltat << endl << endl;
  }else if(params.metric.flatspacetime){
    //int nt = ceil((params.time.tmax-params.time.t0) / params.time.courantfac / dx);
    //deltat = (params.time.tmax - params.time.t0) / nt;
    deltat = params.time.dt;
    cout << "deltat set to dt" << endl;
    cout << dx << " " << deltat << endl;
  }

  }//end while loop

  cout.flush();
  
clean_source();

 delete orb;

 if(params.opts.use_world_tube){
   delete worldtb;
 }
 
}
                         

void initialSchwarzchild(TwoDVectorGridFunction<complex<double>>& uh, Grid& grd, DiffEq& eqn) {
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
          uh.set(n,2,i,j,modeval);

	
	//}else{
	  //uh.set(n,2,i,j,0.0);
	  //}
        
	  //uh.set(n,1,i,j,0.0);
	  fs << setprecision(16);
	  fs << rho.get(i,j) << " " << uh.get(0,2,i,j).real() << endl;
	}//end if
      }//end n loop modes
    }//end j loop

   
    if(params.opts.useSource){
      vector<double> rschw = grd.rschw.get(i);
      vector<double> window = eqn.window.get(i);
      vector<double> dwindow = eqn.dwindow.get(i);
      vector<double> d2window = eqn.d2window.get(i);
      double * r = &rschw[0];
      double * win = &window[0];
      double * dwin = &dwindow[0];
      double * d2win = &d2window[0];
      calc_window(params.grid.elemorder+1,r,win, dwin, d2win);
      vector<double> win3(win, win+params.grid.elemorder+1);
      eqn.window.set(i,win3);
      vector<double> dwin3(dwin, dwin+params.grid.elemorder+1);
      eqn.dwindow.set(i,dwin3);
      vector<double> d2win3(d2win, d2win+params.grid.elemorder+1);
      eqn.d2window.set(i,d2win3);
    }
    
    double dxmin = fabs(nodes.get(0,0)
			-nodes.get(0,2));
    
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

