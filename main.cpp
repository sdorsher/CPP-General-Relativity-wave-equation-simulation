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

using namespace std;
using namespace layers;
using namespace orbit;
using namespace window;
using namespace source_interface;


//Initial condition options
void initialGaussian(TwoDVectorGridFunction<complex<double>>& uh, Grid grd);
void initialSinusoid(TwoDVectorGridFunction<complex<double>>& uh, Grid grd);
void initialSchwarzchild(TwoDVectorGridFunction<complex<double>>& uh, Grid grd);


//Characterization of convergence, error using L2 norm
double LTwoError(Grid thegrid, TwoDVectorGridFunction<complex<double>>& uh0, 
                 TwoDVectorGridFunction<complex<double>>& uhend);

int main()
{

  //params is initialized at the beginning of this file by its 
  //inclusion in the file
  //modify such that hyperboloidal and tortoise boundaries
  // occur at element boundaries

  //setup the modes
  Modes lmmodes(params.modes.lmax);

  double rmin = params.schw.p_orb / (1.0 + params.schw.ecc);
  double rmax = params.schw.p_orb / (1.0 - params.schw.ecc);
  double xip = 0.5*(rmin+rmax);
  Sminus = params.hyperb.Sminus;
  double rstar_orb = rstar_of_r(xip, params.schw.mass);

  double deltar = (rstar_orb - Sminus)* 2.0/params.grid.numelems; 
  Splus = rstar_orb +round(0.5* params.grid.numelems) *deltar;
  Rminus = rstar_orb 
    - round(0.175 * params.grid.numelems) * deltar;
  Rplus = rstar_orb 
    + round(0.125 * params.grid.numelems) * deltar;
  Wminus = Rminus + params.window.noffset * deltar;
  Wplus = Rplus - params.window.noffset * deltar;

  cout << "R_star orbit" << endl;
  cout << rstar_orb << endl << endl;

  cout << "Sminus Rminus Rplus Splus Wminus Wplus" << endl;
  cout << Sminus << " " << Rminus << " " 
       << Rplus <<" " << Splus << " " 
       << Wminus << " " << Wplus << endl;
  cout << endl;

  //initialize orbit
  initialize_orbit();
  cout << "p =" << p << endl;
  cout << "e =" << e << endl;
  cout << "chi =" << chi << endl;
  cout << "phi = " << phi << endl;
  cout << endl;

 
  if(params.opts.useSource) {
    init_source( lmmodes, params.schw.mass);
  }
  R1= invert_tortoise(Rminus, params.schw.mass) + 2.0*params.schw.mass;
  R2 = invert_tortoise(Rplus, params.schw.mass) + 2.0* params.schw.mass;
  w1 = params.schw.p_orb-(invert_tortoise(2.0*deltar, params.schw.mass)
                          +2.0*params.schw.mass)-R1;
  w2 = R2 - (params.schw.p_orb + invert_tortoise(2.0*deltar, params.schw.mass)+2.0*params.schw.mass);
  nmodes = lmmodes.ntotal;

  cout << "R1 R2 w1 w2" << endl;
  cout << R1 << " " << R2 << " " << w1 <<  " " << w2 << endl << endl;
  
  if(params.opts.useSource) {
    set_window(R1, w1, 1.0, 1.5, R2, w2, 1.0, 1.5, lmmodes.ntotal);
  }
  
  
  double lowlim, uplim; 
  
  //setup the grid and the reference element
  if (params.metric.flatspacetime) {
    lowlim = params.grid.lowerlim;
    uplim = params.grid.upperlim;
    
   } else if (params.metric.schwarschild) {
    lowlim = Sminus;
    uplim = Splus;
  }
 
  Grid thegrid(params.grid.elemorder, params.grid.numelems, lmmodes.ntotal,
               lowlim, uplim);

  cout << "grid established " << endl;
  

  //find the indices associated with the radii to extract the solution at
  
  int ifinite, iSplus, jfinite, jSplus;
  thegrid.find_extract_radii(rstar_of_r(params.grid.outputradius,
                                        params.schw.mass), Splus, 
                             ifinite,iSplus, jfinite, jSplus);
  
  cout << "Oribital radius and output radius in Schwarzschild coords" << endl;
  cout << params.schw.p_orb<< " " << params.grid.outputradius << endl << endl;
  cout << "Output indices for finite and scri-plus radii" << endl;
  cout << ifinite << " " << jfinite << " " << iSplus << " " << jSplus << endl << endl;

  GridFunction<double> nodes = thegrid.gridNodeLocations();
  
  
  
  //setup the differential equation
  DiffEq theequation(thegrid, lmmodes, lmmodes.ntotal);

  //Declaration of calculation variables and 
  //Initialization to either zero or value read from file
  TwoDVectorGridFunction<complex<double>> uh(lmmodes.ntotal,
                                   params.grid.pdenum,
                                             params.grid.numelems,
                                             params.grid.elemorder+1,
                                0.0); 
  TwoDVectorGridFunction<complex<double>> uh0(lmmodes.ntotal,
                                              params.grid.pdenum,
                                              params.grid.numelems,
                                              params.grid.elemorder + 1,
                                              0.0);
  //Solution to PDE, possibly a vector 
  TwoDVectorGridFunction<complex<double>> RHStdvgf(lmmodes.ntotal,
                                       params.grid.pdenum,
                                                 params.grid.numelems,
                                                 params.grid.elemorder + 1,
                                    0.0); //right hand side of PDE
  

  //Setup initial conditions and initialize window
  if(params.waveeq.issinusoid){
    initialSinusoid(uh, thegrid);
  } else if(params.waveeq.isgaussian) {
    initialGaussian(uh, thegrid);
  } else if(params.metric.schwarschild) {
    initialSchwarzchild(uh, thegrid);
    //need to write initial swcharzchild
  }
          
  //output window function
  ofstream fs;
  fs.open("window.txt");
  for(int i=0; i<params.grid.numelems; i++) {
    for(int j = 0; j<params.grid.elemorder+1; j++) {
      fs << thegrid.gridNodeLocations().get(i,j) << " " 
         << thegrid.window.get(i,j) << " " << thegrid.dwindow.get(i,j) << " "
         << thegrid.d2window.get(i,j) << endl;
    }
  }
  cout << "window output" << endl;

  uh0 = uh;


  //Set time based on smallest grid spacing
  double dt0 = nodes.get(0, 1) - nodes.get(0, 0);

  int nt = ceil(params.time.tmax / params.time.courantfac / dt0);

  double deltat;
  if(params.time.usefixedtimestep){
    deltat = params.time.dt;
  } else {
    deltat = (params.time.tmax - params.time.t0) / nt; 
    //Make deltat go into tmax an integer number of times
  }
  cout << "set and actual time step, based on courant factor" << endl;
  cout << dt0 << " " << deltat << endl << endl;


  theequation.modeRHS(thegrid, uh, RHStdvgf, 0.0, false);

  

  //Initialize loop variables to determine when output
  //double output = deltat / 2.0;
  int outputcount =0;
  for(double t = params.time.t0; t < params.time.tmax + deltat; t += deltat) {
    //    if(output > 0.0){
    if (outputcount%params.time.outputevery == 0){
      //Output in gnuplot format
      for(int k = 0; k < uh.modesDim(); k++) {
        if(params.file.outputtimefixed) {

          
          ofstream fs;
          
          string solnfilestring;
          ostringstream oss;
          oss << params.file.pdesolution << "." << k << ".txt";
          fs.open(oss.str(), ios::app);
          fs << endl << endl;
          fs << " #time = " << t << endl;
          for (int i = 0; i < uh.gridDim(); i++){
            for(int j = 0; j < uh.pointsDim(); j++){
              //Print out at select time steps
              fs << nodes.get(i, j) << " "
                 << uh.get(k, 0, i, j).real() << " " 
                 << uh.get(k, 1, i, j).real() <<" " 
                 << uh.get(k, 2, i, j).real()<< endl;
            }
          }
            fs.close();

	    //HERE
	    //	    for(int k = 0; k< uh.modesDim(); k++){
	    ofstream fs5;
	    ofstream fs6;
	    ostringstream oss5;
	    ostringstream oss6;
	    oss5 << "source" << "." << k << ".txt";
	    oss6 << "rhs" << "." << k << ".txt";
	    fs5.open(oss5.str(), ios::app);
	    fs6.open(oss6.str(), ios::app);
	    fs5 << endl << endl;
	    fs5 << " #time = " << t << endl;
	    fs6 << endl << endl;
	    fs6 << " #time = " << t << endl;
	    for (int i = 0; i < uh.gridDim(); i++){
	      for(int j = 0; j < uh.pointsDim(); j++){
		//Print out at select time steps
		fs5 << thegrid.gridNodeLocations().get(i, j) << " "
		    << thegrid.source.get(k, i, j).real() << " " 
		    << thegrid.source.get(k, i, j).imag() << endl; 
		fs6 << thegrid.gridNodeLocations().get(i,j) << " "
		    << RHStdvgf.get(k,0,i,j).real() << " "
		    << RHStdvgf.get(k,1,i,j).real() << " "
		    << RHStdvgf.get(k,2,i,j).real() << " "
		    << RHStdvgf.get(k,0,i,j).imag() << " "
		    << RHStdvgf.get(k,1,i,j).imag() << " "
		    << RHStdvgf.get(k,2,i,j).imag() << endl;
	      }
	    }
	    fs5.close();
	    fs6.close();
	    //	}
	


	
	}
      




	if(params.file.outputradiusfixed){
          ofstream fs;
          ostringstream oss;
          oss << params.file.fixedradiusfilename << "." << k << ".txt";
          fs.open(oss.str(), ios::app);
          fs << nodes.get(ifinite, jfinite) << " " 
             << t << " "
             << uh.get(k, 0, ifinite, jfinite).real() << " " 
             << uh.get(k, 1, ifinite, jfinite).real() <<" " 
             << uh.get(k, 2, ifinite, jfinite).real()<< " " 
             << nodes.get(iSplus, jSplus) << " " 
             << uh.get(k, 0, iSplus, jSplus).real() << " " 
             << uh.get(k, 1, iSplus, jSplus).real() <<" " 
             << uh.get(k, 2, iSplus, jSplus).real()<< endl;
          fs.close();
        }
      }//end for
    }else{
    }
    
      //Output the difference in the waveforms between the 
      //oscillation initially and after one period
      /*
      if(outputcount == params.time.comparisoncount){
        cout << t << endl;
        ofstream fs2;
        fs2.open(params.file.oneperioderror);
        for(int i = 0; i < uh.gridDim(); i++){
          for (int j = 0; j < uh.pointsDim(); j++){
            fs2 << thegrid.gridNodeLocations().get(i, j) << " " 
                << uh.get(0, 0, i, j) - uh0.get(0, 0, i, j) << endl;
          }
        }
        fs2.close();
      */
        //Append the L2 error to that file, measured after one period
        /*ofstream fsconvergence;
        fsconvergence.open(params.file.L2error,ios::app);
        
        double L2;
        L2=LTwoerror(thegrid, uh0, uh);
        cout << "Order, deltat, num elems, L2 norm" << endl;
        cout << params.grid.elemorder << " " << deltat << " " 
             << params.grid.numelems << " " << L2 << endl;
        fsconvergence << params.grid.elemorder << " " << deltat 
                      << " " << params.grid.numelems << " " << L2 << endl;
        fsconvergence.close();
       
        }*/
    

    //    cout << "outputcount = " <<outputcount << endl;
    //Increment the timestep
    rk4lowStorage(thegrid, theequation, uh, RHStdvgf, t, deltat);
    //Initial conditions, numerical fluxes, boundary conditions handled inside 
    //Evolution.cpp, in RHS.
    
    //Increment the count to determine whether or not to output
    outputcount++;
    
  }//end for t 
}
                         

void initialSchwarzchild(TwoDVectorGridFunction<complex<double>>& uh, Grid grd) {
  GridFunction<double> rho(uh.gridDim(), uh.pointsDim(), false);
  rho=grd.gridNodeLocations();
  ofstream fs;
  fs.open("initialdata.txt");
  for(int i = 0; i < uh.gridDim(); i++) {
    for (int j = 0; j < uh.pointsDim(); j++) {
      for (int n = 0; n < uh.modesDim(); n++) {
        double modeval = exp(-0.5 * pow((rho.get(i,j) / params.schw.sigma), 2.0));
        uh.set(n,0,i,j,0.0);
        if(!params.opts.useSource){
          uh.set(n,2,i,j,modeval);
        }
        
        uh.set(n,1,i,j,0.0);
        fs << rho.get(i,j) << " " << uh.get(0,2,i,j) << endl;

      }

      if(params.opts.useSource){
        double * r = &grd.rschw.get(i)[0];
        double * win = &grd.window.get(i)[0];
        double * dwin = &grd.dwindow.get(i)[0];
        double * d2win = &grd.d2window.get(i)[0];
        calc_window(params.grid.elemorder+1,r,win, dwin, d2win);
      } 
      double dxmin = fabs(grd.gridNodeLocations().get(0,0)
                          -grd.gridNodeLocations().get(0,2));
    }
    
  }
  fs.close();
}
  

void initialSinusoid(TwoDVectorGridFunction<complex<double>>& uh, Grid grd){
  double omega = 2.0 * PI / params.sine.wavelength;
  
  GridFunction<double> nodes(uh.gridDim(), uh.pointsDim(), false);
  nodes=grd.gridNodeLocations();

  for(int k = 0; k < uh.modesDim(); k++) {
    for(int i = 0; i < uh.gridDim(); i++){
      for (int j = 0; j < uh.pointsDim(); j++){
        double psi = params.sine.amp * sin(omega * nodes.get(i, j)
                                           + params.sine.phase);
        double pivar = omega * params.sine.amp * cos(omega * nodes.get(i, j)
                                                     + params.sine.phase);
        double rho = -params.waveeq.speed * pivar;
        //travelling wave
        //are rho and pi backward?
        uh.set(k, 0, i, j, psi);
        uh.set(k, 1, i, j, rho);
        uh.set(k, 2, i, j, pivar);
      }
    }
  }
}
void initialGaussian(TwoDVectorGridFunction<complex<double>>& uh, Grid grd){
  GridFunction<double> nodes(uh.gridDim(), uh.pointsDim(), false);
  nodes=grd.gridNodeLocations();
  
  for(int k = 0; k < uh.modesDim(); k++) {
    for(int i = 0; i < uh.gridDim(); i++){
      for(int j = 0; j < uh.pointsDim(); j++){
        double gaussian = params.gauss.amp * exp(-pow((nodes.get(i, j)
                                                       - params.gauss.mu), 2.0)
                                                 / 2.0 
                                                 / pow(params.gauss.sigma, 2.0));
        double dgauss = -(nodes.get(i, j) - params.gauss.mu)
          / pow(params.gauss.sigma, 2.0) * gaussian;
        uh.set(k, 0, i, j, gaussian);
        uh.set(k, 1, i, j, 0.0); //time derivative is zero
        //starts at center and splits
        uh.set(k, 2, i, j, dgauss);
      }
    }
  }
}

  double LTwoerror(Grid thegrid, TwoDVectorGridFunction<complex<double>>& uh0, 
                   TwoDVectorGridFunction<complex<double>>& uhend)
{
  //FIX THIS SO IT DEALS WITH SUM OF MODES
  //The square root of the integral of the squared difference
  double L2;
  L2 = 0.0;
  Array1D<double> weights;
  weights = thegrid.refelem.getw();
  GridFunction<double> nodes(uh0.gridDim(), uh0.pointsDim(), false);
  nodes = thegrid.gridNodeLocations();
  for(int i = 0; i < uh0.gridDim(); i++){
    for(int j = 0; j < uh0.pointsDim(); j++){
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

