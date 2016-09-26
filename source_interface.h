#ifndef SOURCE_INTERFACE_H
#define SOURCE_INTERFACE_H

#include <complex>
#include <cassert>
#include <iostream>
#include <vector>
#include "EffectiveSource.h"
#include "Modes.h"
#include "numerics.h"
#include "namespaces.h"
#include "Orbit.h"
#include "Grid.h"
#include "GridFunction.h"
#include "VectorGridFunction.h"
#include <omp.h>
#include "EllipticalOrbit.h"
#include "CircularOrbit.h"
#include "Coordinates.h"

//might require modifying to use Orbit.h instead of orbit.h

/*Interfaces with Barry Wardell's effective source package. Mostly
  copied directly (copy and paste) from source_interface.cc in the
  Fortran code. Modified to not expect all arguments to be passed by
  reference.*/

namespace source_interface {

  using namespace std;
  using namespace window;

  //vector for all lmmodes that contains EffectiveSource objects from Barry Wardell's
  //effective source package.
  extern vector<EffectiveSource*> effsource;


  //constructor for vector
  void init_source ( const Modes& lmmodes, const double &M );


  //Sets the window function over all elements of the vector
    void set_window ( const double& r1, const double& w1, const double& q1,
		      const double& s1, const double& r2, const double& w2,
		      const double& q2, const double& s2, const int& nmodes );

    void set_window_params(Coordinates & coords);
    
    //Calculates the window that is applied to the effective source
    void calc_window ( const int& n, const double r[],
		       double Win[], double dWin[], double d2Win[] );


    //Sets the window in time used to turn on the effective source slowly
    void set_time_window ( const double& T, const double& dT_dt,
		           const double& d2T_dt2, const int& nmodes );

    //Sets the particle position
    void set_particle ( const double& p, const double& e,
		        const double& chi, const double& phi,
                        const int& nmodes );

    //Evaluates the effective source for a mode, and a position in Schwarzschild coordinates
    void eval_source (const int& mode,  const double& r,
                      std::complex<double> &src);


    //Evaluates the effective source for a mode, an order of the nodes, a set of node
    //coordinates in Schwarzschild coordinates,
    // and a window function and its first and second deerivative
    void eval_source_all (const int& mode,  const int& n, const double r[],
                          const double Win[], const double dWin[],
                          const double d2Win[], complex<double> src[]);

    void clean_source ();

    //Find the phi coordinate
    void Phi ( const int* mode, const double* r,
               double* phire, double* phiim );


    //Find dPhi/dr
    void dPhi_dr ( const int* mode, const double* r,
                   double* dphidrre, double* dphidrim ) ;

    //Find dPhi/dt
    void dPhi_dt ( const int* mode, const double* r,
                   double* dphidtre, double* dphidtim );


    //Same thing as fill_source_all, but for ony one node.
    /*void fill_source(Grid& thegrid, double& time, int& nummodes,
		     VectorGridFunction<complex<double>>& source,
		     GridFunction<double>& window,
		     GridFunction<double>& dwindow,
		     GridFunction<double>& d2window);
    */
    //An overarching routine that handles setting particle position, setting up the time
    //window, evaluating the source, and
    // setting the source in the DiffEq object for all modes and nodes within an element.
    void fill_source_all(Grid& thegrid, double time, int nummodes,
			 VectorGridFunction<complex<double>>& source,
			 GridFunction<double>& window,
			 GridFunction<double>& dwindow,
			 GridFunction<double>& d2window, Orbit * orb);
}


#endif
