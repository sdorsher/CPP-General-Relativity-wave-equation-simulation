#ifndef WRITEFILE_H
#define WRITEFILE_H

#include "TwoDVectorGridFunction.h"
#include "Grid.h"
#include "DiffEq.h"
#include <string>
#include <complex>
#include <limits>
#include "Orbit.h"
#include "namespaces.h"
#include "EllipticalOrbit.h"
#include "CircularOrbit.h"

//might need updating to use Orbit.h rather than orbit.h


//Output data at a fixed time
void write_fixed_time(int k, double t,
		      TwoDVectorGridFunction<complex<double>>& uh,
		      TwoDVectorGridFunction<complex<double>>& RHStdvgf,
		      Grid& thegrid, DiffEq& theequation, Modes& lmmodes, bool append, 
                      string filename,
		      int type, Orbit * orb);

//Output data at a fixed radius
void write_fixed_radius(OutputIndices& ijoutput, int& k, double t,
			TwoDVectorGridFunction<complex<double>>& uh,
			TwoDVectorGridFunction<complex<double>>& RHStdvgf,
                        Grid& thegrid, DiffEq& theequation, Modes& lmmodes, bool append, 
                        string filename,
                        int type);

//Output psi, summed over the modes
void write_summed_psi(OutputIndices& ijoutput, int& k, double t,
		      TwoDVectorGridFunction<complex<double>>& uh,
		      TwoDVectorGridFunction<complex<double>>& RHStdvgf,
		      Grid& thegrid, DiffEq& theequation, Modes& lmmodes, bool append, 
		      string filename,
		      int type, Orbit * orb);

#endif
