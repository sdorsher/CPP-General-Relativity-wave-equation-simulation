#ifndef WRITEFILE_H
#define WRITEFILE_H

#include "TwoDVectorGridFunction.h"
#include "Grid.h"
#include "DiffEq.h"
#include <string>
#include <complex>

void write_fixed_time(OutputIndices& ijoutput, int& k, double t, TwoDVectorGridFunction<complex<double>>& uh,
		      TwoDVectorGridFunction<complex<double>>& RHStdvgf,
		      Grid& thegrid, DiffEq& theequation, Modes& lmmodes, bool append, 
                      string filename,
		      int type);

void write_fixed_radius(OutputIndices& ijoutput, int& k, double t, TwoDVectorGridFunction<complex<double>>& uh,
			TwoDVectorGridFunction<complex<double>>& RHStdvgf,
                        Grid& thegrid, DiffEq& theequation, Modes& lmmodes, bool append, 
                        string filename,
                        int type);




#endif
