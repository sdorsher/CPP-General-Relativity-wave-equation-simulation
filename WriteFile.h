#ifndef WRITEFILE_H
#define WRITEFILE_H

#include "TwoDVectorGridFunction.h"
#include "Grid.h"
#include "DiffEq.h"
#include <string>


void write_fixed_time(int& k, TwoDVectorGridFunction<double>& uh,
		      TwoDVectorGridFunction<double>& RHStdvgf,
		      Grid& thegrid, DiffEq& theequation, bool& append, string& filename,
		      int type);

void write_fixed_radius(int& k, TwoDVectorGridFunction<double>& uh,
			TwoDVectorGridFunction<double>& RHStdvgf,
		      Grid& thegrid, DiffEq& theequation, bool& append, string& filename,
		      int type);




#endif
