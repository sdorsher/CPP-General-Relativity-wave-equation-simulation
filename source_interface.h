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
#include "orbit.h"
#include "Grid.h"

namespace source_interface {
  static std::vector<EffectiveSource*> effsource;

  void init_source ( const Modes& lmmodes, const double &M );

    void set_window ( const double& r1, const double& w1, const double& q1,
		      const double& s1, const double& r2, const double& w2,
		      const double& q2, const double& s2, const int& nmodes );

    void calc_window ( const int& n, const double r[],
		       double Win[], double dWin[], double d2Win[] );

    void set_time_window ( const double& T, const double& dT_dt,
		           const double& d2T_dt2, const int& nmodes );

    void set_particle ( const double& p, const double& e,
		        const double& chi, const double& phi,
                        const int& nmodes );

    void eval_source (const int& mode,  const double& r,
                      std::complex<double> &src);

    void eval_source_all (const int& mode,  const int& n, const double r[],
                          const double Win[], const double dWin[],
                          const double d2Win[], complex<double> src[]);

    void clean_source ();

    void Phi ( const int* mode, const double* r,
               double* phire, double* phiim );

    void dPhi_dr ( const int* mode, const double* r,
                   double* dphidrre, double* dphidrim ) ;

    void dPhi_dt ( const int* mode, const double* r,
                   double* dphidtre, double* dphidtim );

    void fill_source(Grid& thegrid, double& time, int& nummodes);

    void fill_source_all(Grid thegrid, double time, int nummodes);
}


#endif
