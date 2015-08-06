#include "source_interface.h"

namespace source_interface
{

  void init_source ( const Modes& lmmodes, const double &M ) {
      assert(effsource.empty());
      for (int i=0; i<lmmodes.ntotal; i++) {
        
        effsource.push_back(new EffectiveSource(lmmodes.ll[i], 
                                                lmmodes.mm[i], M));
      }
    }
    void set_window ( const double& r1, const double& w1, const double& q1,
		      const double& s1, const double& r2, const double& w2,
		      const double& q2, const double& s2, const int& nmodes ) {
      assert(effsource.size()==nmodes);
      for (auto& x : effsource) 
	(*x).set_window ( r1, w1, q1, s1, r2, w2, q2, s2 );
    }
    void calc_window ( const int& n, const double r[],
		       double Win[], double dWin[], double d2Win[] ) {
      assert(!effsource.empty());
      (*effsource.at(0)).calc_window ( n, r, Win, dWin, d2Win );
    }

    void set_time_window ( const double& T, const double& dT_dt,
		           const double& d2T_dt2, const int& nmodes ) {
      assert(effsource.size()==nmodes);
      for (auto& x : effsource)
        (*x).set_time_window ( T, dT_dt, d2T_dt2 );
    }

    void set_particle ( const double& p, const double& e,
		        const double& chi, const double& phi,
                        const int& nmodes ) {
      assert(effsource.size()==nmodes);
      for (auto& x : effsource) 
	(*x).set_particle ( p, e, chi, phi );
    }

    void eval_source (const int& mode,  const double& r,
                      std::complex<double> &src) {
      assert(effsource.size()>=mode);
      src =(*(effsource.at(mode)))(r);
    }
    void eval_source_all (const int& mode,  const int& n, const double r[],
                          const double Win[], const double dWin[],
                          const double d2Win[], complex<double> src[]) {
      assert(effsource.size()>=mode);
      (*(effsource.at(mode)))(n, r, Win, dWin, d2Win, src);
    }


    void clean_source () {
      int s{effsource.size()};
      for (int i=0; i<s; i++) {
        assert(effsource.at(s-1-i)!=nullptr);
        delete effsource.at(s-1-i);
        effsource.pop_back(); 
      }
      assert(effsource.size()==0);
    }
    void Phi ( const int* mode, const double* r,
               double* phire, double* phiim ) {
      assert(effsource.size()>=*mode);
      std::complex<double> phi{(*(effsource.at(*mode))).Phi(*r)};
      *phire = std::real(phi);
      *phiim = std::imag(phi);
    }

    void dPhi_dr ( const int* mode, const double* r,
                   double* dphidrre, double* dphidrim ) {
      assert(effsource.size()>=*mode);
      std::complex<double> dphidr{(*(effsource.at(*mode))).dPhi_dr(*r)};
      *dphidrre = std::real(dphidr);
      *dphidrim = std::imag(dphidr);
    }
    void dPhi_dt ( const int* mode, const double* r,
                   double* dphidtre, double* dphidtim ) {
      assert(effsource.size()>=*mode);
      std::complex<double> dphidt{(*(effsource.at(*mode))).dPhi_dt(*r)};
      *dphidtre = std::real(dphidt);
      *dphidtim = std::imag(dphidt);
    }
    void fill_source(Grid& thegrid, double& time, int& nummodes)
    {

      using namespace orbit;
      complex<double> src; 

      phi = phi_of_t(time);
      set_particle(p,e,chi,phi,nummodes);
      //set particle postion

      for(int i = 0; i<thegrid.numberElements(); i++){
        for(int k=0; k<nummodes; k++) {
          for(int j=0; j<thegrid.nodeOrder()+1; j++) {
            eval_source(k,thegrid.rschw.get(i,j),src);
            if((i==thegrid.numberElements()-1)&&(j==thegrid.nodeOrder())) {
              src = {0.0,0.0};
            }
            thegrid.source.set(i,j,k,src);
          }
        }
      }
    }

    void fill_source_all(Grid thegrid, double time, int nummodes){

      using namespace orbit;

      double tfac, dtfac_dt, d2tfac_dt2;


      phi= phi_of_t(time);
      
      set_particle(p,e,chi,phi,nummodes);
      
      if(params.opts.turn_on_source_smoothly){
        time_window(time, params.timewindow.tsigma, params.timewindow.torder, 
                    tfac, dtfac_dt, d2tfac_dt2);
      }else{
        tfac = 1.0; dtfac_dt = 0.0; d2tfac_dt2 = 0.0;
      }
      //      cout << "tfac = " << tfac << endl;
      // cout << "dtfac_dt = " << dtfac_dt << endl;
      //cout << "d2tfac_dt2 = " << d2tfac_dt2 << endl;
      
      set_time_window(tfac,dtfac_dt, d2tfac_dt2, nummodes);
      
      for(int i=0; i<thegrid.numberElements(); i++) {
        double *r = &thegrid.rschw.get(i)[0];
        double * win = &thegrid.window.get(i)[0];
        double * dwin = &thegrid.dwindow.get(i)[0];
        double * d2win = &thegrid.d2window.get(i)[0];
        for(int k=0; k<nummodes; k++) {
          complex<double> * src = &thegrid.source.get(k,i)[0];
          eval_source_all(k,thegrid.nodeOrder()+1, r, win, dwin, d2win, src);
        }
      }
      for (int k=0; k<nummodes; k++) {
        thegrid.source.set(k,thegrid.numberElements()-1,thegrid.nodeOrder(), 
                           {0.0,0.0});
      }

    }

}
