#include "source_interface.h"

// see source_interface.h for comments on what functions do.

namespace source_interface
{
  vector<EffectiveSource*> effsource;
  
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
	x->set_window ( r1, w1, q1, s1, r2, w2, q2, s2 );
  }

  void set_window_params(Coordinates & coords,Grid& thegrid, Modes & lmmodes){
    double deltar = thegrid.gridNodeLocations().get(0,1)-thegrid.gridNodeLocations().get(0,0);
    R1= coords.invert_tortoise(Rminus, params.schw.mass) + 2.0*params.schw.mass;
    R2 = coords.invert_tortoise(Rplus, params.schw.mass) + 2.0* params.schw.mass;
    w1 = params.schw.p_orb-(coords.invert_tortoise(2.0*deltar, params.schw.mass)
			    +2.0*params.schw.mass)-R1;
    w2 = R2 - (params.schw.p_orb + coords.invert_tortoise(2.0*deltar, params.schw.mass\
						   )+2.0*params.schw.mass);
    nmodes = lmmodes.ntotal;


    cout << "R1 R2 w1 w2" << endl;
    cout << R1 << " " << R2 << " " << w1 <<  " " << w2 << endl << endl;
    if(params.opts.useSource) {
      set_window(R1, w1, 1.0, 1.5, R2, w2, 1.0, 1.5, lmmodes.ntotal);
    }
    

  }
  
  void calc_window ( const int& n, const double r[],
		       double Win[], double dWin[], double d2Win[] ) {
      assert(!effsource.empty());
      auto & temp = effsource.at(0);
      temp->calc_window ( n, r, Win, dWin, d2Win );
      //for(int i=0; i<n; i++){
      //cout <<r[i] << " " <<  Win[i] << " " << dWin[i] << " " << d2Win[i] << endl;
      //}
    }

    void set_time_window ( const double& T, const double& dT_dt,
		           const double& d2T_dt2, const int& nmodes ) {
      assert(effsource.size()==nmodes);
      for (auto& x : effsource)
        x->set_time_window ( T, dT_dt, d2T_dt2 );
    }

    void set_particle ( const double& p, const double& e,
		        const double& chi, const double& phi,
                        const int& nmodes ) {
      assert(effsource.size()==nmodes);
      for (auto& x : effsource) 
	x->set_particle ( p, e, chi, phi );
    }

    void eval_source (const int& mode,  const double& r,
                      std::complex<double> &src) {

      //assert(effsource.size()>=mode);
      auto& temp = effsource.at(mode);
      src =(*temp)(r);
    }
    void eval_source_all (const int& mode,  const int& n, const double r[],
                          const double Win[], const double dWin[],
                          const double d2Win[], complex<double> src[]) {
      assert(effsource.size()>=mode);
      auto & temp = effsource.at(mode);
            //for(int i=0; i<17; i++){
      //	cout << Win[i] << endl;
      //}

      (*temp)(n, r, Win, dWin, d2Win, src);
    }


    void clean_source () {
      vector<complex<double>>::size_type s{effsource.size()};
      for (vector<complex<double>>::size_type i=0; i<s; i++) {
        assert(effsource.at(s-1-i)!=nullptr);
        delete effsource.at(s-1-i);
        effsource.pop_back(); 
      }
      assert(effsource.size()==0);
    }
    void Phi ( const int* mode, const double* r,
               double* phire, double* phiim ) {
      assert(effsource.size()>=*mode);
      auto & temp = effsource.at(*mode);
      auto temp2 = (*temp).Phi(*r);
      std::complex<double> phi{temp2};
      *phire = std::real(phi);
      *phiim = std::imag(phi);
    }

    void dPhi_dr ( const int* mode, const double* r,
                   double* dphidrre, double* dphidrim ) {
      assert(effsource.size()>=*mode);
      auto & temp = effsource.at(*mode);
      auto temp2 = temp->dPhi_dr(*r);
      std::complex<double> dphidr{temp2};
      *dphidrre = std::real(dphidr);
      *dphidrim = std::imag(dphidr);
    }
    void dPhi_dt ( const int* mode, const double* r,
                   double* dphidtre, double* dphidtim ) {
      assert(effsource.size()>=*mode);
      auto & temp = effsource.at(*mode);
      auto temp2 = (*temp).dPhi_dt(*r);
      std::complex<double> dphidt{temp2};
      *dphidtre = std::real(dphidt);
      *dphidtim = std::imag(dphidt);
    }
  /*void fill_source(Grid& thegrid, double& time, int& nummodes,
		   VectorGridFunction<complex<double>>& source,
		   GridFunction<double>& window,
		   GridFunction<double>& dwindow,
		   GridFunction<double>& d2window)
    {

      //using namespace orbit
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
            source.set(i,j,k,src);
          }
        }
      }
    }
  */
  void fill_source_all(Grid& thegrid, double time, int nummodes,
		       VectorGridFunction<complex<double>>& source,
		       GridFunction<double>& window,
		       GridFunction<double>& dwindow,
		       GridFunction<double>& d2window, Orbit * orb, Modes & lmmodes)
  {


    //using namespace orbit;

      double tfac, dtfac_dt, d2tfac_dt2;
      if(orb->orbType()==circular){
	CircularOrbit* corb = dynamic_cast<CircularOrbit*>(orb);
	orb->phi= corb->phi_of_t(time);
      }
      //cout << time << " "<< p << " " << e<< " "<< chi << " " << phi << " " << nummodes << endl;
 
      set_particle(orb->p,orb->e,orb->chi,orb->phi,lmmodes.ntotal);

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

      
      //#pragma omp parallel for
      for(int i=0; i<thegrid.numberElements(); i++) {
	vector<double> rschwv = thegrid.rschw.get(i);
	vector<double> windowv = window.get(i);
	vector<double> dwindowv = dwindow.get(i);
	vector<double> d2windowv = d2window.get(i);
	
        double *r = &rschwv[0];
        double * win = &windowv[0];
        double * dwin = &dwindowv[0];
        double * d2win = &d2windowv[0];
	//for(int ii=0; ii<17; ii++){
	  //cout <<r[ii] << " " <<  win[ii] << " " << dwin[ii] << " " << d2win[ii] << endl;
	//}
	for(int k=0; k<nummodes; k++) {
	  vector<complex<double>> temp = source.get(k,i);
          complex<double> * src = &temp[0];
	  eval_source_all(k,thegrid.nodeOrder()+1, r, win, dwin, d2win, src);
	  vector<complex<double>> src2(src, src+params.grid.elemorder+1);
	  
	  for(int j=0; j<= thegrid.nodeOrder(); j++){
	      double eps11 = 1.e-200;
	      //if(abs(win[j]>eps11)){
	      //cout << win[j] << endl;
	      //}
	      if((i==thegrid.numberElements()-1)&&(j==thegrid.nodeOrder())) {
		source.set(k,i,j,{0.,0.});
	      }else if((i==0)&&(j==0)){
		source.set(k,i,j,{0.,0.});
 	      }else{
		source.set(k,i,j,src2[j]);
	      }
	  }
	}
	
      }
      
      //#pragma omp parallel for if(nummodes>omp_get_max_threads())
      //      for (int k=0; k<nummodes; k++) {
      //source.set(k,thegrid.numberElements()-1,thegrid.nodeOrder(), 
      //	   {0.,0.});
      //}
      
  }
  
  
}//end namespace
