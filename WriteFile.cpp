#include "WriteFile.h"


//write a file at a fixed time step, as a function of computational coordinate
void write_fixed_time(OutputIndices& ijoutput, int& k,double t, TwoDVectorGridFunction<complex<double>>& uh,
		      TwoDVectorGridFunction<complex<double>>& RHStdvgf,
		      Grid& thegrid, DiffEq& theequation, Modes& lmmodes, bool append, 
                      string filename,
		      int type)


{
  ofstream fs;
  ostringstream oss;
  oss << filename << "." << k << ".txt";
  if(append){
    fs.open(oss.str(), ios::app);
  }else{
    fs.open(oss.str());
  }
  fs << endl << endl;
  fs << " #time = " << t << endl;

  switch (type){
  case 1: // uh

    for (int i = 0; i < uh.GFvecDim(); i++){
      for(int j = 0; j < uh.GFarrDim(); j++){
	//Print out at select time steps
	fs << thegrid.gridNodeLocations().get(i, j) << " "
	   << uh.get(k, 0, i, j).real() << " " 
	   << uh.get(k, 1, i, j).real() <<" " 
	   << uh.get(k, 2, i, j).real()<< endl;
      }
    }
    break;
  case 2: //source
    for (int i = 0; i < uh.GFvecDim(); i++){
      for(int j = 0; j < uh.GFarrDim(); j++){
  	//Print out at select time steps
	fs << thegrid.gridNodeLocations().get(i, j) << " "
	    << theequation.source.get(k, i, j).real() << " " 
	    << theequation.source.get(k, i, j).imag() << endl; 
      }
    }
    break;
  case 3: //rh
   for (int i = 0; i < uh.GFvecDim(); i++){
      for(int j = 0; j < uh.GFarrDim(); j++){
	fs << thegrid.gridNodeLocations().get(i,j) << " "
	    << RHStdvgf.get(k,0,i,j).real() << " "
	    << RHStdvgf.get(k,1,i,j).real() << " "
	    << RHStdvgf.get(k,2,i,j).real() << " "
	    << RHStdvgf.get(k,0,i,j).imag() << " "
	    << RHStdvgf.get(k,1,i,j).imag() << " "
	    << RHStdvgf.get(k,2,i,j).imag() << endl;
      }
   }
   break;
  default:
    throw invalid_argument("Ivalid type in write_fixed_time");
    break;
  }
  fs.close();
}



//write a file at a fixed radius (known output indices), as a function of time 
void write_fixed_radius(OutputIndices& ijoutput, int& k, double t, TwoDVectorGridFunction<complex<double>>& uh,
                        TwoDVectorGridFunction<complex<double>>& RHStdvgf,
                        Grid& thegrid, DiffEq& theequation, Modes& lmmodes, bool append, 
                        string filename,
                        int type)
{
  ofstream fs;
  ostringstream oss;
  oss << filename << "." << k << ".txt";
  if(append){
    fs.open(oss.str(),ios::app);
  }else{
    fs.open(oss.str());
  }
  switch(type)
    {
    case 1: //uh
      fs << thegrid.gridNodeLocations().get(ijoutput.ifinite, ijoutput.jfinite) << " " 
	 << t << " "
	 << uh.get(k, 0, ijoutput.ifinite, ijoutput.jfinite).real() << " " 
	 << uh.get(k, 1, ijoutput.ifinite, ijoutput.jfinite).real() <<" " 
	 << uh.get(k, 2, ijoutput.ifinite, ijoutput.jfinite).real()<< " " 
	 << thegrid.gridNodeLocations().get(ijoutput.iSplus, ijoutput.jSplus) << " " 
	 << uh.get(k, 0, ijoutput.iSplus, ijoutput.jSplus).real() << " " 
	 << uh.get(k, 1, ijoutput.iSplus, ijoutput.jSplus).real() <<" " 
	 << uh.get(k, 2, ijoutput.iSplus, ijoutput.jSplus).real()<< endl;
      break;
    default:
      throw invalid_argument("Type out of range in write_fixed_radius");
      break;
    }//end switch case
  fs.close();
}

//write a file for psi summed over modes, possibly with some derivative. 
void write_summed_psi(OutputIndices& ijoutput, int& k, double t, TwoDVectorGridFunction<complex<double>>& uh,
                        TwoDVectorGridFunction<complex<double>>& RHStdvgf,
                        Grid& thegrid, DiffEq& theequation, Modes& lmmodes, bool append, 
                        string filename,
                        int type)
{
  ofstream fs;
  ostringstream oss;
  oss << filename << ".txt";
  if(append){
    fs.open(oss.str(),ios::app);
  }else{
    fs.open(oss.str());
  }
  switch(type)
    {
    case 1: //psi
      //	  if(k==params.modes.lmax){
      //	    oss7 << "psil.txt";
      fs << t << " " << chi <<  " " << phi << " " << p  << " " << e << " ";
      for(int n = 0; n<lmmodes.psil.size(); n++){
	fs << lmmodes.psil.at(n) << " ";
      }
      fs << endl;
      break;
    case 2: //dpsi/dt
      //oss8 << "psitl.txt";
      fs << t << " " << chi <<  " " << phi << " " << p  << " " << e << " ";
      for(int n = 0; n<lmmodes.psitl.size(); n++){
	fs << lmmodes.psitl.at(n) << " ";
      }
      fs << endl;
      break;
    case 3: //dpsi/dphi
      //      oss9 << "psiphil.txt";
      fs << t << " " << chi <<  " " << phi << " " << p  << " " << e << " ";
      for(int n = 0; n<lmmodes.psiphil.size(); n++){
	fs << lmmodes.psiphil.at(n) << " ";
      }
      fs << endl;
      break;
    case 4: //dpsi/dr
      //oss10 << "psirl.txt";
      fs << t << " " << chi <<  " " << phi << " " << p  << " " << e << " ";
      for(int n = 0; n<lmmodes.psirl.size(); n++){
	fs << lmmodes.psirl.at(n) << " ";
      }
      fs << endl;
      break;
    default:
      throw invalid_argument("Type out of range in write_fixed_radius");
      break;
    }//end switch case
  fs.close();
}
