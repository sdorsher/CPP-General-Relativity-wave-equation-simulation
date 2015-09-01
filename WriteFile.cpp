#include "WriteFile.h"

void write_fixed_time(int& k, TwoDVectorGridFunction<double>& uh,
		      TwoDVectorGridFunction<double>& RHStdvgf,
		      Grid& thegrid, DiffEq& theequation, bool& append, string& filename,
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
  case 1:

    for (int i = 0; i < uh.gridDim(); i++){
      for(int j = 0; j < uh.pointsDim(); j++){
	//Print out at select time steps
	fs << thegrid.gridNodeLocations().get(i, j) << " "
	   << uh.get(k, 0, i, j).real() << " " 
	   << uh.get(k, 1, i, j).real() <<" " 
	   << uh.get(k, 2, i, j).real()<< endl;
      }
    }
    break;
  case 2:
    for (int i = 0; i < uh.gridDim(); i++){
      for(int j = 0; j < uh.pointsDim(); j++){
  	//Print out at select time steps
	fs << thegrid.gridNodeLocations().get(i, j) << " "
	    << theequation.source.get(k, i, j).real() << " " 
	    << theequation.source.get(k, i, j).imag() << endl; 
      }
    }
    fs5.close();
    fs6.close();
    break;
  case 3:
   for (int i = 0; i < uh.gridDim(); i++){
      for(int j = 0; j < uh.pointsDim(); j++){
	fs << thegrid.gridNodeLocations().get(i,j) << " "
	    << RHStdvgf.get(k,0,i,j).real() << " "
	    << RHStdvgf.get(k,1,i,j).real() << " "
	    << RHStdvgf.get(k,2,i,j).real() << " "
	    << RHStdvgf.get(k,0,i,j).imag() << " "
	    << RHStdvgf.get(k,1,i,j).imag() << " "
	    << RHStdvgf.get(k,2,i,j).imag() << endl;
      }
   }
  default:
    throw invalid_argument("Ivalid type in write_fixed_time");
    break;
  }
  fs.close();
}


void write_fixed_radius(int& k, TwoDVectorGridFunction<double>& uh,
			TwoDVectorGridFunction<double>& RHStdvgf,
		      Grid& thegrid, DiffEq& theequation, bool& append, string& filename,
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
    case 1:
      //	  oss << params.file.fixedradiusfilename << "." << k << ".txt";
      fs << nodes.get(ifinite, jfinite) << " " 
	 << t << " "
	 << uh.get(k, 0, ifinite, jfinite).real() << " " 
	 << uh.get(k, 1, ifinite, jfinite).real() <<" " 
	 << uh.get(k, 2, ifinite, jfinite).real()<< " " 
	 << nodes.get(iSplus, jSplus) << " " 
	 << uh.get(k, 0, iSplus, jSplus).real() << " " 
	 << uh.get(k, 1, iSplus, jSplus).real() <<" " 
	 << uh.get(k, 2, iSplus, jSplus).real()<< endl;
      break;
    case 2:
      //	  if(k==params.modes.lmax){
      //	    oss7 << "psil.txt";
      fs << t << " " << chi <<  " " << phi << " " << p  << " " << e << " ";
      for(int n = 0; n<lmmodes.psil.size(); n++){
	fs << lmmodes.psil.at(n) << " ";
      }
      fs << endl;
      break;
    case 3:
      //oss8 << "psitl.txt";
      fs << t << " " << chi <<  " " << phi << " " << p  << " " << e << " ";
      for(int n = 0; n<lmmodes.psitl.size(); n++){
	fs << lmmodes.psitl.at(n) << " ";
      }
      fs << endl;
      break;
    case 4:
      //      oss9 << "psiphil.txt";
      fs << t << " " << chi <<  " " << phi << " " << p  << " " << e << " ";
      for(int n = 0; n<lmmodes.psiphil.size(); n++){
	fs << lmmodes.psiphil.at(n) << " ";
      }
      fs << endl;
      break;
    case 5:
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
