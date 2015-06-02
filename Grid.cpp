#include "Grid.h"

Grid::Grid(int elemorder, int numelements,double lowerlim, double upperlim):
  order{elemorder},
  NumElem{numelements},
  nodeLocs{0,elemorder+1}, 
  Amatrices{numelements,elemorder+1},
  Bmatrices{numelements,elemorder+1},
  refelem{elemorder}

{
  for(int i=0; i<=numelements; i++){
    elementBoundaries.push_back(lowerlim + i*(upperlim-lowerlim)/float(numelements));
  }
  
  Array1D<double> physicalPosition(elemorder+1);
  for(int elem=0; elem<numelements; elem++){
    physicalPosition = ((elementBoundaries[elem+1] - elementBoundaries[elem]) / 2.0)
      *refelem.getr()
      +((elementBoundaries[elem+1] + elementBoundaries[elem]) / 2.0);
    nodeLocs.append(physicalPosition);
  }

  calcjacobian();
  Amatrices=setupAmatrix(nodeLocs);
  Bmatrices=setupBmatrix(nodeLocs);

  for (int i=0; i<nodeLocs.gridDim(); i++){
    CharacteristicFlux left(Amatrices.get(i,0));
    CharacteristicFlux right(Amatrices.get(i,nodeLocs.pointsDim()-1));
    AleftBoundaries.push_back(left);
    ArightBoundaries.push_back(right);
  }
  

  duL.resize(NumElem);
  duR.resize(NumElem);

  //need some kind of A matrix as a function of position or function
  //build A matrix Gridfunction
  //build CharacteristicFlux vector-- make that a function of position (how)
  

}

vector<TNT::Array2D<double>> Grid::characteristicflux(VectorGridFunction<double>& uh)
{
  vector<Array2D<double>> du;
  du.resize(NumElem);
  for(int elemnum=0; elemnum<NumElem; elemnum++){
     int indL=0; //index of leftmost node
    int indR=uh.pointsDim()-1; //index of rightmost node
    double nL=-1.0; //normal
    double nR=1.0;
    //vmin and vmax are min and max indices in vector dimension (psi, rho phi)
    int vmaxL=AleftBoundaries[elemnum].getAdim()-1;
    int DdimL = AleftBoundaries[elemnum].getDdim();
    int vminL=vmaxL-DdimL+1;
    int vmaxR=ArightBoundaries[elemnum].getAdim()-1;
    int DdimR=ArightBoundaries[elemnum].getDdim();
    int vminR=vmaxR-DdimR+1;

    Array1D<double> uintL(DdimL);
    Array1D<double> uintR(DdimR);
    Array1D<double> uextL(DdimL);
    Array1D<double> uextR(DdimR);
    
    uintL=uh.getVectorAsArray1D(elemnum,indL,vminL,vmaxL); //internal u
    uintR=uh.getVectorAsArray1D(elemnum,indR,vminR,vmaxR);

    
    if(elemnum>0){
      uextL=uh.getVectorAsArray1D(elemnum-1,indR,vminL,vmaxL); //external u, left boundary
    }else{
      uextL=uh.getVectorAsArray1D(NumElem-1,indR,vminL,vmaxL); //periodic boundary conditions
    }

    if(elemnum<NumElem-1){
      uextR=uh.getVectorAsArray1D(elemnum+1,indL,vminR,vmaxR); //external u, right boundary
    }else{
      uextR=uh.getVectorAsArray1D(0,indL,vminR,vmaxR); //periodic boundary conditions
    }
    

    if(elemnum==5)
      {
        //cout << "LHS0 uext=" << uextL[0] << " uint=" << uintL[0] << endl;
        //cout << "RHS0 uext=" << uextR[0] << " uint=" << uintR[0] <<endl;
        //cout << "LHS1 uext=" << uextL[1] << " uint=" << uintL[1] << endl;
        // cout << "RHS1 uext=" << uextR[1] << " uint=" << uintR[1] <<endl;
        //correct at first time step, subsequently too small
      }
    Array2D<double> lambdaminusL(DdimL,DdimL,0.0);
    Array2D<double> lambdaminusR(DdimR,DdimR,0.0);
    Array2D<double> lambdaplusL(DdimL,DdimL,0.0);
    Array2D<double> lambdaplusR(DdimR,DdimR,0.0);

    Array2D<double> lambdaL= AleftBoundaries[elemnum].getLambda();
    Array2D<double> lambdaR= ArightBoundaries[elemnum].getLambda();


    for(int j=0; j<DdimL; j++){
      if(nL*lambdaL[j][j]<=0){
        lambdaminusL[j][j]=nL*lambdaL[j][j];
      }else{
        lambdaplusL[j][j]=nL*lambdaL[j][j];
      }
      
      if(nR*lambdaR[j][j]<=0){
        lambdaminusR[j][j]=nR*lambdaR[j][j];
      }else{
        lambdaplusR[j][j]=nR*lambdaR[j][j];
      }
    }

    if(elemnum==5)
      {
        //        output2D(lambdaplusR);
        //  output2D(lambdaminusR);
      }

    Array2D<double> sinvL=AleftBoundaries[elemnum].getSinv();
    Array2D<double> sinvR=ArightBoundaries[elemnum].getSinv();
    Array2D<double> SL=AleftBoundaries[elemnum].getS();
    Array2D<double> SR=ArightBoundaries[elemnum].getS();
    

    if(elemnum==5)
      {
        //        output2D(sinvL);
        // output2D(sinvR);
        //output2D(SL);
        // output2D(SR);
      }
    
    Array1D<double> nfluxL=matmult(lambdaplusL,matmult(sinvL,uintL));
    nfluxL+= matmult(lambdaminusL,matmult(sinvL,uextL));
    nfluxL=matmult(SL,nfluxL);

    Array1D<double> nfluxR=matmult(lambdaplusR,matmult(sinvR,uintR));
    nfluxR+= matmult(lambdaminusR,matmult(sinvR,uextR));
    nfluxR=matmult(SR,nfluxR);


    if(elemnum==5)
      {
        //        output1D(nfluxL);
        // output1D(nfluxR);
      }

    Array2D<double> AtrimmedL= AleftBoundaries[elemnum].getAtrimmed();
    Array2D<double> AtrimmedR= AleftBoundaries[elemnum].getAtrimmed();
    Array2D<double> duelem(AtrimmedR.dim1(),2,0.0);


    
    Array1D<double> duL = nL*matmult(AtrimmedL,uintL)-nfluxL; //needs - sign?
    Array1D<double> duR = nR*matmult(AtrimmedR,uintR)-nfluxR; //needs - sign?

    if(elemnum==5){

      //      output1D(duL);
      //output1D(duR);
    }
    insert_1D_into_2D(duelem,duL,0,false);
    insert_1D_into_2D(duelem,duR,1,false);


    du[elemnum]=duelem;


  }
  return du;
}

void Grid::RHS(VectorGridFunction<double>& uh, 
               VectorGridFunction<double>& RHSvgf, double t, vector<Array2D<double>>& du )
{

  //  output2D(du[5]);

  //problem is in this routine HERE

  for(int elemnum=0; elemnum<NumElem; elemnum++){
    int vmaxAB=ArightBoundaries[elemnum].getAdim()-1;
    int vminA=vmaxAB-ArightBoundaries[elemnum].getDdim()+1;
    Array2D<double> RHSB(uh.pointsDim(),ArightBoundaries[elemnum].getAdim());
                         
    for(int nodenum=0; nodenum<uh.pointsDim(); nodenum++){
      Array1D<double> RHSBpernode;
      RHSBpernode=matmult(Bmatrices.get(elemnum,nodenum),
                          uh.getVectorAsArray1D(elemnum,nodenum,0,vmaxAB));
      insert_1D_into_2D(RHSB,RHSBpernode,nodenum,false);

    }//this can be sped up by skipping the insert step and reading directly from per node


    //A contribution:
    
    if(elemnum==5)
      {
        //        output2D(Bmatrices.get(elemnum,0));
        //     output2D(RHSB);
        // output1D(uh.getVectorAsArray1D(elemnum,0,0,vmaxAB));
       //correct based on inputs but for some reason uh[1] is still zero for second timestep
       // output1D(uh.get(1,elemnum));
      }

        Array2D<double> RHSA1=jacobian(elemnum)*(matmult(refelem.getD(),
                                  uh.getVectorNodeArray2D(elemnum,
                                                          vminA,vmaxAB)));
    /*
    Array2D<double> prederiv(uh.pointsDim(),ArightBoundaries[elemnum].getDdim());
    for(int nodenum=0; node<uh.pointsDim(); nodenum++)
      {

        Array1D<double> prederivpernode = matmult(trimmedAmatrices.get(elemnum,nodenum),uh.getVectorAsArray1D(elemnum,nodenum,vminA,vmaxAB));
        
        insert_1D_into_2D(prederiv,prederivpernode,nodenum,false);
      }        

    Array2D<double> RHSA1=jacobian(elemnum)*(matmult(refelem.getD(),prederiv);
    */
    //needs a multiplication by an A matrix between D and vectornodearray,
    //but A is position dependent. not sure how to handle this.

    if(elemnum==5)
      {
        // output1D(uh.get(2,elemnum));
        //output2D(RHSA1);//setting 2 column, not 1 column. wrong
        output2D(uh.getVectorNodeArray2D(elemnum,vminA,vmaxAB));
      }

    Array2D<double> RHSA2=jacobian(elemnum)*matmult(refelem.getLift(),du[elemnum]);

       //RHSA and RHSB will have different sizes due to the different 
    //number of diffeq variables stored in each. sum them using a 
    //for loop while assigning values to the RHSvgf vector grid function

    Array2D<double> RHSA=RHSA1+RHSA2;

    for(int vecnum=0; vecnum<RHSvgf.vectorDim(); vecnum++)
      {
        for(int nodenum=0; nodenum<RHSvgf.pointsDim(); nodenum++)
          {
            if(vecnum<vminA){
              RHSvgf.set(vecnum,elemnum,nodenum,RHSB[nodenum][vecnum]);
            }else{
              RHSvgf.set(vecnum,elemnum,nodenum,RHSB[nodenum][vecnum]
                         +RHSA[nodenum][vecnum-vminA]);
            }
          }
      }
  }

}


/*Grid::Grid(string fileElemBoundaries,int elemorder, int numelements):order{elemorder},NumElem{numelements},nodeLocs{0,elemorder+1,false}, refelem{elemorder}
{
  ifstream fs;
  fs.open(fileElemBoundaries);
  double data;
  fs >> data;
  while(!fs.eof())
    {
      elementBoundaries.push_back(data);
      fs >> data;
    }
  fs.close();
  if(elementBoundaries.size()-1 != numelements)
    {
      throw invalid_argument("Element boundaries too long or too short for number of elements given.");
    }

  //  ReferenceElement refelem(elemord);
  // when we generalize this, use map to store order, element pairs
  // so they do not need to be recalculated with each element of the same
  // order
  Array1D<double> physicalPosition(elemorder+1);
  for(int elem=0; elem<numelements; elem++)
    {
      physicalPosition = ((elementBoundaries[elem+1]-elementBoundaries[elem]) / 2.0)
                          *refelem.getr()
                       + ((elementBoundaries[elem+1]+elementBoundaries[elem]) / 2.0);
      nodeLocs.append(physicalPosition);
    }
  calcjacobian();
}
*/
int Grid::numberElements()
{
  return NumElem;
}

GridFunction<double> Grid::gridNodeLocations()
{
  return nodeLocs;
}

vector<double> Grid::gridBoundaries()
{
  return elementBoundaries;
}

void Grid::calcjacobian()
{
  for(int elem=0; elem<NumElem; elem++)
    {
      double rx= 2.0/(elementBoundaries[elem+1]-elementBoundaries[elem]);
      //cout << elem << " " <<rx << endl;
      drdx.push_back(rx);
    }
}

double Grid::jacobian(int elemnum)
{
  return drdx[elemnum];
}
