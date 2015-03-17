#include "Entropy.hpp"
#include <cmath>
using namespace std;
//calculating the Entropy of the prob array.
//since it is a prob array, we simply do -1*sum(p*logP)
//using base of _base for log
//this assuming a discrete distribution.
//will try to write the samples for entropy like matlab
double Entropy(Matrix<double>& _p, const double& _base)
{
  double base_value=_base;
  
  if(base_value<0 )
    {
      base_value=exp(1);
    }
  if(base_value<1)
    {
      cout<<"WARNING: a value less than 1 is specified as the log base, be careful"<<endl;
    }
  /*if(_p.dim()!=1)
    {
      cout<<"the input matrix is not one dimension vector, please check your input"<<endl;
      cerr<<"the input matrix is not one dimension vector, please check your input"<<endl;
      exit(-1);
      }*/

  //start doing the computing
  unsigned total=_p.nTotal();
  double ret=0;
  for(unsigned i=0;i<total;i++)
    {
      ret+=-1*_p(i)*log(_p.Get1DArrayElement(i))/log(base_value);
    }
  return ret;
}

void CalculateAssignmentEntropy(VDJ_cuts_insertions_dinuc_ntbias_model& _model)
{
  _model.S_V=Entropy(_model.PV, exp(1));
  _model.S_DJ=Entropy(_model.PDJ, exp(1));
  unsigned dim_size[]={0,0,0,0}
  dim_size[0]=-1; 
  _model.S_D=Entropy(sum(_model.PDJ, 1), exp(1));
  _model.S_J=Entropy(sum(_model.PDJ, 0), exp(1));

  _model.S_gene=_model.S_V+_model.S_DJ;
  
}
