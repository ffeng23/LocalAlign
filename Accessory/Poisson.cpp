#include "Poisson.hpp"
#include <iostream>
#include <math.h>

Poisson::Poisson(const double& _lamda)
{
  c_lamda=_lamda;
}

Poisson::~Poisson()
{
  //empty for now
}

double Poisson::Density(const unsigned long& _k, const double& _lamda)
{
  //need to check for input and decide whether we need to set lamda
  if(_lamda>0)
    {
      c_lamda=_lamda;
    }
  if(c_lamda<=0)
    {
      std::cerr<<"ERROR: the Poisson distribution has NOT been set up correctly. Parameter is necessary for this"<<std::endl;
      std::cout<<"ERROR: the Poisson distribution has NOT been set up correctly. Parameter is necessary for this"<<std::endl;
      exit(-1);
    }
  if(_k==0)
    {
      std::cerr<<"ERROR: request k=0 for the Poisson distribution"<<std::endl;
      std::cout<<"ERROR: request k=0 for the Poisson distribution"<<std::endl;
      exit(-1);
    }

  //now we are good
  double p=_k*log(c_lamda)+(-1*c_lamda);
  for(unsigned i=_k;i>=1;i--)
    {
      p-=log(i);
    }

  p=exp(p);
  return p;
}

double Poisson::Probability(const unsigned long& _k, const double& _lamda)
{
  double p=0;
  std::cerr<<"ERROR: in Poisson, function not implemented yet"<<std::endl;
  exit(-1);
  return p;
}

unsigned Quantile (const double& _p, const double& _lamda)
{
  unsigned p=0;
  std::cerr<<"ERROR: in Poisson, function not implemented yet"<<std::endl;
  std::cout<<"ERROR: in Poisson, function not implemented yet"<<std::endl;
  exit(-1);
  return p;
}

unsigned Rand(const double& _lamda)
{
  unsigned p=0;
  std::cerr<<"ERROR: in Poisson, function not implemented yet"<<std::endl;
  std::cout<<"ERROR: in Poisson, function not implemented yet"<<std::endl;
  exit(-1);
  return p;
}
