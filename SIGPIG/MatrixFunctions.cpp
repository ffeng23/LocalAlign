#include "MatrixFunctions.hpp"
#include <iostream>

using namespace std;

double max_mf(const vector<double>& _m)
{
  unsigned size=_m.size();
  double temp=0;
  if(size>0)
    {
      temp=_m.at(0);
      for(unsigned int i=1;i<size;i++)
	{
	  if(_m.at(i)>temp)
	    {
	      temp=_m.at(i);
	    }
	}
    }
  else
    {
      cout<<"WARNING: zero lengthed array , no maximum returned"<<endl;
    }
  return temp;
}
double max_mf(const double*  _m, const unsigned& _len)
{
  
  double temp=0;
  if(_len>0)
    {
      temp=_m[0];
      for(unsigned int i=1;i<_len;i++)
	{
	  if(_m[i]>temp)
	    {
	      temp=_m[i];
	    }
	}
    }
  else
    {
      cout<<"WARNING: zero lengthed array , no maximum returned"<<endl;
    }
  return temp;
}

unsigned max_mf(const vector<unsigned>& _m)
{
  unsigned size=_m.size();
  unsigned temp=0;
  if(size>0)
    {
      temp=_m.at(0);
      for(unsigned int i=1;i<size;i++)
	{
	  if(_m.at(i)>temp)
	    {
	      temp=_m.at(i);
	    }
	}
    }
  else
    {
      cout<<"WARNING: zero lengthed array , no maximum returned"<<endl;
    }
  return temp;
}
unsigned max_mf(const unsigned* _m, const unsigned& _len)

{
  
  unsigned temp=0;
  if(_len>0)
    {
      temp=_m[0];
      for(unsigned int i=1;i<_len;i++)
	{
	  if(_m[i]>temp)
	    {
	      temp=_m[i];
	    }
	}
    }
  else
    {
      cout<<"WARNING: zero lengthed array , no maximum returned"<<endl;
    }
  return temp;
}




double min_mf(const vector<double>& _m)
{
  unsigned size=_m.size();
  double temp=0;
  if(size>0)
    {
      temp=_m.at(0);
      for(unsigned int i=1;i<size;i++)
	{
	  if(_m.at(i)<temp)
	    {
	      temp=_m.at(i);
	    }
	}
    }
  else
    {
      cout<<"WARNING: zero lengthed array , no maximum returned"<<endl;
    }
  return temp;
}
double min_mf(const double*  _m, const unsigned& _len)
{
  
  double temp=0;
  if(_len>0)
    {
      temp=_m[0];
      for(unsigned int i=1;i<_len;i++)
	{
	  if(_m[i]<temp)
	    {
	      temp=_m[i];
	    }
	}
    }
  else
    {
      cout<<"WARNING: zero lengthed array , no maximum returned"<<endl;
    }
  return temp;
}
unsigned min_mf(const vector<unsigned>& _m)
{
  unsigned size=_m.size();
  unsigned temp=0;
  if(size>0)
    {
      temp=_m.at(0);
      for(unsigned int i=1;i<size;i++)
	{
	  if(_m.at(i)<temp)
	    {
	      temp=_m.at(i);
	    }
	}
    }
  else
    {
      cout<<"WARNING: zero lengthed array , no maximum returned"<<endl;
    }
  return temp;
}
unsigned min_mf(const unsigned* _m, const unsigned& _len)
{
  
  unsigned temp=0;
  if(_len>0)
    {
      temp=_m[0];
      for(unsigned int i=1;i<_len;i++)
	{
	  if(_m[i]<temp)
	    {
	      temp=_m[i];
	    }
	}
    }
  else
    {
      cout<<"WARNING: zero lengthed array , no maximum returned"<<endl;
    }
  return temp;
}
