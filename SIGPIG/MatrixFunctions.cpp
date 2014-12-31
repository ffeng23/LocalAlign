#include "MatrixFunctions.hpp"
using std;
double max(vector<double> _m)
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
double max(const double*  _m, const unsigned& _len)
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

unsigned max(vector<unsigned> _m)
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
unsigned max(const unsigned* _m, const unsigned& _len)

{
  
  unsigned temp=0;
  if(_len>0)
    {
      temp=_m[0];
      for(unsigned int i=1;i<_len;i++)
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




double min(vector<double> _m)
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
double min(const double*  _m, const unsigned& _len)
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
unsigned min(vector<unsigned> _m)
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
unsigned min(const unsigned* _m, const unsigned& _len)
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
