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

/**
 * Quicksort.
 * @param a - The array to be sorted.
 * @param first - The start of the sequence to be sorted.
 * @param last - The end of the sequence to be sorted.
 */
void QuickSort(double* _a, const unsigned& _first, const unsigned int& _last, unsigned* _index ) 
{
  int pivotElement;
 
  if(_first < _last)
    {
      pivotElement = Pivot(_a, _first, _last, _index);
      QuickSort(_a, _first, pivotElement-1, _index);
      QuickSort(_a, pivotElement+1, _last, _index);
    }
}
 
/**
 * Find and return the index of pivot element.
 * @param a - The array.
 * @param first - The start of the sequence.
 * @param last - The end of the sequence.
 * @return - the pivot element
 */
//here we use a new way to choose pivot. a median between [first], [middle] and [last]
unsigned Pivot(double* _a, const unsigned int& _first, const unsigned int& _last, unsigned* _index) 
{
  
  unsigned  p = _first;
  unsigned middle=(_first+_last)/2;
  double pivotElement;
  unsigned pivotIndex=_first;
  
  if(middle==_first||middle==_last)
    {
      pivotIndex = GetMedianIndex(_a, _first, _last, middle) ;
      //swap the index
      if(pivotIndex!=_first)
	{
	  swap(_a[_first], _a[pivotIndex]);
	}
    }
  //now we have put the pivot element in the [_first]
  pivotElement=_a[_first];
  cout<<"8888888888888>>pivote value:"<<_a[_first]<<endl;
  //if(_index!=NULL)
  //  _index[_first]=pivotIndex;
  Print(_a, _last-_first+1);
  for(unsigned int i = _first+1 ; i <= _last ; i++)
    {
      cout<<"\tround i:"<<i<<"--";
      
      /* If you want to sort the list in the other order, change "<=" to ">" */
      if(_a[i] <= pivotElement)
        {
	  p++;
	  cout<<"p:"<<p<<";i:"<<i<<endl;
	  swap(_a[i], _a[p]);
	  cout<<"********************SWAP*************!!!"<<endl;
        }
Print(_a, _last-_first+1);
    }
  cout<<"\t>>before pivot"<<endl;
  Print(_a, _last-_first+1);
  swap(_a[p], _a[_first]);
  cout<<"\t>>after"<<endl;
  Print(_a,_last-_first+1);
  
  return p;
}
 
//return 0,1,2 indicating which one is the median
//input:
// m is the array,
// a, b, c is the index of the element that need to be considered
// return index of the median value
unsigned GetMedianIndex(const double* m, const unsigned& a, const unsigned& b, const unsigned& c)
{
  //unsigned index=0;
  if(m[a]>m[b])
    {
      if(m[c]>m[a]) //a is the median
	{
	  return a;
	}
      else //c<=a//still not sure, compare b and c
	{
	  if(m[b]>m[c])
	    return b;
	  else
	    return c;
	}
    }
  else //a<=b
    {
      if(m[c]<m[a]) 
	return a;
      else //c>a need to compare b and c
	{
	  if(m[c]<m[b])
	    return c;
	  else //c>b
	    return b;
	}
    }
  return a;
}
 
/**
 * Swap the parameters.
 * @param a - The first parameter.
 * @param b - The second parameter.
 */
void Swap(int& _a, int& _b)
{
  int temp = _a;
  _a = _b;
  _b = temp;
}
 
/**
 * Print an array.
 * @param a - The array.
 * @param N - The size of the array.
 */
void Print(const double* a, const int& N)
{
  for(int i = 0 ; i < N ; i++)
    cout << "array[" << i << "] = " << a[i] << endl;
} 
