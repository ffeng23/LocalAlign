#include "MatrixFunctions.hpp"
#include <iostream>
#include <limits>
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


//this is the function to find max value among the 2 D array
//input: T** _m1: 2D array pointer
//       _len1: is the size of the first diment
//      _len2: is the size of the second dimenstion.
//           by doing this, we allow rectangle arrays
template<class T>
T max_mf2(const T* const* _m1, const unsigned& _len1, 
		const unsigned& _len2)
{
  
  T max_temp=_m1[0][0];
  //go through the loop to get the max value
  for(unsigned i=0;i<_len1;i++)
    {
      for(unsigned j=0;j<_len2;j++)
	{
	  if(max_temp<_m1[i][j])
	    {
	      max_temp=_m1[i][j];
	    }
	}
    }
  return max_temp;
}

template
double max_mf2(const double* const* _m1, const unsigned& _len1, 
	  const unsigned& _len2);
template
unsigned max_mf2(const unsigned* const* _m1, const unsigned& _len1, 
	  const unsigned& _len2);


//this is the function to find max value among the 2 D array
//input: T** _m1: 2D array pointer
//       _len1: is the size of the first diment
//      unsigned* _m1_size is the size of the second dimenstion.
//           by doing this, we allow non_square arrays
template<class T>
T max_mf2(const T* const* _m1, const unsigned& _len1, 
		const unsigned* _m1_size)
{
  
  T max_temp=_m1[0][0];
  //go through the loop to get the max value
  for(unsigned i=0;i<_len1;i++)
    {
      for(unsigned j=0;j<_m1_size[i];j++)
	{
	  if(max_temp<_m1[i][j])
	    {
	      max_temp=_m1[i][j];
	    }
	}
    }
  return max_temp;
}

template double max_mf2<double>(const double* const* _m1, const unsigned& _len1, const unsigned* _m1_size);
template unsigned max_mf2<unsigned>(const unsigned* const* _m1, const unsigned& _len1, const unsigned* _m1_size);


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
 * 
 * _index is the index of element of the _a array, this is required only when
 * the user want to return a sorted index of the _a array. it has to be initialized 
 * and prefilled with 1 through N in order.!!!!
 *
 */
template <class T> void QuickSort(T* _a, const unsigned& _first, const unsigned int& _last, unsigned* _index, unsigned* _b ) 
{
  //cout<<"&&&&&&&&&&call QuickSort"<<endl;
  unsigned int pivotElement;
 
  if(_first < _last)
    {
      pivotElement = Pivot(_a, _first, _last, _index, _b);
      //cout<<"pivotE,ement:"<<pivotElement<<endl;
      if(pivotElement>0)
	QuickSort(_a, _first, pivotElement-1, _index, _b);
      if(pivotElement< numeric_limits<unsigned int>::max())
	QuickSort(_a, pivotElement+1, _last, _index, _b);
      /*else
	{
	  cout<<"*************WARNING: OVERFLOW for the 
	  }*/
    }
}
 
template  void QuickSort<double>(double* _a,const unsigned& _first, const unsigned int& _last, unsigned* _index, unsigned* _b ); 
template  void QuickSort<unsigned>(unsigned* _a,const unsigned& _first, const unsigned int& _last, unsigned* _index, unsigned* _b ); 
/**
 * Find and return the index of pivot element.
 * @param a - The array.
 * @param first - The start of the sequence.
 * @param last - The end of the sequence.
 * @return - the pivot element
 */
//here we use a new way to choose pivot. a median between [first], [middle] and [last]
template <class T> unsigned Pivot(T* _a, const unsigned int& _first, const unsigned int& _last, unsigned* _index, unsigned* _b) 
{
  //cout<<"%%%%%%%%%%%%call pivot:[_last, _first]:["<<_first<<","<<_last<<"];"<<endl;
  
  unsigned  p = _first;
  unsigned middle=(_first+_last)/2;
  T pivotElement;
  unsigned pivotIndex=_first;
  
  if(middle!=_first&&middle!=_last) //middle==first==middel, we don't do anything, simply use the first one
    {
      //cout<<"\t\t\t[first, last, middle]:["<<_first<<","<<_last<<","<<middle<<"]."<<endl;
      pivotIndex = GetMedianIndex(_a, _first, _last, middle, _index) ;
      //cout<<"\t\t\t>>pivotIndex:"<<pivotIndex<<endl;
      //swap the index
      if(pivotIndex!=_first)
	{
	  //cout<<"#$$$$$$$$$$$$$$$$$$ SWAP the Median Starting Value:"<<endl;
	  swap(_a[_first], _a[pivotIndex]);
	  if(_index!=NULL)
	    {
	      swap(_index[_first], _index[pivotIndex]);
	    }
	  if(_b!=NULL)
	    {
	      swap(_b[_first], _b[pivotIndex]);
	      
	    }
	  pivotIndex=_first;
	}
    }
  //we should have already put the pivot element in the [_first] no matter whether we get the median
  //just use the first one as pivot anyway
  pivotElement=_a[_first];
  pivotIndex=_first;
  //cout<<"8888888888888>>pivote value:"<<_a[_first]<<endl;
  /*if(_index!=NULL)
    {
      _index[_first]=pivotIndex;
      }
  */
  //Print(_a, _last-_first+1);
  for(unsigned int i = _first+1 ; i <= _last ; i++)
    {
      //cout<<"\tround i:"<<i<<"--";
      
      /* If you want to sort the list in the other order, change "<=" to ">" */
      if( _a[i] < pivotElement ||
	  (_a[i]==pivotElement&&(_b!=NULL&&_b[i]<_b[pivotIndex]))||
	  (_a[i]==pivotElement&&_b==NULL&&_index!=NULL&&_index[i]<_index[pivotIndex])
	 )
        {
	  
	  p++;
	  //cout<<"p:"<<p<<";i:"<<i<<endl;
	  if(p!=i)
	    {
	      swap(_a[i], _a[p]);
	      //cout<<"********************SWAP*************!!!"<<endl;
	      if(_index!=NULL)
		{
		  swap(_index[p], _index[i]);
		}
	      if(_b!=NULL)
		{
		  swap(_b[p], _b[i]);
		}
	    }
	  else
	    {
	      //cout<<"******swap step but no swap actions, since p==i"<<endl;
	    }
        }
     
      //Print(_a, _last+1,_index);
    }
  //cout<<"\t>>before pivot"<<endl;
  
  //Print(_a, _last+1,_index);
  swap(_a[p], _a[_first]);
  if(_index!=NULL)
    {
      swap(_index[p],_index[_first]);
    }
  if(_b!=NULL)
    {
      swap(_b[p],_b[_first]);
    }
  //cout<<"\t>>after"<<endl;
  
  //Print(_a, _last+1, _index);
  //cout<<"********end of pivot, p :"<<p<<endl;
  return p;
}
template unsigned Pivot<double>(double* _a, const unsigned int& _first, const unsigned int& _last, unsigned* _index, unsigned* _b) ;
template unsigned Pivot<unsigned>(unsigned* _a, const unsigned int& _first, const unsigned int& _last, unsigned* _index, unsigned* _b) ;
//return 0,1,2 indicating which one is the median
//input:
// m is the array,
// a, b, c is the index of the element that need to be considered
// return index of the median value
template <class T> unsigned GetMedianIndex(const T* m, const unsigned& a, const unsigned& b, const unsigned& c, const unsigned* _index)
{
  if((m[a]==m[b])&&(m[a]==m[c]))
    {
      if(_index!=NULL)
	{
	  if(_index[a]<_index[b]&&_index[a]<_index[c])
	    return a;
	  if(_index[b]<_index[a]&&_index[b]<_index[c])
	    return b;
	  if(_index[c]<_index[a]&&_index[c]<_index[b])
	    return c;
	}
      else
	{
	  if(a<b &&a<c)
	    return a;
	  if(b<a&&b<c)
	    return b;
	  if(c<a&&c<b)
	    return c;
	}
    }
  if(m[a]==m[b])
    {
      if(_index!=NULL)
	{
	  if(_index[a]<_index[b])
	    return a;
	  else
	    return b;
	}
      else
	return a<b?a:b;
    }
  if(m[a]==m[c])
    {
      if(_index!=NULL)
	{
	  if(_index[a]<_index[c])
	    return a;
	  else
	    return c;
	}
      else
	return a<c?a:c;
    }
  if(m[b]==m[c])
    {
      if(_index!=NULL)
	{
	  if(_index[b]<_index[c])
	    return b;
	  else
	    return c;
	}
      else
	return b<c?b:c;
    }
  //for unequal cases
  //unsigned index=0;
  if(m[a]>m[b])
    {
      if(m[c]>m[a]) //a is the median
	{
	  return a;
	}
      else //c<a still not sure, compare b and c
	{
	  if(m[b]<m[c])//b<c
	    {
	      return c;
	    }
	  else //b>c
	    {
	      return b;
	    }	  
	}	
    }
  else //a<b
    {
      if(m[c]<m[a])
	{
	  return a;
	}
      else //c>a need to compare b and c
	{	    
	  if(m[c]>m[b]) //c>b
	    {
	      return b;
	    }
	  else //c<b
	    {	      
	      return c;
	    }
	}      
    }
//  return a;
}
template unsigned GetMedianIndex<double>(const double* m, const unsigned& a, const unsigned& b, const unsigned& c, const unsigned* _index); 
template unsigned GetMedianIndex<unsigned>(const unsigned* m, const unsigned& a, const unsigned& b, const unsigned& c, const unsigned* _index); 
/**
 * Swap the parameters.
 * @param a - The first parameter.
 * @param b - The second parameter.
 */
template <class T> void Swap(T& _a, T& _b)
{
  int temp = _a;
  _a = _b;
  _b = temp;
}
template  void Swap<unsigned>(unsigned& _a, unsigned & _b);
template  void Swap<int>(int& _a, int & _b);
template  void Swap<double>(double& _a, double & _b);
/**
 * Print an array.
 * @param a - The array.
 * @param N - The size of the array.
 */
void Print(const double* a, const int& N, const unsigned* _index)
{
  for(int i = 0 ; i < N ; i++)
    {
      cout << "array[" << i << "] = " << a[i]; 
      if(_index!=NULL)
	{
	  cout<<"-"<<_index[i];
	}
      cout<< endl;
    }
}

 
/**
 * Print an array.
 * @param a - The array.
 * @param N - The size of the array.
 */
void Print(const unsigned* a, const int& N, const unsigned* _index)
{
  for(int i = 0 ; i < N ; i++)
    {
      cout << "array[" << i << "] = " << a[i];// <<"-"<<_index[i]<< endl;
      if(_index!=NULL)
	{
	  cout<<"-"<<_index[i];
	}
      cout<< endl;
    }
} 
 

void Reverse( unsigned* _array, const unsigned& _N)
{
  unsigned start=0,end=_N-1;
  unsigned temp;
  for(;start<end;)
    {
      temp=_array[end];
      _array[end]=_array[start];
      _array[start]=temp;
      start++;
      end--;
    }

}
void Reverse(double* _array, const unsigned& _N)
{
  unsigned start=0,end=_N-1;
  double temp;
  for(;start<end;)
    {
      temp=_array[end];
      _array[end]=_array[start];
      _array[start]=temp;
      start++;
      end--;
    }

}
void Reverse(vector<unsigned> _v)
{
  unsigned start=0,end=_v.size()-1;
  unsigned temp;
  for(;start<end;)
    {
      temp=_v.at(end);
      _v[end]=_v[start];
      _v[start]=temp;
      start++;
      end--;
    }
}
void Reverse(vector<double> _v)
{
  //unsigned size=_v.size();
  unsigned start=0,end=_v.size()-1;
  double temp;
  for(;start<end;)
    {
      temp=_v.at(end);
      _v[end]=_v[start];
      _v[start]=temp;
      start++;
      end--;
    }

}

void Reverse(unsigned* _v, unsigned _size)
{
  //unsigned size=_v.size();
  unsigned start=0,end=_size-1;
  double temp;
  for(;start<end;)
    {
      temp=_v[end];
      _v[end]=_v[start];
      _v[start]=temp;
      start++;
      end--;
    }  

}

/*CopyElements only return false when the caller
 *messes up the sizes of the input arrays. 
 *so if we are calling it right, there will not be errors
 */
/*we want to copy elements from _source to _target based on the indice in the _indexOfElementToCopy
 *target has to be initialized by the caller. The size of target has to be at least the size of the index array
 *
 *
 */
bool CopyElements(const unsigned* _source, const unsigned& _s_size, unsigned* _target, const unsigned& _t_size,
		  const unsigned* _indexOfElementToCopy, const unsigned& _i_size)
{
  //cout<<"**&&^^^%%%%inside copy"<<endl;
  //first, need to check to make sure the target size has to be larger or equal to index array
  if (_t_size<_i_size||_s_size<_i_size)
    return false;
  //cout<<"loop before"<<endl;
  for(unsigned i=0;i<_i_size;i++)
    {
      //cout<<"\tloop "<<i<<endl;
      _target[i]=_source[_indexOfElementToCopy[i]];
      //cout<<"\t\tend loop"<<endl;
    }
  return true;
}

bool CopyElements
  (unsigned** _source, const unsigned& _s_size1, const unsigned& _s_size2, 
   unsigned** _target, const unsigned& _t_size1, const unsigned& _t_size2,
   const unsigned* _indexOfElementToCopy, const unsigned& _i_size)
{
  //first, need to check to make sure the target size has to be larger or equal to index array
  if (_t_size1<_i_size||_s_size1<_i_size)
    return false;
  if(_t_size2<_s_size2)
    return false;
  
  for(unsigned i=0;i<_i_size;i++)
    {
      for(unsigned j=0;j<_s_size2;j++)
	{
	  _target[i][j]=_source[_indexOfElementToCopy[i]][j];
	}
    }
  return true;
}

void Unique(const unsigned* _in, const unsigned& _iSize, 
	    /*output*/ unsigned* _out, unsigned* _out_index, 
	     unsigned& _oSize)
{
  //prepare the index, we assuming in the beginning, the index
  //is in order
  for(unsigned i=0;i<_iSize;i++)
    {
      _out_index[i]=i;
    }
  //first copy over the output array
  std::memcpy(_out, _in, sizeof(unsigned)/sizeof(char)*_iSize);
  //we first sort the array
  QuickSort<unsigned>(_out, 0, _iSize-1, _out_index);
  /*cout<<"After sorting:"<<endl;
  for(unsigned i=0;i<_iSize;i++)
    {
      cout<<_out[i]<<"-"<<_out_index[i]<<",";
    }
    cout<<endl;*/
  //now we need to go through the array to pick the unique ones only
  _oSize=0;

  unsigned runningValue=_out[0];
  //cout<<"initially runningValue:"<<runningValue;  
  _out_index[_oSize]=_out_index[0];
  
  for(unsigned i=1;i<_iSize; i++)
    {
      if(_out[i]!=runningValue)
	{
	  _oSize++;
	  runningValue=_out[i];
	  _out_index[_oSize]=_out_index[i];
	  
	}
    }
  _oSize++;//one more because it is the size, not index
  //we are done with index, just need to populate the output array
  for(unsigned i=0;i<_oSize;i++)
    {
      _out[i]=_in[_out_index[i]];
    }

  //done!!!
} 




