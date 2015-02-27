#include <stdexcept>
#include <exception>
#include <iostream>
#include "Matrix.hpp"

using namespace std;

//constructor
template<class T>
Matrix<T>::Matrix():c_dim(-1),c_dim_size(NULL),c_data(NULL)
{
  //default constructor

}

//data input now is one data.
//we will restructure it
//could be zero dimension, which is like a scalar
template<class T>
Matrix<T>::Matrix(const unsigned& _dim, unsigned _dim_size[], T _data[]):
  c_dim(_dim)
{
  c_dim_size=new unsigned[c_dim];
  unsigned totalNumData=1;
  //need to copy over the data and size
  for(unsigned i=0;i<c_dim;i++)
    {
      if(_dim_size[i]==0)
	{
	  this->c_dim=-1;
	  delete[] c_dim_size;
	  c_dim_size=NULL;
	  cerr<<"constructing an empty/unintialized matrix"<<endl;
	  return;
	}
      c_dim_size[i]=_dim_size[i];
      totalNumData*=c_dim_size[i];
    }
  
  //now we need to figure out the data structure later,
  //for now we say things in 1D
  c_data=new T[totalNumData];
  for(unsigned i=0;i<totalNumData;i++)
    {
      c_data[i]=_data[i];
    }
  //done
}


//data input now is 1D data.
//we will restructure it
//could be zero dimension, which is like a scalar
template<class T>
Matrix<T>::Matrix(const unsigned& _dim, const unsigned* _dim_size, const T* _data):
  c_dim(_dim)
{
  c_dim_size=new unsigned[c_dim];
  unsigned totalNumData=1;
  //need to copy over the data and size
  for(unsigned i=0;i<c_dim;i++)
    {
      if(_dim_size[i]==0)
	{
	  this->c_dim=-1;
	  delete[] c_dim_size;
	  c_dim_size=NULL;
	  cerr<<"constructing an empty/unintialized matrix"<<endl;
	  return;
	}
      c_dim_size[i]=_dim_size[i];
      totalNumData*=c_dim_size[i];
    }
  
  //now we need to figure out the data structure later,
  //for now we say things in 1D
  c_data=new T[totalNumData];
  if(_data!=NULL)
    {
      memcpy(c_data, _data, sizeof(T)*totalNumData);
    }
  else
    {//in case of null data input, set the default value to zero
      std::memset(c_data, 0, sizeof(T)*totalNumData);
    }
  /*  for(unsigned i=0;i<totalNumData;i++)
    {
      
      c_data[i]=_data[i];
      }*/
  //done
}

//destructor
template<class T>
Matrix<T>::~Matrix()
{
  
  if(c_dim_size!=NULL)
    {
      delete[] c_dim_size;
      c_dim_size=NULL;
    }
  if(c_data!=NULL)
    {
      delete[] c_data;
      c_data=NULL;
    }
}

//copy constructor
template<class T>
Matrix<T>::Matrix(const Matrix<T>& _m):
  c_dim(_m.c_dim)
{
  if((signed)(_m.c_dim)==-1)
    {
      //this is an unitialized one,
      //we will keep everything unintialized for this current one
      c_dim_size=NULL;
      c_data=NULL;
      return;
    }
  c_dim_size=new unsigned[c_dim];
  unsigned totalNumData=1;
  //need to copy over the data and size
  for(unsigned i=0;i<c_dim;i++)
    {
      c_dim_size[i]=_m.c_dim_size[i];
      totalNumData*=c_dim_size[i];
    }

  //now we need to figure out the data structure later,
  //for now we say things in 1D
  c_data=new T[totalNumData];
  for(unsigned i=0;i<totalNumData;i++)
    {
      c_data[i]=_m.c_data[i];
    }
  //done
}

//assignment operator
template<class T>
Matrix<T>& Matrix<T>::operator = (const Matrix<T>& _m)
{
  if(this==&_m)
    return *this;

  //deep copy
  //first check to see whether the one is initialized or not
  if((signed)(_m.c_dim)==-1)
    {
      //a uninitialized one
      this->c_dim=-1;
      this->c_dim_size=NULL;
      this->c_data=NULL;
      return *this;
    }
  //a good one with data, could be zero dim'ed (means a scalar like matrix);
  c_dim=_m.c_dim;
  c_dim_size=new unsigned[c_dim];
  unsigned totalNumData=1;
  //need to copy over the data and size
  for(unsigned i=0;i<c_dim;i++)
    {
      c_dim_size[i]=_m.c_dim_size[i];
      totalNumData*=c_dim_size[i];
    }
  
  //now we need to figure out the data structure later,
  //for now we say things in 1D
  c_data=new T[totalNumData];
  for(unsigned i=0;i<totalNumData;i++)
    {
      c_data[i]=_m.c_data[i];
    }
  //done
  return *this;
}

//subscript
//scalar-like matrix
template<class T>
T& Matrix<T>::operator ()()
{
  if(this->c_dim!=0)
    {
      throw runtime_error("unsupported array subscription (dimension not compatible), please check your data! (NOTE:this empty access mehtod only supported on a scalar likematrix)");
    }
  //check for out bound accessing
  
  return this->c_data[0];
}

  
//vector
template<class T>
T& Matrix<T>::operator ()(const unsigned& _d0)
{
  if(this->c_dim!=1)
    {
      throw runtime_error("unsupported array subscription (dimension not compatible), please check your data!");
    }
  //check for out bound accessing
  if(_d0>=c_dim_size[0])
    {
      throw std::out_of_range("index out of range"); 
    }
  //good return 
  return this->c_data[_d0];
}
  
//2d
template<class T>
T& Matrix<T>::operator ()(const unsigned& _d0, const unsigned& _d1)
{
  if(this->c_dim!=2)
    {
      throw std::runtime_error("unsupported array subscription (dimension not compatible), please check your data!");
    }
  //check for out range
  if(_d0>=c_dim_size[0]||_d1>=c_dim_size[1])
    {
      throw std::out_of_range("index out of range");  
    }
  return this->c_data[_d0*c_dim_size[1]+_d1];
  //return temp;
} 

//3d
template<class T>
T& Matrix<T>::operator ()(const unsigned& _d0, const unsigned& _d1, const unsigned& _d2)
{
  if(this->c_dim!=3)
    {
      throw runtime_error("unsupported array subscription (dimension not compatible), please check your data!");
    }
  //check for out range
  if(_d0>=c_dim_size[0]||_d1>=c_dim_size[1]||_d2>=c_dim_size[2])
    {
      throw std::out_of_range("index out of range");  
    }
  return this->c_data[_d0*(_d1*c_dim_size[2])+_d2];
  //return this->temp;
}
  
//4d
template<class T>
T& Matrix<T>::operator ()(const unsigned& _d0, const unsigned& _d1, const unsigned& _d2, const unsigned& _d3)
{
  if(this->c_dim!=4)
    {
      throw std::runtime_error("unsupported array subscription (dimension not compatible), please check your data!");
    }
  //check for out range
  if(_d0>=c_dim_size[0]||_d1>=c_dim_size[1]||_d2>=c_dim_size[2])
    {
      throw std::out_of_range("index out of range");  
    }
  return this->c_data[_d0*_d1*(_d2*c_dim_size[3])+_d3];
  //return this->temp;
}
  
  //NOTE: only support 4d or below

//==========================
  //other operator, so all .dot operation
  //add
template<class T>
Matrix<T> Matrix<T>::operator + (const T& _t)
{
  if((signed)(this->c_dim)==-1)
    {
      throw std::runtime_error("unitialized matrix");
    }
  Matrix<T> temp(*this);
  //determine total number of elements
  unsigned totalNumOfData=1;
  for(unsigned i=0;i<this->c_dim;i++)
    {
      totalNumOfData*=this->c_dim_size[i];
    }
  for(unsigned i=0;i<totalNumOfData;i++)
    {
      temp.c_data[i]=temp.c_data[i]+_t;
    }
  return temp;
}

//subtract
template<class T>
Matrix<T> Matrix<T>::operator - (const T& _t)
{
  if((signed)(this->c_dim)==-1)
    {
      throw std::runtime_error("unitialized matrix");
    }
  Matrix<T> temp(*this);
  //determine total number of elements
  unsigned totalNumOfData=1;
  for(unsigned i=0;i<this->c_dim;i++)
    {
      totalNumOfData*=this->c_dim_size[i];
    }
  for(unsigned i=0;i<totalNumOfData;i++)
    {
      temp.c_data[i]=temp.c_data[i]-_t;
    }
  return temp;
}

//multiplication
template<class T>
Matrix<T>  Matrix<T>::operator * (const T& _t)
{
  if((signed)(this->c_dim)==-1)
    {
      throw std::runtime_error("unitialized matrix");
    }
  Matrix<T> temp(*this);
  //determine total number of elements
  unsigned totalNumOfData=1;
  for(unsigned i=0;i<this->c_dim;i++)
    {
      totalNumOfData*=this->c_dim_size[i];
    }
  for(unsigned i=0;i<totalNumOfData;i++)
    {
      temp.c_data[i]=temp.c_data[i]*_t;
    }
  return temp;
}
//division
template<class T>
Matrix<T>  Matrix<T>::operator / (const T& _t)
{
  if((signed)(this->c_dim)==-1)
    {
      throw std::runtime_error("unitialized matrix");
    }
  Matrix<T> temp(*this);
  //determine total number of elements
  unsigned totalNumOfData=1;
  for(unsigned i=0;i<this->c_dim;i++)
    {
      totalNumOfData*=this->c_dim_size[i];
    }
  for(unsigned i=0;i<totalNumOfData;i++)
    {
      temp.c_data[i]=temp.c_data[i]/_t;
    }
  return temp;
}
/*  
//size
template<class T>
Matrix<unsigned> Matrix<T>::size()
{
  if((unsigned)(this->c_dim)==-1)
    cerr<<"calling on an unitialized matrix object"<<endl;
  unsigned temp_array[] = {this->c_dim};
  Matrix<unsigned> temp(1, temp_array, NULL);
  
  for(unsigned i=0;i<this->c_dim;i++)
    {
      temp.c_data[i]=this->c_dim_size[i];
    }
  return temp;
}
*/
//size
//_dim starting at index
template<class T>
unsigned Matrix<T>::size(unsigned _dim)
{
  if((signed)(this->c_dim)==-1)
    cerr<<"calling on an unitialized matrix object"<<endl;
 
  if(_dim>=this->c_dim)
    {
      throw std::out_of_range("out of range error:dimension is too big");
    }

    return  this->c_dim_size[_dim];

}
  
//
template<class T>
unsigned Matrix<T>::dim()
{
  if((signed)(this->c_dim)==-1)
    cerr<<"calling on an unitialized matrix object"<<endl;
  return this->c_dim;
}
  //
//template<class T>
//unsigned Matrix<T>::dim_size(const unsigned& _dim)
//{
//}

//get submatrix
//here we did not pass the _dim_pos[] by reference
//that means we don't have to be const anyway.
//if we want to pass by reference, we simply pass by pointer
//
//Note: currently don't support resize/reshape, _n (the size of the input array _dim_pos[])
//		has to be the same as the matrix dim. 
//		In the future, we will implement reshape and the size of the input array 
//       can be different from matrix dim.
//		the dim of the output array is determined by the _dim_pos array
//      the input position array starting from zero and can be -1, mean all 
//      elements of this dim are picked, similar to Matlab : notation
//      eg. {-1,2,-1} ~ (:,2,:)
//input:
//
//output: Note: this will return a Matrix, or a scalar like Matrix
//		in any case
template<class T>
Matrix<T> Matrix<T>::SubMatrix(const unsigned& _n, int _dim_pos[])
{
	//first we need decide whether the input are valid
	if((signed)(this->c_dim)==-1)
	{
	  throw std::runtime_error("calling the getsubmatrix on a uninitialized matrix.");
	}
	if(_n!=this->dim())
	{
	  throw std::runtime_error ("unsupport format, the input array has to be same size as the matrix");
	}
	
	for(unsigned i=0;i<_n;i++)
	{
	  if(_dim_pos[i]>0&&((unsigned)_dim_pos[i])>=this->size(i))
		  throw std::out_of_range("dimension size is out of range");	
	}
	
	//start doing the job
	//first determine the output matrix dimension
	unsigned new_dim=0;
	unsigned* new_dim_size;
	for(unsigned i=0;i<_n;i++)
	{
		if(_dim_pos[i]>0)
		{
			new_dim++;
		}
	}
	//determine demension size
	new_dim_size=new unsigned[new_dim];
	unsigned counter=0;
	for(unsigned i=0;i<_n;i++)
	{
		if(_dim_pos[i]>0)
		{
			new_dim_size[counter]=_dim_pos[i];
			counter++;
		}
	}
	Matrix<T> temp_m(new_dim, new_dim_size, NULL);
	
	//first determine block number and offset and block size
	unsigned* blockNumber=new unsigned[new_dim];
	unsigned* blockSize_index_in_dim_size=new unsigned[new_dim];
	unsigned* offset=new unsigned[new_dim];
	unsigned* blockSize=new unsigned [new_dim];

	for(unsigned i=0;i<new_dim;i++)
	{
		blockNumber[i]=0;
		blockSize_index_in_dim_size[i]=1;
		offset[i]=0;
		blockSize[i]=1;
	}
	unsigned curr_block_index=0;
	for(unsigned i=0;i<_n;i++)
	{
		if(_dim_pos[i]<0) //-1 means all
		{
			blockNumber[curr_block_index]=new_dim_size[curr_block_index];
			blockSize_index_in_dim_size[curr_block_index]=i+1;
			//curr_block_index++;
		}
		else //positive means very specific element in this dimension
		{
			offset[curr_block_index]+=_dim_pos[i]*this->c_dim_size[i];
		}
		if(_dim_pos[i]<0)
		{
			curr_block_index++;
		}
	}
	//now we need to determine the block size based on the blockSize_index_in_dim_size
	//unsigned* blockSize=new unsigned [new_dim];
	for(unsigned i=0;i<new_dim;i++)
	{
		for(unsigned j=blockSize_index_in_dim_size[i];j<this->c_dim;j++)
		{
			blockSize[i]*=this->c_dim_size[j];
		}
	}
	
	//now we kind of get everything ready and just need to copy over the data
	unsigned totalNumOfOutputData=1;
	for(unsigned i=0;i<new_dim;i++)
	{
		totalNumOfOutputData*=new_dim_size[i];
	}
	unsigned* data_index=new unsigned[totalNumOfOutputData];//used to rember where the data will be copied to the data out array

	unsigned sizeOfCurrentBlock_running=totalNumOfOutputData;
	
	for(unsigned i=0;i<new_dim;i++)
	{
		
		for(unsigned j=0;j<blockNumber[i];j++)
		{
			sizeOfCurrentBlock_running=sizeOfCurrentBlock_running/new_dim_size[i];
			for(unsigned k=0;k<sizeOfCurrentBlock_running;k++)
			{
				data_index[k+j*sizeOfCurrentBlock_running]+=j*blockSize[i]+offset[i];
			}	
		}
	}

	//now based on the indices, copy over the element
	for(unsigned i=0;i<totalNumOfOutputData;i++)
	  {
	    temp_m.c_data[i]=this->c_data[data_index[i]];
	  }

	//done!!
	return temp_m;
}
/*  Matrix<T> GetSubMatrix(int d1, int d2);

  Matrix<T> GetSubMatrix(int d1, int d2, int d3);

  Matrix<T> GetSubMatrix(int d1, int d2, int d3, int d3);
*/

//c++ equivalent to Matlab sum(A), return
//the sum of the elements of A along the first array 
//dimension whose size does not equal 1:
template<class T>
Matrix<T> sum(const Matrix<T>& _m)
{
  if((signed)(_m.dim())==-1)
    {
      throw std::runtime_error("calling on a uninitialized Matrix object for sum");
    }
  if(_m.dim()==0)
    {
      return _m;
    }
  //we decide which dimension to return
  for(unsigned int i=0;i<_m.dim();i++)
    {
      if(_m.size(i)>1)
	{
	  return sum(_m, i);
	}
      
    }
  //if we are here the dimension is 1 for all dimensions
  //so return the original matrix
  return _m;
}

//c++ equivalent to Matlab sum(A,dim), return
//sums the elements of A along dimension dim. 
//The dim input is a Matrix.
//dim starts at zero
template<class T>
Matrix<T> sum(const Matrix<T>& _m, const unsigned& _dim)
{
  if((signed)(_m.dim())==-1)
    throw runtime_error("call on an unitialized Matrix object for sum()");

  if(_m.dim()==0)
    {
      cerr<<"calling on an scalar like matrix"<<endl;
      //no sum to do
      return _m;
    }
  //first determine dimension and size of each dim
  if(_dim>=_m.dim())
    {
      throw std::out_of_range("dimension is out of range");
    }
  unsigned new_dim=_m.dim()-1;
  unsigned* new_dim_size=new unsigned [new_dim];
  unsigned totalNumberOfElements=1;
  unsigned new_index_counter=0;
  for(unsigned i=0;i<new_dim;i++)
    {
      if(i!=_dim)
	{
	  new_dim_size[new_index_counter]=_m.size(i);
	  new_index_counter++;
	}
    }
  Matrix<T> temp_m(new_dim, new_dim_size, NULL);

  //start getting the sums
  //first need to figuire out how many blocks we need
  T* p_firstElementEachBlock_dst;
  T* p_firstElementEachBlock_src;
  
  unsigned numberOfBlocks=1;
  unsigned sizeOfEachBlock=1;
  bool flag=false;//indicating whether the deterimin the block number action is done
  //start to determine the numOfBlocks
  for(unsigned i=0;i<_m.dim();i++)
    {
      if(i==_dim) //we are doing with determine the block number
	flag=true;//we are done
      
      if(!flag)
	{
	  numberOfBlocks*=_m.size(i);
	}
      else
	{//we skip this current if this is the flag round, since it will collapsed on this _dim dimension
	  if(i!=_m.c_dim-1)
	    sizeOfEachBlock *= _m.size(i+1);
	}

      totalNumberOfElements*=_m.size(i);
    }
  
  //temp
  //now we need to determine the sizeOfEachBlock
  p_firstElementEachBlock_dst=temp_m.c_data; //starting point destination
  p_firstElementEachBlock_src=_m->c_data;//starting point source
  for(unsigned i=0;i<numberOfBlocks;i++)
    {
      p_firstElementEachBlock_dst += i*sizeOfEachBlock;
      p_firstElementEachBlock_src += i*(sizeOfEachBlock*_m.size(_dim));
      
      for(unsigned k=0;k<sizeOfEachBlock;k++)
	{
	  p_firstElementEachBlock_dst[k]=0;
	  for(unsigned j=0;j<_m.size(_dim);j++)
	    {
	      p_firstElementEachBlock_dst[k]+=p_firstElementEachBlock_src[k+j*sizeOfEachBlock];
	    }
	}//end k for loop
    }//end of i for loop

  return temp_m;
  //done here
}

//======template instantiation
//class
template class Matrix<int>;
template class Matrix<unsigned>;
//template class Matrix<bool>;
template class Matrix<double>;
template class Matrix<float>;
//template class Matrix<byte>;

//function
template Matrix<int> sum(const Matrix<int>& _m);
template Matrix<unsigned> sum(const Matrix<unsigned>& _m);
template Matrix<float> sum(const Matrix<float>& _m);
template Matrix<double> sum(const Matrix<double>& _m);
/*

//function
template Matrix<int> sum(const Matrix<int>& _m, const unsigned& _dim);
template Matrix<unsigned> sum(const Matrix<unsigned>& _m, const unsigned& _dim);
template Matrix<float> sum(const Matrix<float>& _m, const unsigned& _dim);
template Matrix<double> sum(const Matrix<double>& _m, const unsigned& _dim);
*/
