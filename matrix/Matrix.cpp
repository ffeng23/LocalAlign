#include <stdexcept>
#include <exception>
#include <iostream>
#include <string>
#include <sstream>
//#include <stdio.h>
#include <cstring>
#include <stdlib.h>

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
Matrix<T>::Matrix(const unsigned& _dim, const unsigned _dim_size[], const T& _data):
  c_dim(_dim)
{
  //cout<<"#####build with array the matrix"<<endl;
  //for now only support 
  if(c_dim>4)
    {
      cerr<<"ERROR: so far, matrix larger than 4 dimensions is not support, quit...."<<endl;
      exit(-1);
    }
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
  for(unsigned int i=0;i<totalNumData;i++)
    {
      c_data[i]=_data;
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
  //cout<<"^^^^^^^^building an object"<<endl;
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
  //cout<<"*******calling destructor*********"<<endl;
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

  //now we need to check to see whether the current one is empty or not.
  //we could allow the current to be a filled one, but just need to
  //delete/clear the data first
  if((signed)(this->c_dim)!=-1)
    {
      if(c_dim_size!=NULL)
	{
	  delete[] c_dim_size;
	  c_dim_size=NULL;
	}
      if(c_data!=NULL)
	{
	  delete [] c_data;
	  c_data=NULL;
	}
    }

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
  memcpy(c_data, _m.c_data, sizeof(T)*totalNumData);
  
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
      cerr<<"unsupported array subscription (dimension not compatible), please check your data! exception thrown (NOTE:this empty access mehtod only supported on a scalar likematrix)"<<endl;
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
      cerr<<"unsupported array subscription (dimension not compatible), please check your data! exception thrown"<<endl;
      throw runtime_error("unsupported array subscription (dimension not compatible), please check your data!");
    }
  //check for out bound accessing
  if(_d0>=c_dim_size[0])
    {
      cerr<<"index out of range, exception thrown"<<endl;
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
      cerr<<"unsupported array subscription (dimension not compatible), please check your data! exception thrown"<<endl;
      throw runtime_error("unsupported array subscription (dimension not compatible), please check your data!");
    }
  //check for out range
  if(_d0>=c_dim_size[0]||_d1>=c_dim_size[1]||_d2>=c_dim_size[2])
    {
      cerr<<"index out of range, exception thrown"<<endl;
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
      cerr<<"unsupported array subscription (dimension not compatible), please check your data!"<<endl;
      throw std::runtime_error("unsupported array subscription (dimension not compatible), please check your data!");
    }
  //check for out range
  if(_d0>=c_dim_size[0]||_d1>=c_dim_size[1]||_d2>=c_dim_size[2])
    {
      cerr<<"index out of rangen, exception thrown"<<endl;
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
      cerr<<"unitialized matrix, exception thrown"<<endl;
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
      cerr<<"unitialized matrix, exception thrown"<<endl;
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
      cerr<<"unitialized matrix, exception thrown"<<endl;
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
      cerr<<"unitialized matrix, exception thrown"<<endl;
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
  
//size
template<class T>
Matrix<unsigned> Matrix<T>::size() const
{
  if((signed)(this->c_dim)==-1)
    cerr<<"calling on an unitialized matrix object"<<endl;
  unsigned temp_array[] = {this->c_dim};
  Matrix<unsigned> temp(1, temp_array, (unsigned*)NULL);
  
  for(unsigned i=0;i<this->c_dim;i++)
    {
      temp(i)=this->c_dim_size[i];
    }
  return temp;
}
//*/
//size
//_dim starting at index
template<class T>
unsigned Matrix<T>::size(const unsigned& _dim) const
{
  if((signed)(this->c_dim)==-1)
    cerr<<"calling on an unitialized matrix object"<<endl;
 
  if(_dim>=this->c_dim)
    {
      cerr<<"out of range error:dimension is too big, exception thrown"<<endl;
      throw std::out_of_range("out of range error:dimension is too big");
    }

    return  this->c_dim_size[_dim];

}
  
//
template<class T>
unsigned Matrix<T>::dim() const
{
  if((signed)(this->c_dim)==-1)
    cerr<<"calling on an unitialized matrix object"<<endl;
  return this->c_dim;
}

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
Matrix<T> Matrix<T>::SubMatrix(const unsigned& _n, const int _dim_pos[]) const
{
	//first we need decide whether the input are valid
	if((signed)(this->c_dim)==-1)
	{
	  cerr<<"calling the getsubmatrix on a uninitialized matrix.exception thrown"<<endl;
	  throw std::runtime_error("calling the getsubmatrix on a uninitialized matrix.");
	}
	if(_n!=this->dim())
	{
	  cerr<<"unsupport format, the input array has to be same size as the matrix, exception thrown"<<endl;
	  throw std::runtime_error ("unsupport format, the input array has to be same size as the matrix");
	}
	
	for(unsigned i=0;i<_n;i++)
	{
	  if(_dim_pos[i]>0&&((unsigned)_dim_pos[i])>=this->size(i))
	    {
	      cerr<<"the specified index is out of range of the sie of the dimension in submatrix, please check"<<endl;
	      throw std::out_of_range("dimension size is out of range");
	    }
	}
	//cout<<"do 1"<<endl;	
	//start doing the job
	//first determine the output matrix dimension
	unsigned new_dim=0;
	unsigned* new_dim_size;
	for(unsigned i=0;i<_n;i++)
	{
		if(_dim_pos[i]<0)
		{
			new_dim++;
		}
	}
	//cout<<"do 2"<<endl;	
	//determine demension size
	new_dim_size=new unsigned[new_dim];
	unsigned counter=0;
	unsigned totalNumOfInputData=1;
	for(unsigned i=0;i<_n;i++)
	{
	  totalNumOfInputData*=c_dim_size[i];
		if(_dim_pos[i]<0)
		{
			new_dim_size[counter]=c_dim_size[i];
			counter++;
			
		}
	}
	//cout<<"do 3"<<endl;	
	Matrix<T> temp_m(new_dim, new_dim_size,(T*) NULL);
	
	//first determine block number and offset and block size
	unsigned* blockNumber=new unsigned[new_dim];
	unsigned* blockSize_index_in_dim_size=new unsigned[new_dim];
	unsigned* offset=new unsigned[new_dim+1];
	/*NOTE: offset is a array with new_dim+1 elements!!! one more than the others
	 *the rationale is that the offset[0] is the starting offset, then the dimension i offset
	 *is offset[i+1]. this is necessary when there are trailing offsets.
	 *for example, when there are 1D, the offset is of 2elements. offset[0]
	 *is the starting offset (could be zero). we will add the starting offset to all the elements in
	 * the begining. then for dim 1 (dim[0]), we will add offset[1] (i+1, where i=0).
	 * for another example, when 0D, the offset is of 1 element long, which is the 
	 *the starting offset, will be add to all element in the begining.
	 *for any offset elements, 0 offset is possible. it depends one what is following
	 *to the specific dim. if a dimension is followed with a another dimension, specified by -1, 
	 *then the offset is zero at this level. if followed with a specified one, offset is something
	 *pointed to the elements within this specified block
	 */
	unsigned* blockSize=new unsigned [new_dim];

	for(unsigned i=0;i<new_dim;i++)
	{
		blockNumber[i]=0;
		blockSize_index_in_dim_size[i]=1;
		offset[i]=0;
		blockSize[i]=1;
	}
	offset[new_dim]=0;//since offset is 1+new_dim, see above 
	//cout<<"do 4"<<endl;	
	unsigned curr_block_index=0;
	unsigned enclosedInputBlockSize_running=totalNumOfInputData;
	//bool roundFlag=false;//one round means, we see

		
	for(unsigned i=0;i<_n;i++)
	{
	  //cout<<"i:"<<i<<endl;
	  enclosedInputBlockSize_running/=c_dim_size[i];
	  //cout<<"\tencloseInputblocksize:"<<enclosedInputBlockSize_running<<endl;
		if(_dim_pos[i]<0) //-1 means all
		{
			blockNumber[curr_block_index]=new_dim_size[curr_block_index];
			blockSize_index_in_dim_size[curr_block_index]=i+1;
			//curr_block_index++;
			//cout<<"\tin negative case!!"<<endl;
		}
		else //positive means very specific element in this dimension
		{
		  offset[curr_block_index]+=_dim_pos[i]*enclosedInputBlockSize_running;//_dim_pos[i]*this->c_dim_size[i];
		  //cout<<"\tpositive case:_dim_pos[i]"<<_dim_pos[i]<<endl;
		  
		}
		if(_dim_pos[i]<0)
		{
			curr_block_index++;
		}
	}
	/*
	//now determine the offset
	unsigned current_offset_index=0;
	for(unsigned i=0;i<_n;i++)
	  {
	    //the first offset
	    if(_dim_size[0]<0) //-1 means all, offset is zero
	      {
		
	      }
	    else  //positive menas very specific element
	      {
		offset[curr_block_index]+=_dim_pos[i]*enclosedInputBlockSize_running;//_dim_pos[i]*this->c_dim_size[i];
		  cout<<"\tpositive case:_dim_pos[i]"<<_dim_pos[i]<<endl;
	      }
	  }
	*/
	//cout<<"do 5"<<endl;	
	//now we need to determine the block size based on the blockSize_index_in_dim_size
	//unsigned* blockSize=new unsigned [new_dim];
	unsigned totalNumOfOutputData=1;
	
	for(unsigned i=0;i<new_dim;i++)
	{
	  totalNumOfOutputData*=new_dim_size[i];
		for(unsigned j=blockSize_index_in_dim_size[i];j<this->c_dim;j++)
		{
			blockSize[i]*=this->c_dim_size[j];
		}
	}
	//cout<<"do 6"<<endl;	
	//now we kind of get everything ready and just need to copy over the data
	
	unsigned* data_index=new unsigned[totalNumOfOutputData];//used to rember where the data will be copied to the data out array

	//cout<<"do 7"<<endl;	
	unsigned sizeOfCurrentBlock_running=totalNumOfOutputData;
	//cout<<"new_dim:"<<new_dim<<";totalNumOfOutputData:"<<totalNumOfOutputData<<";blockNumber:"
	//    <<blockNumber<<";blockSize:"<<blockSize<<";offset:"<<offset<<endl;

	//debugging
	/*cout<<"blockNumber:";
	for(unsigned i=0;i<new_dim;i++)
	  {
	    cout<<blockNumber[i]<<",";
	  }
	cout<<endl;
	
	cout<<"blockSize:";
	for(unsigned i=0;i<new_dim;i++)
	  {
	    cout<<blockSize[i]<<",";
	  }
	cout<<endl;

	cout<<"offset:";
	for(unsigned i=0;i<new_dim+1;i++)
	  {
	    cout<<offset[i]<<",";
	  }
	cout<<endl;
	*/
	//cout<<"debugging index determining:"<<endl;
	//the first step is to add the starting offset (offset[0]) to all the entries
	for(unsigned i=0;i<totalNumOfOutputData;i++)
	  {
	    data_index[i]=offset[0];
	  }
	
	
	unsigned totalNumOfOutputBlock_running;
	unsigned totalNumOfOutputBlock_size_running;
	for(unsigned i=0;i<new_dim;i++)
	{
	  //cout<<"i-"<<i<<":";
	   sizeOfCurrentBlock_running=sizeOfCurrentBlock_running/new_dim_size[i];
	   totalNumOfOutputBlock_running=totalNumOfOutputData/sizeOfCurrentBlock_running/new_dim_size[i];
	   totalNumOfOutputBlock_size_running=totalNumOfOutputData/totalNumOfOutputBlock_running;
	   //cout<<"repeat block number:"<<totalNumOfOutputBlock_running<<endl;
	   for(unsigned p=0;p<totalNumOfOutputBlock_running;p++)
	     {
	       //cout<<"sizeofcurrentBlock_running:"<<sizeOfCurrentBlock_running<<":"<<endl;
	       for(unsigned j=0;j<blockNumber[i];j++)
		 {
		   //cout<<"\tj-"<<j<<",";
		   
		   for(unsigned k=0;k<sizeOfCurrentBlock_running;k++)
		    {
		      //cout<<"\t\tk-"<<k<<",entry:[k+j*sizeOfCurrentBlock_running]"<<k+j*sizeOfCurrentBlock_running+p*totalNumOfOutputBlock_size_running;
		      data_index[k+j*sizeOfCurrentBlock_running+p*totalNumOfOutputBlock_size_running]+=j*blockSize[i]+offset[i+1];
		      //cout<<"=>data_index:"<<data_index[k+j*sizeOfCurrentBlock_running+p*totalNumOfOutputBlock_size_running]<<endl;
		      
		    }	
		}
	     }
	}
	//cout<<"do 8"<<endl;	
	//now based on the indices, copy over the element
	//cout<<"data index:"<<endl;
	for(unsigned i=0;i<totalNumOfOutputData;i++)
	  {
	    temp_m.c_data[i]=this->c_data[data_index[i]];
	    //cout<<data_index[i]<<",";
	  }
	//cout<<endl;
	//	cout<<"do 10"<<endl;	

	//done!!
	return temp_m;
}

//this is  a function only used to initialize/fill in data for an empty matrix
  //
template<class T>
void Matrix<T>::initialize(const unsigned& _dim, const unsigned _dim_size[], const T _data[])
{
  //first we need to check whether this current one is an empty matrix
  if((signed)(this->c_dim)!=-1)
    {
      cerr<<"Run time error:calling to reinitialize on an unempty matrix, exception thrown"<<endl;
      throw runtime_error("Run time error:calling to reinitialize on an unempty matrix");
    }
  
  //now do the data filling
  //for now only support 
  if(_dim>4)
    {
      cerr<<"ERROR: so far, matrix larger than 4 dimensions is not support, quit...."<<endl;
      exit(-1);
    }
  c_dim=_dim;
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
  
  c_data=new T[totalNumData];
  if(_data==NULL)
    {
      //intialize it with zeros
      std::memset(c_data, 0, totalNumData*sizeof(T));
      return;
    }
  
  //now we need to figure out the data structure later,
  //for now we say things in 1D
  
  for(unsigned i=0;i<totalNumData;i++)
    {
      c_data[i]=_data[i];
    }
  //done
}

//this is  a function only used to initialize/fill in data for an empty matrix
  //
template<class T>
void Matrix<T>::initialize(const unsigned& _dim, const unsigned _dim_size[], const T& _data)
{
  //first we need to check whether this current one is an empty matrix
  if((signed)(this->c_dim)!=-1)
    {
      cerr<<"Run time error:calling to reinitialize on an unempty matrix, exception thrown"<<endl;
      throw runtime_error("Run time error:calling to reinitialize on an unempty matrix");
    }
  
  //now do the data filling
  //for now only support 
  if(_dim>4)
    {
      cerr<<"ERROR: so far, matrix larger than 4 dimensions is not support, quit...."<<endl;
      exit(-1);
    }
  c_dim=_dim;
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
  
  c_data=new T[totalNumData];
  
  //now we need to figure out the data structure later,
  //for now we say things in 1D
  
  for(unsigned i=0;i<totalNumData;i++)
    {
      c_data[i]=_data;
    }
  //done
}


/*  Matrix<T> GetSubMatrix(int d1, int d2);

  Matrix<T> GetSubMatrix(int d1, int d2, int d3);

  Matrix<T> GetSubMatrix(int d1, int d2, int d3, int d3);
*/

//used for printing
template<class T>
std::string Matrix<T>::toString() const
{
  //depending on the dimension, we need to do different output
  //first check for a valid/initialized matrix
  stringstream ss;
  if((signed)(this->c_dim)==-1)
    {
      ss<<"empty/uninitialized matrix []\n";
      return ss.str();
    }
  
  //0D, scalar-like matrix
  if(this->c_dim==0)
    {
      ss<<"\t"<<this->c_data[0]<<"\n";
      return ss.str();
    }
  
  //1D, vector
  if(this->c_dim==1)
    {
      for(unsigned i=0;i<c_dim_size[0];i++)
	{
	  ss<<c_data[i];
	  if(i!=c_dim_size[0]-1)
	    {
	      ss<<" ";
	    }
	}
      ss<<"\n";
      return ss.str();
    }

  //2D, matrix
  if(this->c_dim==2)
    {
      for(unsigned i=0;i<c_dim_size[0];i++)
	{
	  for(unsigned j=0;j<c_dim_size[1];j++)
	    {
	      ss<<c_data[i*c_dim_size[1]+j];
	      if(j!=c_dim_size[1]-1)
		{
		  ss<<" ";
		}
	      else
		{
		  ss<<"\n";
		}
	    }
	}
      ss<<"\n";
      return ss.str();
    }


  //3D, matrix
  if(this->c_dim==3)
    {
      for(unsigned i=0;i<c_dim_size[0];i++)
	{
	  //by section here or block
	  ss<<"("<<i<<",:,:):\n";
	  for(unsigned j=0;j<c_dim_size[1];j++)
	    {
	      for(unsigned k=0;k<c_dim_size[2];k++)
		{
		  ss<<c_data[i*c_dim_size[2]*c_dim_size[1]+j*c_dim_size[2]+k];
		  if(k!=c_dim_size[2]-1)
		    {
		      ss<<" ";
		    }
		  else
		    {
		      ss<<"\n";
		    }
		}
	    }
	}
      ss<<"\n";
      return ss.str();
    }  

  //4D, matrix
  if(this->c_dim==4)
    {
      for(unsigned i=0;i<c_dim_size[0];i++)
	{
	  //by section here or block
	  
	  for(unsigned p=0;p<c_dim_size[1];p++)
	    {
	      ss<<"("<<i<<":"<<p<<",:,:):\n";
	      for(unsigned j=0;j<c_dim_size[2];j++)
		{
		  for(unsigned k=0;k<c_dim_size[3];k++)
		    {
		      ss<<c_data[i*c_dim_size[1]*c_dim_size[2]*c_dim_size[3]+
				 p*c_dim_size[2]*c_dim_size[3]+
				 j*c_dim_size[3]+k];
		      if(k!=c_dim_size[3]-1)
			{
			  ss<<" ";
			}
		      else
			{
			  ss<<"\n";
			}
		    }
		}
	    }
	}
      ss<<"\n";
      return ss.str();
    }  
  //done
  cerr<<"Run time error:matrix with more than 4 dimensions is not supported so far, exception thrown"<<endl;
  throw runtime_error("Run time error:matrix with more than 4 dimensions is not supported so far");

  return "";
}



//c++ equivalent to Matlab sum(A), return
//the sum of the elements of A along the first array 
//dimension whose size does not equal 1:
template<class T>
Matrix<T> sum(const Matrix<T>& _m)
{
  if((signed)(_m.dim())==-1)
    {
      cerr<<"call on an unitialized Matrix object for sum(),exception thrown"<<endl;
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
    {
      cerr<<"call on an unitialized Matrix object for sum(), exception thrown"<<endl;
      throw runtime_error("call on an unitialized Matrix object for sum()");
    }
  if(_m.dim()==0)
    {
      cerr<<"calling on an scalar like matrix"<<endl;
      //no sum to do
      return _m;
    }
  //first determine dimension and size of each dim
  if(_dim>=_m.dim())
    {
      cerr<<"dimension is out of range, exception thrown"<<endl;
      throw std::out_of_range("dimension is out of range");
    }
  //cout<<"calling sum 1"<<endl;
  unsigned new_dim=_m.dim()-1;
  unsigned* new_dim_size=new unsigned [new_dim];
  unsigned totalNumberOfElements=1;
  unsigned new_index_counter=0;
  for(unsigned i=0;i<_m.c_dim;i++)
    {
      if(i!=_dim)
	{
	  new_dim_size[new_index_counter]=_m.size(i);
	  //cout<<"dim size:"<<new_dim_size[new_index_counter]<<endl;
	  new_index_counter++;
	  
	}
    }
  //cout<<"calling sum 2, new_dim:"<<new_dim<<",new_dim_size:"<<new_dim_size[0]<<endl;

  Matrix<T> temp_m(new_dim, new_dim_size, (T*)NULL);
  //cout<<temp_m.toString()<<endl;
  //start getting the sums
  //first need to figuire out how many blocks we need
  T* p_firstElementEachBlock_dst;
  T* p_firstElementEachBlock_src;
  
  unsigned numberOfBlocks=1;
  unsigned sizeOfEachBlock_output=1;
  unsigned sizeOfEachBlock_input=1;
  bool flag=false;//indicating whether the deterimin the block number action is done
  //cout<<"calling sum 3"<<endl;

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
	    {
	      sizeOfEachBlock_output *= _m.size(i+1);
	    }
	  sizeOfEachBlock_input*=_m.size(i);
	}

      totalNumberOfElements*=_m.size(i);
    }
  /*
  cout<<"totalNumberOfElements:"<<totalNumberOfElements
      <<"number of blocks:"<<numberOfBlocks
      <<";sizeOfEachBlock_output:"<<sizeOfEachBlock_output
      <<";sizeOfEachBlack_input:"<<sizeOfEachBlock_input<<endl;
  */
  //cout<<temp_m.toString()<<endl;
  //temp
  //now we need to determine the sizeOfEachBlock
  p_firstElementEachBlock_dst=temp_m.c_data; //starting point destination
  p_firstElementEachBlock_src=_m.c_data;//starting point source
  //cout<<"p_first dst:"<<p_firstElementEachBlock_dst<<";p_first src:"<<p_firstElementEachBlock_src<<endl;
  for(unsigned i=0;i<numberOfBlocks;i++)
    {
      //cout<<"==>loop i:"<<i<<endl;
      p_firstElementEachBlock_dst =temp_m.c_data+i* sizeOfEachBlock_output;
      p_firstElementEachBlock_src = _m.c_data+i*sizeOfEachBlock_input;
      //cout<<"\ti*sizeOfEachBlock_output:"<<i*sizeOfEachBlock_output<<endl;
      //cout<<"\ti*sizeOfEachBlock_input:"<<i*sizeOfEachBlock_input<<endl;
      for(unsigned k=0;k<sizeOfEachBlock_output;k++)
	{
	  //cout<<"\tk:"<<k<<endl;
	  p_firstElementEachBlock_dst[k]=0;
	  //cout<<"\t_m.size(_dim):"<<_m.size(_dim)<<endl;
	  for(unsigned j=0;j<_m.size(_dim);j++)
	    {
	      //cout<<"\t\tj:"<<j<<";k+j*sizeOfEachBlock_output:"<<k+j*sizeOfEachBlock_output<<endl;
	      //cout<<"\t\telement:"<<p_firstElementEachBlock_src[k+j*sizeOfEachBlock_output]<<endl;
	      p_firstElementEachBlock_dst[k]+=p_firstElementEachBlock_src[k+j*sizeOfEachBlock_output];
	    }
	}//end k for loop
    }//end of i for loop
  
  return temp_m;
  
}


//this is the one linearize the matrix data and sum all the entries together to 
//get the result
template<class T>
T sum_all(const Matrix<T>& _m)
{
  unsigned dim=_m.dim();
  //unsigned dim
  if(dim==0)
    {
      return _m.c_data[0];
    }
  unsigned totalNumOfData=1;
  for(unsigned i=0;i<dim;i++)
    {
      totalNumOfData *=_m.size(i);
    }
  T ret=0;
  for(unsigned i=0;i<totalNumOfData;i++)
    {
      ret+=_m.c_data[i];
    }
  return ret;
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


//function
template Matrix<int> sum(const Matrix<int>& _m, const unsigned& _dim);
template Matrix<unsigned> sum(const Matrix<unsigned>& _m, const unsigned& _dim);
template Matrix<float> sum(const Matrix<float>& _m, const unsigned& _dim);
template Matrix<double> sum(const Matrix<double>& _m, const unsigned& _dim);

//function
template int sum_all(const Matrix<int>& _m);
template unsigned sum_all(const Matrix<unsigned>& _m);
template float sum_all(const Matrix<float>& _m);
template double sum_all(const Matrix<double>& _m);
//*/
