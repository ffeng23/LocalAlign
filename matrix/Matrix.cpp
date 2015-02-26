#include <stdexcept>
#include <iostream>
#include <Matrix.hpp>

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
Matrix<T>::Matrix(unsigned _dim, unsigned []_dim_size, T[] _data):
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
Matrix<T>::Matrix(unsigned _dim, unsigned* _dim_size, T* _data):
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
*/
//destructor
template<class T>
Matrix<T>::~Matrix()
{
  
  if(c_dim_size!=NULL)
    {
      delete[] c_dim_size;
      c_dim_size=NULL:
    }
  if(c_data!=NULL)
    {
      delete[] c_data;
      c_data=NULL:
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
Matrix<T>::Matrix<T>& operator = (const Matrix<T>& _m)
{
  if(this==&_m)
    return *this;

  //deep copy
  //first check to see whether the one is initialized or not
  if((signed)(_m.c_dim)=-1)
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
T Matrix<T>::operator []()
{
  if(this->c_dim!=0)
    {
      throw new exception("unsupported array subscription (dimension not compatible), please check your data! (NOTE:this empty access mehtod only supported on a scalar likematrix)")
    }
  //check for out bound accessing
  
  return this->c_data[0];
}

  
//vector
template<class T>
T Matrix<T>::operator [](const unsigned& _d0)
{
  if(this->c_dim!=1)
    {
      throw new exception("unsupported array subscription (dimension not compatible), please check your data!")
    }
  //check for out bound accessing
  if(d0>=c_dim_size[0])
    {
      throw std::out_range("index out of range"); 
    }
  //good return 
  return this->c_data[d0];
}
  
//2d
template<class T>
T& Matrix<T>::operator [](const unsigned& _d0, const unsigned& _d1)
{
  if(this->c_dim!=2)
    {
      throw new exception("unsupported array subscription (dimension not compatible), please check your data!")
    }
  //check for out range
  if(d0>=c_dim_size[0]||d1>=c_dim_size[1])
    {
      throw std::out_range("index out of range");  
    }
  T temp=this->c_data[d0*c_dim_size[1]+d1];
  return this->temp;
} 

//3d
template<class T>
T& Matrix<T>::operator [](const unsigned& _d0, const unsigned& _d1, const unsigned& _d2)
{
  if(this->c_dim!=3)
    {
      throw new exception("unsupported array subscription (dimension not compatible), please check your data!")
    }
  //check for out range
  if(d0>=c_dim_size[0]||d1>=c_dim_size[1]||d2>=c_dim_size[2])
    {
      throw std::out_range("index out of range");  
    }
  T temp=this->c_data[d0*(d1*c_dim_size[2])+d2];
  return this->temp;
}
  
//4d
template<class T>
T& Matrix<T>::operator [](const unsigned& _d0, const unsigned& _d1, const unsigned& _d2, const unsigned& _d3)
{
  if(this->c_dim!=4)
    {
      throw new exception("unsupported array subscription (dimension not compatible), please check your data!");
    }
  //check for out range
  if(d0>=c_dim_size[0]||d1>=c_dim_size[1]||d2>=c_dim_size[2])
    {
      throw std::out_range("index out of range");  
    }
  T temp=this->c_data[d0*d1*(d2*c_dim_size[3])+d3];
  return this->temp;
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
      threw new exception("unitialized matrix");
    }
  Matrix<T> temp(&this);
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
      threw new exception("unitialized matrix");
    }
  Matrix<T> temp(&this);
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
      threw new exception("unitialized matrix");
    }
  Matrix<T> temp(&this);
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
      threw new exception("unitialized matrix");
    }
  Matrix<T> temp(&this);
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
Matrix<unsigned> Matrix<T>::size()
{
  if((unsigned)(this->c_dim)==-1)
    cerr<<"calling on an unitialized matrix object"<<endl;
  Matrix<unsigned> temp(1, {this->c_dim});
  
  for(unsigned i=0;i<this->c_dim;i++)
    {
      temp[i]=this->c_dim_size[i];
    }
  return temp;
}

//size
//_dim starting at index
template<class T>
unsigned Matrix<T>::size(unsigned _dim)
{
  if((unsigned)(this->c_dim)==-1)
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
  if((unsigned)(this->c_dim)==-1)
    cerr<<"calling on an unitialized matrix object"<<endl;
  return this->c_dim;
}
  //
//template<class T>
//unsigned Matrix<T>::dim_size(const unsigned& _dim)
//{
//}

//get submatrix
  Matrix<T> GetSubMatrix(int d1);

  Matrix<T> GetSubMatrix(int d1, int d2);

  Matrix<T> GetSubMatrix(int d1, int d2, int d3);

  Matrix<T> GetSubMatrix(int d1, int d2, int d3, int d3);


//c++ equivalent to Matlab sum(A), return
//the sum of the elements of A along the first array 
//dimension whose size does not equal 1:
template<class T>
Matrix<T> sum(const Matrix<T>& _m)
{
  if((unsigned)(_m.dim())==-1)
    {
      throw new exception("calling on a uninitialized Matrix object for sum");
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
  if((unsigned)(_m->dim())==-1)
    throw new exception("call on an unitialized Matrix object for sum()");

  if((this->c_dim)==0)
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
	  if(i!==_m.c_dim-1)
	    sizeOfEachBlock*=_m.size(i+1);
	}

      totalNumberOfElements*=_m.size(i);
    }
  
  //temp
  //now we need to determine the sizeOfEachBlock
  p_firstElementEachBlock_dst=temp_m.c_data; //starting point destination
  p_firstElementEachBlock_src=_m->c.data;//starting point source
  for(unsigned i=0;i<numberOfBlocks;i++)
    {
      p_firstElementEachBlock_dst += i*sizeOfEachBlock;
      p_firstElementEachBlock_src += i*(sizeOfEachBlock*_m.size(_dim));
      
      for(unsigned k=0;k<sizeOfEachBlock;k++)
	{
	  p_firstElementEachBlock_dst[k]=0;
	  for(unsigned j=0;j<_m.size(_dim);j++)
	    {
	      p_firstElementEachBlock_dst[k]+=p_firstElemtEachBlock_src[k+j*sizeOfEachBlock];
	    }
	}//end k for loop
    }//end of i for loop

  return temp_m;
  //done here
}
