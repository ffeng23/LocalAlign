#include <Matrix.hpp>

//constructor
template<class T>
Matrix<T>::Matrix():c_dim(0),c_dim_size(NULL),c_data(NULL)
{
  //default constructor
}

//data input now is one data.
//we will restructure it
template<class T>
Matrix<T>::Matrix(unsigned _dim, unsigned []_dim_size, T[] _data):
  c_dim(_dim)
{
  c_dim_size=new unsigned[c_dim];
  unsigned totalNumData=1;
  //need to copy over the data and size
  for(unsigned i=0;i<c_dim;i++)
    {
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
//vector
template<class T>
T Matrix<T>::operator [](const unsigned& _d0)
{
  if(this->c_dim!=1)
    {
      throw new exception("unsupported array subscription (dimension not compatible), please check your data!")
    }
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
  T temp=this->c_data[d0*(d1*c_dim_size[2])+d2];
  return this->temp;
}
  
//4d
template<class T>
T& Matrix<T>::operator [](const unsigned& _d0, const unsigned& _d1, const unsigned& _d2, const unsigned& _d3)
{
if(this->c_dim!=4)
    {
      throw new exception("unsupported array subscription (dimension not compatible), please check your data!")
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
  Matrix<unsigned> temp(1, {this->c_dim});
  
  for(unsigned i=0;i<this->c_dim;i++)
    {
      temp[i]=this->c_dim_size[i];
    }
  return temp;
}
  
//
template<class T>
unsigned Matrix<T>::dim()
{
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
Matrix<T> sum(const Matrix<T>& _m)
{


}

//c++ equivalent to Matlab sum(A,dim), return
//sums the elements of A along dimension dim. 
//The dim input is a positive integer scalar.
Matrix<T> sum(const Matrix<T>& _m, const unsigned& _dim)
{


}
