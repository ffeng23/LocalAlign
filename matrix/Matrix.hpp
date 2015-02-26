#ifndef MATRIX_HPP
#define MATRIX_HPP

//user defined matrix class for matrix manipulation
//so far we only support up to 4 Dimensions,
template <class T>
class Matrix
{
public:

  //constructor
  Matrix();

  //data input now is one data.
  //we will restructure it
  Matrix(unsigned _dim, unsigned []_dim_size, T[] _data);

  //destructor
  ~Matrix();

  //copy constructor
  Matrix(const Matrix<T>& _m);

  //assignment operator
  Matrix<T>& operator = (const Matrix<T>& _m);

  //subscript
  //vector
  T operator [](const unsigned& _d1);
  
  //2d
  T& operator [](const unsigned& _d1, const unsigned& _d2);
  
  //3d
  T& operator [](const unsigned& _d1, const unsigned& _d2, const unsigned& _d3);
  
  //4d
  T& operator [](const unsigned& _d1, const unsigned& _d2, const unsigned& _d3, const unsigned& _d4);
  
  //NOTE: only support 4d or below

  //other operator, so all .dot operation
  //add
  Matrix<T> operator + (const T& _t);

  //subtract
  Matrix<T> operator - (const T& _t);

  //multiplication
  Matrix<T> operator * (const T& _t);

  //division
  Matrix<T> operator / (const T& _t);
  
  //size 
  unsigned size();
  
  //
  unsigned dim();

  //
  unsigned dim_size(const unsigned& _dim);

  //get submatrix
  Matrix<T> GetSubMatrix(int d1);

  Matrix<T> GetSubMatrix(int d1, int d2);

  Matrix<T> GetSubMatrix(int d1, int d2, int d3);

  Matrix<T> GetSubMatrix(int d1, int d2, int d3, int d3);

  //
protected:
  unsigned c_dim;
  unsigned* c_dim_size;

  T* c_data; //the real data container, is allocated in one d

};
//c++ equivalent to Matlab sum(A), return
//the sum of the elements of A along the first array 
//dimension whose size does not equal 1:
Matrix<T> sum(const Matrix<T>& _m);

//c++ equivalent to Matlab sum(A,dim), return
//sums the elements of A along dimension dim. 
//The dim input is a positive integer scalar.
Matrix<T> sum(const Matrix<T>& _m, const unsigned& _dim);

#endif
