#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <string>


//user defined matrix class for matrix manipulation
//so far we only support up to 4 Dimensions,
//NOTE: a dimension of zero, means a scalar !!!
//    1D is vector
//

//ToDo: need to do one bulk load of data into the Matrix
template <class T>
class Matrix
{
public:

  //constructor
  Matrix();

  //data input now is one data.
  //we will restructure it
  Matrix(const unsigned& _dim, unsigned _dim_size[],  T _data[]);

  //
  Matrix(const unsigned& _dim, const unsigned* _dim_size, const T* _data=NULL);

  //destructor
  ~Matrix();

  //copy constructor
  Matrix(const Matrix<T>& _m);

  //assignment operator
  Matrix<T>& operator = (const Matrix<T>& _m);

  //subscript
  //scalar-like matrix
  T& operator ()();
  //vector
  T& operator ()(const unsigned& _d1);
  
  //2d
  T& operator ()(const unsigned& _d1, const unsigned& _d2);
  
  //3d
  T& operator ()(const unsigned& _d1, const unsigned& _d2, const unsigned& _d3);
  
  //4d
  T& operator ()(const unsigned& _d1, const unsigned& _d2, const unsigned& _d3, const unsigned& _d4);
  
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
  Matrix<unsigned> size() const;
  
  //return unsigned dim_size of specific dimension
  unsigned size(const unsigned& _dim) const;

  //
  unsigned dim() const ;

  //get submatrix
  Matrix<T> SubMatrix(const unsigned& _n, const int _dim_pos[]) const;

  //get submatrix
  //Matrix<T> SubMatrix(const unsigned& _n, const int (&_dim_pos)[]) const;

  std::string toString() const;

  //this is  a function only used to initialize/fill in data for an empty matrix
  //
  void initialize(const unsigned& _dim, const unsigned _dim_size[], const T _data[]);
  /*
  Matrix<T> GetSubMatrix(int d1, int d2);

  Matrix<T> GetSubMatrix(int d1, int d2, int d3);

  Matrix<T> GetSubMatrix(int d1, int d2, int d3, int d3);
  */
  template<class U>
  friend Matrix<U> sum(const Matrix<U>& _m, const unsigned& _dim);
  //
protected:
  unsigned c_dim;
  unsigned* c_dim_size;

  T* c_data; //the real data container, is allocated in one d

};
//c++ equivalent to Matlab sum(A), return
//the sum of the elements of A along the first array 
//dimension whose size does not equal 1:
template<class T>
Matrix<T> sum(const Matrix<T>& _m);


//c++ equivalent to Matlab sum(A,dim), return
//sums the elements of A along dimension dim. 
//The dim input is a positive integer scalar.
template<class T>
Matrix<T> sum(const Matrix<T>& _m, const unsigned& _dim);

#endif
