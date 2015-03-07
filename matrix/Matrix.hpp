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
  Matrix(const unsigned& _dim, const unsigned _dim_size[],  const T& _data);

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

  //this allows specify the indices by an array or the array has to be of the same
  //size as the dim of the matrix, so _dim input is the size of the indices array
  //
  T& operator()(const unsigned* _indices, const unsigned& _dim);
  
  //NOTE: only support 4d or below

  //other operator, so all .dot operation
  //add
  Matrix<T> operator + (const T& _t);

  //other operator, so all .dot operation
  //add
  Matrix<T> operator + (const Matrix<T>& _t);

  //subtract
  Matrix<T> operator - (const T& _t);

  //multiplication
  Matrix<T> operator * (const T& _t);

  //division
  Matrix<T> operator / (const T& _t);

  //Matrix<T> operator /(const Matrix<T>& _m);
  
  //size 
  Matrix<unsigned> size() const;
  
  //return unsigned dim_size of specific dimension
  unsigned size(const unsigned& _dim) const;

  //return total number of elements
  unsigned nTotal()const;
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
  //initialize the object with identical elements
  void initialize(const unsigned& _dim, const unsigned _dim_size[], const T& _data);

  //normalize by dimension, in situ
  //carry out the division along the specified dimension
  //input matrix: _m (*this), must be larger or equal dimensions than _dim_vector
  //              _dim_vector must be a vector
  //             the size of _dim_vector must be identical to the size of _m at _dim dimension
  //return a bool indicating whether everything is all right.
  
  bool divide_by_dimension(const Matrix<T>& _dim_vector, const unsigned& _dim); 

  //this will use the internal one D array to retrieve elements _index
  //
  T Get1DArrayElement(const unsigned& _index)const;

  //now try to think the c_data as a 1D array, and then set its value
  void Set1DArrayElement(const unsigned& _index, const T& _t);
  
  //dot divide
  
  template<class U>
  friend Matrix<U> sum(const Matrix<U>& _m, const unsigned& _dim);

  template<class U>
  friend U sum_all(const Matrix<U> & _m);
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

//this is the one linearize the matrix data and sum all the entries together to 
//get the result
template<class T>
T sum_all(const Matrix<T>& _m);


#endif
