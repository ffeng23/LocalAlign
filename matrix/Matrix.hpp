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

  //clear the matrix data, make it unitialized, but keep the object
  void clear();
  
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

  //the input _m is bool matrix and we return a vector containning the elements
  //that are true in _m entries. vector is following matlab style
  Matrix<T> GetElements(const Matrix<bool>& _m)const;

  //this allows specify the indices by an array or the array has to be of the same
  //size as the dim of the matrix, so _dim input is the size of the indices array
  //
  T& operator()(const unsigned* _indices, const unsigned& _dim)const;
  
  //NOTE: only support 4d or below

  //other operator, so all .dot operation
  //add
  Matrix<T> operator + (const T& _t) const;

  //other operator, so all .dot operation
  //add
  Matrix<T> operator + (const Matrix<T>& _t) const;

  //subtract
  Matrix<T> operator - (const T& _t) const;

  //multiplication
  Matrix<T> operator * (const T& _t) const;

  //division
  Matrix<T> operator / (const T& _t) const;

  Matrix<T> operator /(const Matrix<T>& _m)const;

  //comparision
  Matrix<bool> operator >(const T& _t)const;
  Matrix<bool> operator <(const T& _t)const;
  Matrix<bool> operator <=(const T& _t)const;
  Matrix<bool> operator >=(const T& _t)const;

  //bit wise operation
  Matrix<T> operator & (const Matrix<T>& _m)const;

   //size 
  Matrix<unsigned> size() const;
  
  //return unsigned dim_size of specific dimension
  unsigned size(const unsigned& _dim) const;

  //return total number of elements
  unsigned nTotal()const;
  //
  unsigned dim() const ;

  //get submatrix
  // _n is the size of following _dim_pos array. number of element in that array
  // _dim_pos[] is the arraying holding the positions/indices of the sub matrix
  //       it could be 0~some valid indecies or -1. -1 mean in that dimension take all
  Matrix<T> SubMatrix(const unsigned& _n, const int _dim_pos[]) const;

  //Set submatrix
  //this is the special case for set up the submatrix with one-d array.
  //we made this to make the data initalization easier
  //_d, is the first dimension position for the sub matrix to go. in this 
  //    case the target matrix to be set the submatrix must have dimension of 2
  //_size, the size of input _matrix, also equals to the size of 
  //     dimension 2 of target matrix
  //_matrix, array of input
  void SetSubMatrix(const unsigned& _d, const unsigned& _size, const T _matrix[]); 

  //this is generic set submatrix function
  //currently it only allows the first dimention to be specified for submatrix
  //so the target matrix and input matrix have only one dimension difference
  //for example, target matrix of dimension 3, then input dimension has to be 2
  //     we specified where in terms of first dimension position for which the 
  //     input dimesion can be set
  //input
  // _d, demsion 1 position to set the submatrix
  //_matrix, submatrix of input, need to have identical dimension
  //        to the dimension 2 and above of the target dimentsion
  void SetSubMatrix(const unsigned& _d, const Matrix<T>& _matrix);

  //this is generic set submatrix function, overloade method
  //it only allows the first and second dimention to be specified for submatrix
  //so the target matrix and input matrix have 2 dimension difference
  //for example, target matrix of dimension 3, then input dimension has to be 1
  //     we specified where in terms of first and second dimension position for which the 
  //     input dimesion can be set
  //input
  // _d1, dimension 1 position to set the submatrix
  // _d2, dimension 2 position to set the submatrix
  //_size, the size of input _matrix, also equals to the size of 
  //     dimension 3 of target matrix
  //_matrix, array of input
  
  void SetSubMatrix(const unsigned& _d1, const unsigned& _d2, 
		    const unsigned& _size, const T _matrix[]);



  //===>PENDING issues, need to add function to take care more dimension cases
  //    currently only support set sub at dimension one. need to take care more later


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

  //here, we only copy over the valid elements (is not a nan), otherwise
  //leave the old value unchanged
  //return false if the size and dimension are not compatible
  bool CopyValidElements( const Matrix<T>& _src);

  //to vector, arranged like what matlab does. it is (1,1), (2,1), (3,1)...(1,2),(2,2)
  //..(1,3),(2,3)...not like we I have in here inside c_data, eg. (1,1), (1,2),
  //(1,3)...(2,1),(2,2)....
  //
  Matrix<T> m2vec() const;

  //-------friend zone----------
  //dot divide
  
  template<class U>
  friend Matrix<U> sum(const Matrix<U>& _m, const unsigned& _dim);

  template<class U>
  friend U sum_all(const Matrix<U> & _m);
  //

  friend unsigned sum_all_bool(const Matrix<bool>& _m);
  
  template<class U>
  friend U max(const Matrix<U> & _m);

  template<class U>
  friend Matrix<U> max(const Matrix<U> & _m, const unsigned& _dim);

  //take the logarithm
  template<class U>
  friend Matrix<U> matrix_log(const Matrix<U>& _m);

  //do multiplication for 2 matrices
  template<class U>
  friend double matrix_multiply_1D(const Matrix<double>& _m1, const Matrix<U>& _m2);

  friend Matrix<double> matrix_multiply_by_scalar(Matrix<unsigned>& _m, const double& _s);
  //==================
protected:
  unsigned c_dim;
  unsigned* c_dim_size;

  T* c_data; //the real data container, is allocated in one d
  //in this internal data block, we arranged data in blocks. The block order is
  //different from Matlab, the last dimension got positioned together
  //.......then the first dimenstion is the outer most block
  //for eg, 2D, d0 is x(3) and d1 is y(2)
  // y0y1, y0y1, y0y1

  //it is like the order x0y0 x0y1 x1y0 x1y1 x2y0 x2y1.
  //to access it, we need xi yj
  //i*dim_y+j

};

//sepcialization, defining boolean matrix bitwise operation
template<>
Matrix<bool> Matrix<bool>::operator &(const Matrix<bool>& _m) const;

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

//specificall for summing bool matrix
unsigned sum_all_bool(const Matrix<bool>& _m);

//find the max element in the matrix
template<class T>
Matrix<T> max(const Matrix<T>& _m, const unsigned& _dim);

//find the max element in the matrix
template<class T>
T max(const Matrix<T>& _m);

//take the logarithm
template<class T>
Matrix<T> matrix_log(const Matrix<T>& _m);

//multiplication 1D equal size
template<class T>
double matrix_multiply_1D(const Matrix<double>& _m1, const Matrix<T>& _m2);

Matrix<double> matrix_multiply_by_scalar(Matrix<unsigned>& _m, const double& _s);

#endif
