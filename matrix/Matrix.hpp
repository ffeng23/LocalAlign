#ifndef MATRIX_HPP
#define MATRIX_HPP

//user defined matrix class for matrix manipulation
template <class T>
class Matrix
{
public:

  //constructor
  Matrix();

  //
  Matrix(unsigned _dim, unsigned []_dim_size, unsigned[] _data);

  //destructor
  ~Matrix();

  //copy constructor
  Matrix(const Matrix<T>& _m);

  //assignment operator
  

  //


protected:
  unsigned c_dim;
  unsigned* c_dim_size;

  T* c_data; //the real data container

};

#endif
