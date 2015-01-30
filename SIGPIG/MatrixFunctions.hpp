#ifndef MATRIX_FUNCTIONS_HPP
#define MATRIX_FUNCTIONS_HPP

#include <vector>
#include <iostream>
#include <string>
#include <cstring>

//this is the module to define some basic and necessary functions for matrix/array manipulations
using namespace std;

double max_mf(const vector<double>& _m);
double max_mf(const double*  _m, const unsigned& _len);
unsigned max_mf(const vector<unsigned>& _m);
unsigned max_mf(const unsigned* _m, const unsigned& _len);


double min_mf(const vector<double>& _m);
double min_mf(const double*  _m, const unsigned& _len);
unsigned min_mf(const vector<unsigned>& _m);
unsigned min_mf(const unsigned* _m, const unsigned& _len);

//sort 
//void sort_mf(const double* _m, const unsigned& _len)
//Quick sort, _a is the array to be sorted,
//_first, _last are the indices of the elements to be sorted in the array
//_b is the array that used to break the ties of the elements in the _a array
//_index is the output array of the indice of the sorted elements
template <class T> void QuickSort(T* _a, const unsigned& _first, const unsigned& _last, unsigned* _index=NULL, unsigned* _b=NULL);
 
/**
 * Find and return the index of pivot element.
 * @param a - The array.
 * @param first - The start of the sequence.
 * @param last - The end of the sequence.
 * @return - the pivot element
 */
//here we use a new way to choose pivot. a median between [first], [middle] and [last]
template <class T> unsigned Pivot(T* _a, const unsigned& _first, const unsigned& _last, unsigned* _index=NULL, unsigned* _b=NULL);

template <class T> void Swap(T& _a, T& _b);
//void swapNoTemp(int& a, int& b);
void Print(const double* _array, const int& _N, const unsigned* _index);
// double* sorted, unsigned * sorted_index);

template <class T> unsigned GetMedianIndex(const T* m, const unsigned& a, const unsigned& b, const unsigned& c, const unsigned* _index);

void Print(const unsigned* _array, const int& _N, const unsigned* _index);

void Reverse( unsigned* _array, const unsigned& _N);
void Reverse(double* _array, const unsigned& _N);
void Reverse(vector<unsigned> _v);
void Reverse(vector<double> _v);


bool CopyElements
  ( unsigned** _source, const unsigned& _s_size1, const unsigned& _s_size2, 
   unsigned** _target, const unsigned& _t_size1, const unsigned& _t_size2,
   const unsigned* _indexOfElementToCopy, const unsigned& _i_size);

bool CopyElements
  (const unsigned* _source, const unsigned& _s_size, 
   unsigned* _target, const unsigned& _t_size,
   const unsigned* _indexOfElementToCopy, const unsigned& _i_size);

/* Find unique values of the input array.
 * it requires the caller to initialize the output array of the size 
 * identical to the input array. the return output array will 
 * be smaller or equal size specified by _oSize
 *-------------
 * input:
 *    _in, input array with size of _iSize
 *    _out, output array with size _iSize but only the
 *         first _oSize are used. Again the caller has
 *         to initialize the array.
 *    _out_index, similar to _out, but holding the index
 *         of element in the original input array. it has
 *         of _iSize but again only the first _oSize is 
 *         occupied.
 */
void Unique(const unsigned* _in, const unsigned& _iSize, 
	    /*output*/ unsigned* _out, unsigned* _out_index, 
	    unsigned& _oSize); 

#endif
