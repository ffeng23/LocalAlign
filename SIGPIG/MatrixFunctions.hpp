#ifndef MATRIX_FUNCTIONS_HPP
#define MATRIX_FUNCTIONS_HPP

#include <vector>
#include <iostream>
#include <string>

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
void QuickSort(double* _a, const unsigned& _first, const unsigned& _last, unsigned* _index=NULL, unsigned* _b=NULL);
unsigned Pivot(double* _a, const unsigned& _first, const unsigned& _last, unsigned* _index=NULL, unsigned* _b=NULL);
void Swap(double& _a, double & _b);
//void swapNoTemp(int& a, int& b);
void Print(const double* _array, const int& _N);
// double* sorted, unsigned * sorted_index);
unsigned GetMedianIndex(const double* m, const unsigned& a, const unsigned& b, const unsigned& c);

void Print(const unsigned* _array, const int& _N);

void Reverse( unsigned* _array, const unsigned& _N);
void Reverse(double* _array, const unsigned& _N);
void Reverse(vector<unsigned> _v);
void Reverse(vector<double> _v);


#endif
