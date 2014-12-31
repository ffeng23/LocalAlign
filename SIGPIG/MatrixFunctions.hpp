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




#endif
