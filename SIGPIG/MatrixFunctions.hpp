#ifndef MATRIX_FUNCTIONS_HPP
#define MATRIX_FUNCTIONS_HPP

#include <vector>
#include <iostream>
#include <string>
//this is the module to define some basic and necessary functions for matrix/array manipulations

double max(vector<double> _m);
double max(const double*  _m, const unsigned& _len);
unsigned max(vector<unsigned> _m);
unsigned max(const unsigned* _m, const unsigned& _len);


double min(vector<double> _m);
double min(const double*  _m, const unsigned& _len);
unsigned min(vector<unsigned> _m);
unsigned min(const unsigned* _m, const unsigned& _len);




#endif
