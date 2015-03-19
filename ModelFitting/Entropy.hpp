#ifndef ENTROPY_HPP
#define ENTROPY_HPP

#include <cmath>
#include <iostream>
#include "../matrix/Matrix.hpp"
#include "VDJ_cuts_insertion_dinuc_ntbias_model.hpp"

double Entropy( const Matrix<double>& _p, const double& _base);

//int feng(int x);
//void GetAssignmentEntropy(VDJ_cuts_insertion_dinuc_ntbias_model& _model);

#endif
