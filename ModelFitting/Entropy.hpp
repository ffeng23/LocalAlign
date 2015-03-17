#ifndef ENTROPY_HPP
#define ENTROPY_HPP

#include <cmath>
#include <iostream>
#include "../matrix/Matrix.hpp"
#include "VDJ_cuts_insertions_dinuc_ntbias_model.hpp"

double Entropy( Matrix<double>& _p, const double& _base);

void CalculateAssignmentEntropy(VDJ_cuts_insertions_dinuc_ntbias_model& _model);

#endif
