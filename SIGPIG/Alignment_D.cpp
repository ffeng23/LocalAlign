#include <cstring>
#include <sstream>
#include "../SequenceString.hpp"
#include "genomicSegments.hpp"
#include "Alignment.hpp"
#include "MatrixFunctions.hpp"
#include "../string_ext.hpp"
#include "AlignmentSettings.hpp"

#include "Alignment_D.hpp"


unsigned Alignment_D::allele_order []={ 0,
					  1, 2, 3, 4, 5, 6, 7, 8, 9,10,
					  11,12,13,14,15,16,17,18,19,20,
					  21,22,23,24,25,26,27,28,29,30,
					  31,32,33
					  };

unsigned Alignment_D:n_D_alleles=AlignmentSettings::n_D_alleles;


Alignment
