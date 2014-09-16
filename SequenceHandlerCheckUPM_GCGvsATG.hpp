#ifndef SEQUENCEHANDLERCHECKUPM_GCGVSATG_HPP
#define SEQUENCEHANDLERCHECKUPM_GCGVSATG_HPP

#include <vector>
#include <string>
#include "SequenceString.hpp"
#include "score.hpp"

using namespace std;


void CheckingUPMGCGvsATC(vector<SequenceString>& _vecReverse, vector<SequenceString>& _vecSeq, 
		     ScoreMatrix* _sm, const double& _gapOpen, const double& _gapExtension,
		      const double& _mismatchRateThreshold, const unsigned _minmumOverlapLength, 
		     
			 const string& _mapReverse_fname, const string& _mapNone_fname
		      ); 

#endif

