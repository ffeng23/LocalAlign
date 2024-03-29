#ifndef SEQUENCEHANDLER_CONSTANT_HPP
#define SEQUENCEHANDLER_CONSTANT_HPP

#include <vector>
#include <string>
#include "Accessory/SequenceString.hpp"
#include "score.hpp"

using namespace std;

//this translation unit is to define some function to take care of sequence mapping
//in this function, we combine the adaptor+mid+primer
/*void ProcessAdaptorSequences(const string& _adaptorFile, const string& _barcodeFile, const string& _forwardFile, 
			     const string& _reverseFile, vector<SequenceString>& _vecForward, 
			     vector<SequenceString>& _vecReverse );*/

void MappingConstants(vector<SequenceString>& _vecForward, vector<SequenceString>& _vecReverse, vector<SequenceString>& _vecSeq, 
		     ScoreMatrix* _sm, const double& _gapOpen, const double& _gapExtension,
		      const double& _mismatchRateThreshold, const unsigned _minmumOverlapLength, const unsigned int& _offsetForward, const unsigned int& _offsetReverse, 
		     const string& _mapBoth_fname, const string& _mapForward_fname,
		      const string& _mapReverse_fname, const string& _mapNone_fname,
		      const string& _mapCrossOver_fname, const string& _mapBreakOutsideCon_fname
		      ); 


/*void SetUpTrimFlag(const bool& _f);*/
/*void SetUpByIsotypeOutputFlag(const bool& _f);*/
#endif

