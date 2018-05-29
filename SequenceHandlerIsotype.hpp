#ifndef SEQUENCEHANDLER_ISOTYPE_HPP
#define SEQUENCEHANDLER_ISOTYPE_HPP

#include <vector>
#include <string>
#include "SequenceString.hpp"
#include "score.hpp"

using namespace std;

enum mapType{ 5prime, 3prime};
//IN this function, we take care to do the isotype identification
void MappingIsotypes(vector<SequenceString>& _vecSeq, /*this is the sequence data that we want to find isotypes*/
		     vector<SequenceString>& _vecIsotype, /*this is the isotype sequences*/
		     const mapType& type, /*indicating whether it is 5'prime or 3' prime*/
		     ScoreMatrix* _sm, const double& _gapOpen, const double& _gapExtension,
		     const double& _mismatchRateThreshold, const unsigned _minmumOverlapLength,
		     const unsigned int& _offset, 
		     const string& _map_fname,
		       const string& _unmap_fname//,
		     //const string& _mapReverse_fname, const string& _mapNone_fname,
		     //const string& _mapCrossOver_fname, const string& _mapBreakOutsideCon_fname
		      ); 


/*void SetUpTrimFlag(const bool& _f);*/
/*void SetUpByIsotypeOutputFlag(const bool& _f);*/
#endif

