#ifndef SEQUENCEHANDLER_HPP
#define SEQUENCEHANDLER_HPP

#include <vector>
#include <string>
#include "SequenceString.hpp"
#include "score.hpp"

using namespace std;

//this translation unit is to define some function to take care of sequence mapping
//in this function, we combine the adaptor+mid+primer
void ProcessAdaptorSequences(const string& _adaptorFile, const string& _barcodeFile, const string& _forwardFile, 
			     const string& _reverseFile, vector<SequenceString>& _vecForward, 
			     vector<SequenceString>& _vecReverse );
void MappingAdaptors(vector<SequenceString>& _vecForward, vector<SequenceString>& _vecReverse, vector<SequenceString>& _vecSeq, 
		     ScoreMatrix* _sm, const double& _gapOpen, const double& _gapExtension, const unsigned int& _trim,
		     const double& _mismatchRateThreshold, const unsigned _minmumOverlapLength, const unsigned int& _offsetForward, const unsigned int& _offsetReverse, 
		     const string& _mapBoth_fname, const string& _mapForward_fname,
		     const string& _mapReverse_fname, const string& _mapNone_fname); 

SequenceString ReverseComplement(SequenceString& seq);

//compare two strings character by character, return # of chars that are different. if two strings are of different size, 
//the ones longer are also counted as different chars
//for example: "abc" vs. "acc" return 1;
//"abc" vs "acb" return 2
//"abc" vs "ab" reurn 1
//"abc" vs "cab" return 3

unsigned int CompareStrings(const string& str1, const string& str2);

#endif

