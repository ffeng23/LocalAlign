#ifndef SEQUENCEHANDLER_BARCODE_HPP
#define SEQUENCEHANDLER_BARCODE_HPP

#include <vector>
#include <string>
#include "SequenceString.hpp"
#include "score.hpp"

using namespace std;

enum mapType{ FivePrime, ThreePrime};
//IN this function, we take care to do the isotype identification

void MappingBarcodes(vector<SequenceString>& _vecSeq1, /*this is the sequence data r1*/
		    vector<SequenceString>& _vecSeq2, /*this is the sequence data r2*/
		     vector<SequenceString>& _vecIndex1, /*this is the sequence index data r1*/
		     vector<SequenceString>& _vecIndex2, /*this is the sequence index data r2*/
		     vector<SequenceString>& _vecBarSeq1, /*this is the barcode sequences R1*/
		     vector<SequenceString>& _vecBarSeq2, /*this is the barcode sequences R2*/
		     const unsigned& misMatchNum, /*total number of mismatches allowed*/
		     const bool& pairEnd, /*indicating whether to do pairEnd*/
		     const bool& dualIndex, /*indicating whether to do dual indexes*/
		     const bool& indexFromSeq, /*indicating whether to do index from header of sequenences*/
		     const mapType& type, /*indicating whether it is 5'prime or 3' prime*/
		     const bool& demux,/*indicating whether to do dumux as out, otherwise only output stats*/
		     const string& _out_fname//, const string& _map_fname
		      ); 

unsigned int ReadIndexFromSequenceName(vector<SequenceString>& _vecSeq,  /*this is the sequence data r1*/
			      vector<SequenceString>& _vecIndex1, /*this is the sequence index data r1*/
			      vector<SequenceString>& _vecIndex2, /*this is the sequence index data r2*/
			      const unsigned& lenOfBarcode,
			      const bool& dualIndex /*indicating whether to do dual indexes*/
				       );

#endif


