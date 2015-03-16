#ifndef CONCENSUS_HPP
#define CONCENSUS_HPP

#include <vector>
#include "../SequenceString.hpp"
#include "../FastaHandler.hpp"

//start defining the functions to run the Concensus sequences

//function

//run the concensus for each individual file, which will give out one concensus sequence
//_numOfSeq is the number of seq in the file to generate the concensus
SequenceString GetConcensus(const string& _fileName, unsigned& _numOfSeq);

//go through a bunch of input file and generate one sequence from each file,
//we also want to remember how many sequences for each concensus sequence.
//the outside caller has to initialize the ss array.
void GenerateConcensus(const string* _ifileNames, const unsigned& _numOfFiles
		       /*output*/ SequenceString* _ss, unsigned* countOfConcensus,
		       );

//based one the concensus sequences, write up the input sequence file for our analys
//it has specific format, like 
//header
//sequence count
//see the ../sample.data for the detail format.
//this format is adopted from the Matlab code.
void GenerateSequenceFile(const SequenceString* _ss, const unsigned* countOfConcensus,
			  const unsigned& _numOfSeq, const string& _ofname);

#endif
