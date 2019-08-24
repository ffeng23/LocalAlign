#ifndef FASTQHANDLER_HPP
#define FASTQHANDLER_HPP
#include <string>
#include <vector>
#include <iostream>
#include <string.h>
#include "SequenceString.hpp"
#include "GzTools.hpp"
#include "FASTQ.hpp"

//here, we won't define a class, instead we simply write up a few functions to
//read, write, etc
using namespace std;

//return total number of sequeces read in.
//toUpper =true : convert the character to upper case. false: keep the original letter
//
//the file could be gz'ed or regular fastq. then the sequence data are read into a vector of 
//fastq objects. 
//return string::npos upon error 
size_t ReadFastq(const string& _fname, vector<Fastq>& _seqStrVec, bool toUpper=false);

size_t ReadFastq(const string& _fname, vector<SequenceString>& _seqStrVec, vector<string>& _vecQ, bool toUpper=false);


void WriteFastq(const string& _fname, vector<Fastq>& _seqStrVec,  ios_base::openmode mode=ios_base::out);
//updated 8/21/2019, add this to write the fastq based on the separated sstring and quality.
void WriteFastq(const string& _fname, vector<SequenceString>& _seqStrVec, 
	vector<string>& _qVec,
	ios_base::openmode mode=ios_base::out);


#endif
