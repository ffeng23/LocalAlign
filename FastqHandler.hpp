#ifndef FASTQHANDLER_HPP
#define FASTQHANDLER_HPP
#include <string>
#include <vector>
#include <iostream>

#include "SequenceString.hpp"
#include "Accessory/GzTools.hpp"
#include "FASTQ.hpp"

//here, we won't define a class, instead we simply write up a few functions to
//read, write, etc
using namespace std;

//return total number of sequeces read in.
//toUpper =true : convert the character to upper case. false: keep the original letter
unsigned int ReadFastq(const string& _fname, vector<Fastq>& _seqStrVec, bool toUpper=false);


void WriteFastq(const string& _fname, vector<Fastq>& _seqStrVec,  ios_base::openmode mode=ios_base::out);


#endif
