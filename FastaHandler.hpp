#ifndef FASTAHANDLER_HPP
#define FASTAHANDLER_HPP
#include <string>
#include <vector>
#include <iostream>

#include "SequenceString.hpp"

//here, we won't define a class, instead we simply write up a few functions to
//read, write, etc
using namespace std;

//return total number of sequeces read in.
//toUpper =true : convert the character to upper case. false: keep the original letter
unsigned int ReadFasta(const string& _fname, vector<SequenceString>& _seqStrVec, bool toUpper=false);


void WriteFasta(const string& _fname, vector<SequenceString>& _seqStrVec, const unsigned int _width=50, ios_base::openmode mode=ios_base::out);

void WriteTextFile(const string& _fname, vector<int>& _seqStrVec, const char& c='\t', const unsigned int _width=1, ios_base::openmode mode=ios_base::out);




#endif
