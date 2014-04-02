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
unsigned int ReadFasta(const string& _fname, vector<SequenceString>& _seqStrVec);


void WriteFasta(const string& _fname, vector<SequenceString>& _seqStrVec, const unsigned int _width=50 );





#endif
