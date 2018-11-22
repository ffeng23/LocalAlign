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


void WriteFasta(const string& _fname, vector<SequenceString>& _seqStrVec, const unsigned int& _width=50, ios_base::openmode mode=ios_base::out);

//writing a text file, mainly for one vector.
//input: _fname, file name
//       _seqStrVec, vector of numbers to be written
//       _c, a char to delimite the columns
//       _width, how many columns to write
//       mode, to open the files, truc (new) or append
void WriteTextFile(const string& _fname, vector<unsigned int>& _seqStrVec, const char& c='\t', const unsigned int& _width=1, ios_base::openmode mode=ios_base::out);

//Writing a text table file, with header or not
//input: _fname, file name
//       _seqStrVec, vector of vector of  numbers to be written
//       _c, a char to delimite the columns
//       _header, whether to write the header
//       mode, to open the files, truc (new) or append
//       _headerStr, the header names to be written.
void WriteTextTableFile(const string& _fname, vector<vector<double> >& _seqStrVec, const char& c='\t', const bool& _header=true, ios_base::openmode mode=ios_base::out,vector<string> _headerStr=vector<string>());


#endif
