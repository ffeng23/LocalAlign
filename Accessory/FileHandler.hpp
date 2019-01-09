#ifndef FILEHANDLER_HPP
#define FILEHANDLER_HPP
//in here we define some general functions to take care of files
#include <vector>
#include <string>
#include <iostream>
#include <sys/stat.h>
#include "SequenceString.hpp"


using namespace std;

enum FileType {GZ, FASTA, FASTQ, TXT, GZ_FASTA, GZ_FASTQ, GZ_TXT, DIR, UNKNOWN};

//get file type, note FileType is user-defined enum.
FileType getFileType(const string& fname);

bool is_file(const char* path) ;

bool is_dir(const char* path) ;

bool exist(const char* path);

//a file handler to read file into a fasta vector vector
//we try to detect the following thing in this function 
// 1, file exist?
// 2, is it a file or directory
// 3, is it a fasta, fastq or gziped fasta, gziped fastq. No other type supported so far
// 
// We will return a vector holding the squenences of the file (SquenceStrings)
// we will return the total number of sequences read in. If none or error, we will return 
// string::npos.
size_t readFile2SeqStrVector(string _fname, vector<SequenceString>& _vec);

#endif