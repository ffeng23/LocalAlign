#ifndef FILEHANDLER_HPP
#define FILEHANDLER_HPP
//in here we define some general functions to take care of files

#include <string>
#include <iostream>
#include <sys/stat.h>

using namespace std;

enum FileType {GZ, FASTA, FASTQ, TXT, GZ_FASTA, GZ_FASTQ, GZ_TXT, DIR, UNKNOWN};

//get file type, note FileType is user-defined enum.
FileType getFileType(const string& fname);

bool is_file(const char* path) ;

bool is_dir(const char* path) ;

bool exist(const char* path);

#endif