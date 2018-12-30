#ifndef GZTOOLS_HPP
#define GZTOOLS_HPP

//this is a C++ accessory library to take care of gzip file input and output
//we rely on the zlib to do uncompression

#include "zlib.h"
#include <string>
#include <vector> 

using namespace std;

//in this library, do we need to do objects??? not sure
//not now!!! at this point (12/9/2018) we simply define a 
//helper functions to take care of reading from gzipped files
//just reading. not doing writtings for now.

//one thing to note is that gzFile is a pointer 

//to read one line from an opened gzFile stream 
gzFile getline(gzFile _if, string& _line, const char& _delimit='\n');

//In there, just for fun, we will create two different ways of get line
//but will only be called (used) by getline for better typedef, these 
//two will be moved to cpp, so not be called directly by other files
bool  getline_A(gzFile _if, string& _line); //, const char& _delimit='\n');
gzFile getline_B(gzFile _if, string& _line, const char& _delimit='\n');

//for the following two, we simply calling the correct gzopen and gzclose
//function, we have them here to make holes for the future usage if we
//we want to do more

gzFile gzOpen(const string& _fname, const string& _mode="rb");

//return null upon return, in case the gzclose were called twice on
//the same file handler.
void gzClose(gzFile _if);

//assuming the file has been opened 
bool getline_B(FILE* _f, string& l);

FILE* gzOpen_B(const string& _fname, const string & _mode="rb");

void gzClose_B(FILE* _f);

int init_gzip_stream(bool full=false/*FILE* file,char* out*/);
bool inflate_gzip( z_stream& strm, const size_t& bytes_read, const size_t& bytes_avail, int& code);


void debugZstream(const string& s="");


//a helper function to find a char in a char buffer, note here we don't assume the char buffer is 
//a null-terminated char c string. so we have to require a size of char buffer.
//therefore we stop either null or size whichever comes first.
//so this is why we don't use the build in strstr, since it assume null terminator.

char* strcnstr(char* source, size_t size, const char& c);
#endif

