#ifndef LOAD_DATA_HPP
#define LOAD_DATA_HPP

#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "../string_ext.hpp"
#include "../SequenceString.hpp"

using namespace std;
unsigned LoadData(const string& _fileName, vector<string>& _header, vector<SequenceString>& _seq, vector<unsigned>& _count, 
		  const char& _commentChar='#', const char& _delimiter='\t'
		  ,const bool& _data_header_lin=true);


unsigned ParseField(vector<string>& _header, const string& _field);

#endif
