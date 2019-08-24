#include "FileHandler.hpp"
#include "string_ext.hpp"
#include "FASTQ.hpp"
#include "FastaHandler.hpp"
#include "FastqHandler.hpp"
#include <string>
#include <string.h>


bool is_file(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISREG(buf.st_mode);
}

bool is_dir(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISDIR(buf.st_mode);
}

bool exist(const char* path) {
    struct stat buf;
	
    int x= stat(path, &buf);
	cout<<"inside the exist:x="<<x<<endl;
    return x==0;
}

//get file type, note FileType is user-defined enum.
FileType getFileType(const string& fname)
{
	const char* f=fname.c_str();
	//
	if(is_dir(f))
	{
		return DIR;
	}
	//
	if(!is_file(f))
	{
		return UNKNOWN;
	}
	//here it must be a regular file 
	size_t dot_pos;
	string suffix;
	//bool tryAgain=false;
	size_t size=fname.size()-1;
	bool gz_flag=false;
	
	//cout<<"*****deteting....."<<endl;
	do{
		//cout<<"\ttrying......."<<endl;
		//tryAgain=false;
		//Check for the suffix.
		dot_pos=fname.find_last_of('.',size);
		//cout<<"dot-Pos:"<<dot_pos<<";string:npos"<<string::npos<<endl;
		if(dot_pos==string::npos)
		{
			//can not find it, then 
			if(gz_flag)
				return GZ;
			return UNKNOWN;
		}
		
		
		//suffix without '.'
		suffix=fname.substr(dot_pos+1, size-dot_pos);
		
		suffix=to_lower_str(suffix);
		//cout<<"suffix:"<<suffix<<endl;
		
		size=dot_pos-1;
		//cout<<"size:"<<size<<endl;
		if(suffix.compare("txt")==0)
		{
			//cout<<"here"<<endl;
			if(gz_flag)
				return GZ_TXT;
			return TXT;
		}
		if(suffix.compare("fasta")==0||suffix.compare("fa")==0||suffix.compare("fas")==0)
		{
			if(gz_flag)
				return GZ_FASTA;
			return FASTA;
		}
		if(suffix.compare("fq")==0||suffix.compare("fastq")==0)
		{
			if(gz_flag)
				return GZ_FASTQ;
			return FASTQ;
		}
		
		if(suffix.compare("gz")==0||suffix.compare("gzip")==0)
		{
			
			if(gz_flag)
			{//multi-gz'ed
				cout<<"WARNING: In detecting file type, multiple compression found!!"<<endl;
				return UNKNOWN;
			}
			gz_flag=true;
			//tryAgain=true;
		}
	}	
	while (gz_flag);
	
	return UNKNOWN;
}


//a file handler to read file into a fasta vector vector
//we try to detect the following thing in this function 
// 1, file exist?
// 2, is it a file or directory
// 3, is it a fasta, fastq or gziped fasta, gziped fastq. No other type supported so far
// 
// We will return a vector holding the squenences of the file (SquenceStrings)
// we will return the total number of sequences read in. If none or error, we will return 
// string::npos.
size_t readFile2SeqStrVector(string _fname, vector<SequenceString>& _vec, vector<string>* _vec_Q)
{
	//check for files
	if(!exist(_fname.c_str()))
	{
		cout<<"ERROR: the specified file ("<<_fname<<") does not exist!!"<<endl;
		return string::npos;
	}
	//check for files
	if(is_dir(_fname.c_str()))
	{
		cout<<"ERROR: the specified file ("<<_fname<<") is a directory!!"<<endl;
		return string::npos;
	}
	//check file type 
	FileType ft=getFileType(_fname);
	size_t numReads;
	switch(ft)
	{
		case FASTA:
	
		case GZ_FASTA:
			numReads=ReadFasta(_fname, _vec);
			break;
			
		case FASTQ:
		case GZ_FASTQ:
		{
			/*if(_vec_Q==NULL)
			{
				cout<<"ERROR: Quality vector(s) for fastq data has not been initalizated correctly. check !!"<<endl;
				exit(-1);
			}*/
			vector<Fastq> _vec_fq;
			numReads=ReadFastq(_fname, _vec_fq);
			if(numReads==string::npos)
				return numReads;
			//rewrite it into fasta vector
			for(unsigned i=0;i<numReads;i++)
			{
				_vec.push_back(_vec_fq.at(i).GetSequenceString());
				_vec_Q->push_back(_vec_fq.at(i).GetQualityString());
			}
			numReads=_vec.size();
			break;
		}	
		case GZ:
		case GZ_TXT:
		case TXT:
		case UNKNOWN:
		case DIR:
		default:
			cout<<"ERROR: unsupported file format for this action"<<endl;
			numReads= string::npos;
			break;
	}
	
	return _vec.size();
}	

