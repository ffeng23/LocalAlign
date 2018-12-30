#include "FileHandler.hpp"
#include "string_ext.hpp"
#include <string>
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


