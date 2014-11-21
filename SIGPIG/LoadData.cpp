#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include "LoadData.hpp"


unsigned LoadData(const string& _fileName, vector<string>& _header, vector<SequenceString>& _seq, vector<unsigned>& _count, 
		  const char& _commentChar, const char& _delimiter, const bool& _data_header_line
		  )
{
  ifstream ifs_p(_fileName.c_str());
      
  if(!ifs_p.is_open())
    {
      cout<<">>>>>>ERROR:the input file \""<<_fileName<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  

  //so far, we have open files, now we need to read the file first.
  string line; //out_line;

  string comment;
  string gene_sequence;
  string gene_count;
  
  int gene_number=0;
  
  int line_count=0;
  string temp_seq;

  cout<<"loading the data file........."<<endl;
  cout<<"line#: ";
  //  double d_storage, d_current;
  string temp_string;
  bool header_is_read=false;

  while(!ifs_p.eof())
    {
      line_count++;
      if(line_count/5000*5000==line_count)
	{
	  cout<<"..... "<<line_count;
	  cout.flush();
	}
      //found_one=false;
      getline(ifs_p, line);
      if(line.compare("\n")==0||line.compare("\r")==0||line.compare("\r\n")==0||
	 line.compare("\n\r")==0||line.compare(" ")==0||line.compare("  ")==0||
	 line.length()==0||line.compare("\t")==0||line.compare("\t\t")==0)
	{
	  cout<<"...read an empty row (skip!)"<<endl;
	  continue;
	}
            //cout<<"here......"<<endl;
      if(line[0]==_commentChar)  //header file
	{
	  //this is a new refSeq header line,
	  temp_string=line.substr( 1,line.length());
	  _header.push_back(temp_string);
	}
      else //sequence line, each line contains the full sequence plus the count
	//atata...ata\t50
	{	
	  if(_data_header_line&&(!header_is_read))
	    {
	      header_is_read=true;
	      continue;
	    }
	  //need to split in order to parse the sequence and count
	  string buffer[100];
	  int num= split_ext(line, buffer, _delimiter);
	  if(num>2||num<1)
	    {
	      cout<<"****ERROR:Unknow format for squence data file, exit!"<<endl;
	      exit(-1);
	    }
	  stringstream sname; sname<<"seq "<<_seq.size()+1;
	  SequenceString ss(sname.str(), to_upper_str(buffer[0]));
	  _seq.push_back(ss);
	  if(num==2)
	    {
	      _count.push_back(atoi(buffer[1].c_str()));
	    }
	  else
	    {
	      _count.push_back(0);
	    }
	  //_seqStrVec.push_back(ss);
	      
	
	  //gene_info=temp_string;
	  //temp_seq.clear();
	  //gene_number++;
	  
	}
      
    }
  
  //here we have to update/add the last gene
  //read_gene_sequence(temp_seq, gene_sequence)==-1)
  //SequenceString ss(gene_info,to_upper_str(temp_seq));
  //_seqStrVec.push_back(ss);
  

  //we are done with reading the promoter sequence file.
  cout<<"\nfinish reading the promoter file........\n"
      <<"\tsummary: total "<<line_count<<" line read in and \n"
      <<"\t\t"<<gene_number<<" sequences store in the vector...."<<endl;
  
  ifs_p.close();
  return _count.size();
}
unsigned ParseField(vector<string>& _header, const string& _field)
{
  //go through the vector and select the field and return the value in the header field
  string headerString="";
  for(unsigned i=0;i<_header.size();i++)
    {
      //string compare 
      if(_header.at(i).find(_field)!=string::npos)
	{
	  //found one, we can jump out
	  headerString=_header.at(i);
	  break;
	}
    }
  if(headerString.size()==0)
    {
      //can not find any,
      return 0;
    }
  //now split and parse the field number
  string buffer[10];
  unsigned num=split_ext(headerString, buffer,'=' );
  if(num!=2)
    {
      cout<<"****WARNING: unknown format for the header in data file, return 0 by default"<<endl;
      return 0;
    }
  return atoi(buffer[1].c_str());
}
