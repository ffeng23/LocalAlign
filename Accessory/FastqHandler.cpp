#include "FastqHandler.hpp"
#include <stdlib.h>
//#include <stdio.h>
#include <fstream>
#include <zlib.h>
#include "string_ext.hpp"
using namespace std;

//not very efficient. might need to rewrite. but how to avoid the code duplication?
size_t ReadFastq(const string& _fname, vector<SequenceString>& _seqStrVec, vector<string>& _vecQ, bool toUpper)
{
		vector<Fastq> v_fastq;
		size_t n=ReadFastq(_fname, v_fastq, toUpper);
		_seqStrVec.clear();
		_vecQ.clear();
		if(v_fastq.size()>0)
		{
			for(unsigned i =0;i<v_fastq.size();i++)
			{
				_seqStrVec.push_back(v_fastq.at(i).GetSequenceString());
				_vecQ.push_back(v_fastq.at(i).GetQualityString());
			}
		}
		return n;
}
//return string::npos upon error 
size_t ReadFastq(const string& _fname, vector<Fastq>& _seqStrVec, bool toUpper)
{
	//need to know whether this is compressed file;
	//we are using the ".gz" suffix as a criterion for the compressed file.
	//we only do gzip compression.not others.
	bool gzflag=_fname.substr(_fname.size()-3, 3)==".gz";
	ifstream ifs_p;
	
	FILE* fb=NULL;
	//cout<<"file name is "<<_fname<<endl;
	if(gzflag)	
	{
		//cout<<"....1"<<endl;
		fb=gzOpen_B(_fname, "rb");
		if(!fb)
		{
			cout<<">>>>>>ERROR:input file\""<<_fname<<"\" can not be opened, quit....\n"<<endl;
			return string::npos;
		}
	}
	else
	{
		ifs_p.open(_fname);
		//cout<<"......2"<<endl;
		if(!ifs_p.is_open())
		{
			cout<<">>>>>>ERROR:the input file \""<<_fname<<"\" can not be opened, quit....\n"<<endl;
			return string::npos;
		}
	}  
	//cout<<"....4 : ifs good :"<<ifs_p.good()<<endl;
	//cout<<".....3:ifs_p_open:"<<ifs_p.is_open()<<endl;
  //so far, we have open files, now we need to read the file first.
  string line; //out_line;

  string gene_info;
  string gene_sequence;
  string quality;
  int gene_number=0;
  
  int line_count=0;
  string temp_seq;

  cout<<"read the fastq file.........";
  //cout<<"line#: ";
  //  double d_storage, d_current;
  string temp_string;
  bool ok=true;
  while(ok)
    {
      
	  //start reading lines and parsing
	  //fq files line 1: name
	  //	line 2: seqence
	  //	line 3: +
	  //	line 4: quality
	  //
	  //cout<<"start the loop......."<<endl;
	  if(gzflag) //compressed file 
	  {
		//cout<<"\tready to read........"<<endl;
		ok=getline_B(fb, line);
	  }
	  else
		getline(ifs_p, line);  
      
	  chomp_ext(line);
	  //cout<<"i:"<<line_count<<"--line:"<<line<<endl;

      if(line.compare("\n")==0||line.length()==0||line.compare("\t")==0||line.compare("\t\t")==0)
		{
		  cout<<"...read an empty row (skip!)"<<endl;
		  //last, check the file end 
		  if(!gzflag)
			ok=!ifs_p.eof();
		  continue;
		}
        //first line has to be "@xxxxxxx    
		//cout<<"here......"<<endl;
      if(line[0]=='@')
		{
		  //this is a new refSeq header line,
		  gene_info=line.substr( 1,line.length());
		  gene_number++;
		}
      else
		{//
			cout<<"ERROR: corrupted file, no lines with '@...' please check! Reading is stopped"<<endl;
			
			/*if(gzflag)
			  {
				gzClose_B(fb);
			  }
			  else
			  {  
				ifs_p.close();
			  }*/
			break;
			
		}
		
	//last, check the file end
		if(!gzflag)	
			ok=!ifs_p.eof();
		
	  //now doing the line 2
	  //sequence
	  if(gzflag) //compressed file 
	  {
		ok=getline_B(fb, line);
	  }
	  else
		getline(ifs_p, line);  
      
	  chomp_ext(line);
	  if(toUpper)
		gene_sequence.assign(to_upper_str(line));
	  else
		gene_sequence.assign(line);

	//cout<<"\t"<<gene_sequence<<endl;
	//last, check the file end 
		if(!gzflag)
			ok=!ifs_p.eof();
		
	  //reading the third line
	  if(gzflag) //compressed file 
	  {
		ok=getline_B(fb, line);
	  }
	  else
		getline(ifs_p, line);  
	
	  chomp_ext(line);
	  
	  if(line[0]!='+'||line.size()!=1)
		{
		  //there is something wrong. we need to quit
		  cout<<"ERROR: something wrong. the file is corrupted!! Reading is stopped!!"<<endl;
		  break;
    	}
	//cout<<line<<endl;
	//last, check the file end 
		if(!gzflag)
			ok=!ifs_p.eof();	
	
	//reading the 4th line
		//reading the third line
	  if(gzflag) //compressed file 
	  {
		ok=getline_B(fb, line);
	  }
	  else
	  {
		getline(ifs_p, line);  
	  }
	  chomp_ext(line);
	  
	  quality.assign(line);
		
	 //cout<<line<<endl; 
	
	  //done, now we need to push to the vector 
	  Fastq fq(gene_info,SequenceString(gene_info, gene_sequence),quality);
	  _seqStrVec.push_back(fq);
	  line_count++;
      if(line_count%10000==0)
		{
		  cout<<"..... "<<line_count;//<<endl;
		  cout.flush();
		}
	  //last, check the file end 
		if(!gzflag)
			ok=!ifs_p.eof();   //we don't have to check for gzipped case, since it will return false anyway.
    }
  

  //we are done with reading the promoter sequence file.
  cout<<"\nfinish reading the file........\n"
      <<"\tsummary: total "<<line_count<<" records read in and \n"
      <<"\t\t"<<gene_number<<" sequences store in the vector...."<<endl;
  if(gzflag)
  {
	gzClose_B(fb);
  }
  else
  {  
	ifs_p.close();
  }
  return gene_number;
}

void WriteFastq(const string& _fname, vector<Fastq>& _fqVec,  ios_base::openmode mode)
{
  ofstream ofs((_fname).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  
  for(unsigned int i=0;i<_fqVec.size();i++)
    {
      ofs<<"@"<<_fqVec.at(i).GetName()<<"\n";
      
	  ofs<<_fqVec.at(i).GetSequence()<<"\n";
	  
	  ofs<<"+\n"<<_fqVec.at(i).GetQualityString()<<"\n";
	
    }

  ofs.close();
}

void WriteFastq(const string& _fname, vector<SequenceString>& _seqStrVec, 
	vector<string>& _qVec,
	ios_base::openmode mode)
{
  ofstream ofs((_fname).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  if(_qVec.size()!=_seqStrVec.size())
  {
	  cout<<"ERROR: the sizes of the sequence string and qualtiy string don't match. please check........"<<endl;
  }
  for(unsigned int i=0;i<_qVec.size();i++)
    {
      ofs<<"@"<<_seqStrVec.at(i).GetName()<<"\n";
      
	  ofs<<_seqStrVec.at(i).GetSequence()<<"\n";
	  
	  ofs<<"+\n"<<_qVec.at(i)<<"\n";
	
    }

  ofs.close();	
		
}
