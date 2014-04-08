#include "FastaHandler.hpp"
#include <stdlib.h>
//#include <stdio.h>
#include <fstream>
#include "string_ext.hpp"
using namespace std;
unsigned int ReadFasta(const string& _fname, vector<SequenceString>& _seqStrVec, bool toUpper)
{
  ifstream ifs_p(_fname.c_str());
      
  if(!ifs_p.is_open())
    {
      cout<<">>>>>>ERROR:the input file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  

  //so far, we have open files, now we need to read the file first.
  string line; //out_line;

  string gene_info;
  string gene_sequence;
  
  int gene_number=0;
  
  int line_count=0;
  string temp_seq;

  cout<<"read the fasta file........."<<endl;
  cout<<"line#: ";
  //  double d_storage, d_current;
  string temp_string;
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
      if(line.compare("\n")==0||line.length()==0||line.compare("\t")==0||line.compare("\t\t")==0)
	{
	  cout<<"...read an empty row (skip!)"<<endl;
	  continue;
	}
            //cout<<"here......"<<endl;
      if(line[0]=='>')
	{
	  //this is a new refSeq header line,
	  temp_string=line.substr( 1,line.length());
	  	  	  
	  if(gene_number>0) //we already read something, so we have to store it to the vector
	    {

	      SequenceString ss(gene_info,to_upper_str(temp_seq));
	      _seqStrVec.push_back(ss);
	      
	    }//otherwise, we are doing the first one, so do bother
	  gene_info=temp_string;
	  temp_seq.clear();
	  gene_number++;
	  
	}
      else
	{//append the seqence to it
	  temp_seq.append(line);
	  //temp_seq.append("\n");
	}

    }
  
  //here we have to update/add the last gene
  //read_gene_sequence(temp_seq, gene_sequence)==-1)
  SequenceString ss(gene_info,to_upper_str(temp_seq));
  _seqStrVec.push_back(ss);
  

  //we are done with reading the promoter sequence file.
  cout<<"\nfinish reading the promoter file........\n"
      <<"\tsummary: total "<<line_count<<" line read in and \n"
      <<"\t\t"<<gene_number<<" sequences store in the vector...."<<endl;
  
  ifs_p.close();
  return gene_number;
}

void WriteFasta(const string& _fname, vector<SequenceString>& _seqStrVec, const unsigned int _width,  ios_base::openmode mode)
{
  ofstream ofs((_fname).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  
  for(unsigned int i=0;i<_seqStrVec.size();i++)
    {
      ofs<<">"<<_seqStrVec.at(i).GetName()<<"\n";
      for(unsigned int j=0;j*_width<_seqStrVec.at(i).GetLength();j++)
	{
	  ofs<<_seqStrVec.at(i).GetSequence().substr(j*_width, _width)<<endl;
	}
    }

  ofs.close();
}

void WriteTextFile(const string& _fname, vector<int>& _seqStrVec, const char& c, const unsigned int _width, ios_base::openmode mode)
{
  ofstream ofs((_fname).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  
  for(unsigned int i=0;i<_seqStrVec.size();i++)
    {
      //ofs<<">"<<_seqStrVec.at(i).GetName()<<"\n";
      ofs<<_seqStrVec.at(i);
      if((i+1)%_width==0||i==_seqStrVec.size()-1)
	{
	  ofs<<"\n";
	}
      else
	{
	  ofs<<c;
	}
    }

  ofs.close();
  
}

