#include <iostream>
#include <fstream>
#include <vector>
#include "../string_ext.hpp"

//from on you need to specify the c libraries explicitly, 11/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
//#include <pthread.h>

#include "../SequenceString.hpp"
#include "../FastaHandler.hpp"
//#include "FileManipulator.hpp"
//#include "Concensus.hpp"
#include "../SIGPIG/LoadData.hpp"

using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);

string seqFileName("");
string outFileName("");

int main(int argc, char* argv[])
{
  //for parsing commandline arguement

  const char *opts = "hs:o:";
  parseArguments(argc, argv, opts);

  if(seqFileName.length()==0)
    {
      cout<<"please specify the input file........quit."<<endl;
      exit(-1);
    }
  if(outFileName.length()==0)
    {
      outFileName=seqFileName+".fasta";
    }

  //now check for validity
  ifstream ifs_ss(seqFileName.c_str(), std::ios::in);
  if(!ifs_ss.is_open())
    {
      cout<<"******ERROR: can not open file, quit..."<<endl;
      exit(-1);
    }
  ifs_ss.close();
  ofstream ofs_ss(outFileName.c_str(), std::ios::out);
  if(!ofs_ss.is_open())
    {
      cout<<"********Error: can not open out file, quit......"<<endl;
      exit(-1);
    }
  ofs_ss.close();

  //now start loading input
  vector<SequenceString> data_vec;
  vector<string> header_vec;
  vector<unsigned> count_vec;
  unsigned totalNumSequences=LoadData(seqFileName, header_vec, data_vec, count_vec);
  cout<<"Load Data: "<<totalNumSequences<<" sequences were read"<<endl;
  
  cout<<"header vector"<<endl;
  for(unsigned int i=0;i<header_vec.size();i++)
    {
      cout<<"\t"<<i<<":"<<header_vec.at(i)<<endl;
    }

  cout<<"count vector"<<endl;
  for(unsigned int i=0;i<count_vec.size();i++)
    {
      cout<<"\t"<<i<<":"<<count_vec.at(i)<<endl;
    }
  
  cout<<"sequence vector"<<endl;
  for(unsigned int i=0;i<data_vec.size();i++)
    {
      cout<<"\t"<<i<<":"<<data_vec.at(i).toString()<<endl;
    }

  //now start saving into fasta format
  WriteFasta(outFileName, data_vec);

  cout<<"Done..!!"<<endl;

  return 0;
}


static void parseArguments(int argc, char **argv, const char *opts)
{
  int optchar;

  while ((optchar = getopt(argc, argv, opts)) != -1)
    {
      switch(char(optchar))
	{
	case 's':
	  seqFileName=optarg;
	  break;
	  
	case 'o':
	  outFileName=optarg;
	  break;

	
	case '?':
	  /*if(optopt == 't')
	    ;//cout<<"option"<<optopt<<" requires an argument.\n"<<endl;
	    else*/
	  if(isprint(optopt))
	    {
	      cout<<"Unknown option "<<optopt<<endl;
	    }
	  else
	    cout<<"Unknown option character"<<endl;
	
	case 'v':  
	case 'h': // usage or unknown
	default:
	  printUsage(argc, argv);
	  exit(-1);
	}
    }
}


static void printUsage(int argc, char* argv[])
{
  cout<<"Usage:\n";
  cout<<"\t"<<argv[0]<<"-s sequence inputfile \n"
      <<"\t [-o output filename]\n"
    ;

  cout<<"\t\t-s filename -- the sequence file as the input \n"
      <<"\n";

  cout<<"\t\t-o filename -- the output file name\n"
      <<"\n";


  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to reformate the sequence inputfile to fasta format\n"
      <<"\t\t\t so that the input sequence file can be read/used by other software,\n"
      <<"\t\t\t mainly by clonanyst for alignment.\n" 
    <<"\t\t\t @3/20/2015\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********by Feng @ BU 2015\n"; 

  exit(-1);
}

