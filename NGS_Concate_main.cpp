#include <iostream>
#include <fstream>
#include <vector>
#include <zlib.h>
#include "Accessory/string_ext.hpp"

//from now on you need to specify the c libraries explicitly
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "Accessory/SequenceString.hpp"
#include "Accessory/FastaHandler.hpp"
#include "SequenceHandlerCommon.hpp"

#include "Accessory/FileHandler.hpp"
#include "Accessory/FASTQ.hpp"
#include "Accessory/FastqHandler.hpp"

using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);
void printCallingCommand(int argc, char* argv[]);

static string fileR1_name;//default input file for single index, R1
static string fileR2_name;
static string outFile_name;//for index R1 output file 

//enum mapType {FivePrime, ThreePrime}; //defined in SequenceHandlerBarcode.hpp
static bool rc=false; //only work for R2, to reverse complement the sequence  R2,
                                 //we never reverse complement the sequence R1. 

int main(int argc, char* argv[])
{
  printCallingCommand(argc, argv);
  //for parsing commandline arguement
  const char *opts = "hvs:t:ro:";
  //r: reverse complement the read2 barcode or not to do demux
  //s: sequence data file, contains the real sequence data, having the same order and name.
  //t: sequence data file, optional, could be specified or not. if dumux, then needed to check for this.
  // v and h: version and help

  cout<<"Start parsing the input commands......."<<endl;
  parseArguments(argc, argv, opts);
  
  cout<<"***Input parameters Summary:\n";
  
  if(fileR1_name.size()==0)
    {
      cout<<"please specify the input file name (R1).......\n";
      //printUsage(argc, argv);
      cout<<"type "<<argv[0]<<" -h \" for usage help\n";
      exit(-1);
    }
  if(fileR2_name.size()==0)
    {
      cout<<"please specify the input file name (R2).......\n";
      //printUsage(argc, argv);
      cout<<"type "<<argv[0]<<" -h \" for usage help\n";
      exit(-1);
    }
  cout<<"file (Read 1):\""<<fileR1_name<<"\".\n";
  cout<<"file (Read 2):\""<<fileR2_name<<"\".\n";
  
  cout<<"reverse complement of R2:"<<rc<<endl;
  
  //*********get output files names
  if(outFile_name.size()==0)
	outFile_name.assign(fileR1_name+".Conc");
  
  //check the input file types
  FileType ft1=getFileType(fileR1_name,true);
  FileType ft2=getFileType(fileR2_name,true);
  if(ft1!=ft2)
  {
	  cout<<"Error: file types are not match. Please check...."<<endl;
	  exit(-1);
  }
  cout<<"File Type:";
  bool unknownType=false;
  switch (ft1)
  {
	  case FASTA:
		cout<<"FASTA";
		outFile_name.append(".fasta");
		break;
	  case GZ_FASTA:
		cout<<"gzip FASTA (compressed)";
		outFile_name.append(".fasta");
		break;
	  case FASTQ:
		cout<<"FASTQ";
		outFile_name.append(".fastq");
		break;
	  case GZ_FASTQ:
		cout<<"gzip FASTQ (compressed)";
		outFile_name.append(".fastq");
		break;
	  case GZ:
		cout<<"compressed file";
		unknownType=true;
		break;
	  case TXT:
		cout<<"text file";
		unknownType=true;
		break;
	  case UNKNOWN:
	  default:
		cout<<"UNKNOWN file type";
		unknownType=true;
		break;
  }
  cout<<endl;
  if(unknownType)
  {
	  cout<<"The file type is not supported, please input file of (gzip) fasta/fastq files"<<endl;
	  exit(-1);
  }
  if(outFile_name.compare(fileR1_name)==0)
  {
	  size_t dot_pos=outFile_name.find_last_of('.');
	  if(dot_pos==string::npos)
	  {
		  outFile_name.append("_2");
	  }
	  else
	  {
		  outFile_name.insert(dot_pos, 1,'2');
	  }
  }
  // the output file names  
  //outFileR1_name=outFileR1_name.append(".conc.
  cout<<"The output file name : \""<<outFile_name<<"\"\n";
  
  //start doing the reading....assuming we only deal with fasta, fastq and fastq gzipped. no others
  //the file type only  
  cout<<"--------------------\n";
  cout<<"Starting reading from input file ad concatenating.....\n";
  size_t pt=0;

  if(ft1==GZ_FASTA||ft1==FASTA)
  {
	 pt=concatenateFasta(fileR1_name, fileR2_name,ft1,
	     outFile_name,ios_base::out,rc);
  }
  else 
  {
	  if(ft1==GZ_FASTQ||ft1==FASTQ)
	  {
		  //cout<<"calling here......"<<endl;
		 pt=concatenateFastq(fileR1_name, fileR2_name, ft1,outFile_name,ios_base::out,rc);
	  }
  }
  cout<<"Total "<<pt<<" records processed"<<endl;
  cout<<"Done!!!"<<endl;
  
  cout<<"Thanks for using our program and have fun!!"<<endl;
  return 0;
}


static void parseArguments(int argc, char **argv, const char *opts)
{
  int optchar;
  //int temp;
  while ((optchar = getopt(argc, argv, opts)) != -1)
    {
      switch(char(optchar))
	{	  
	case 's':
	  fileR1_name=optarg;
	  break;
	case 't':
	  fileR2_name=optarg;
	  break;
	case 'r':
	  rc=true;
	  break;
	case 'o':
		outFile_name=optarg;
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
  cout<<argv[0]<<", the program used to concatenate two fasta/fastq files."
	  <<"Assume the two files are read 1 and read 2 of the same type. "
	  <<"We also reverse complement the read 2 file by deafault."
	  <<"Currently support (fastq/fasta/fq/fa).(gz/gzip).\n";

  cout<<"\tOption string: \"hvs:t:r\" \n";
  cout<<"Usage:\n";
  cout<<"\t"<<argv[0]<<" -s sequence file (R1) -t sequence file r2\n"
	  <<"\t\t [-r (reverse implement r2 sequence index)]\n"
	  <<"\t\t [-o output file name prefix]\n"
      <<"\t\t [-h] [-v]\n\n"
    ;
  cout<<"\t-s: the input barcode file name (for r1) \n";
  cout<<"\t-t: the input file barcode name (for r2) \n"
      <<"\t\tassuming has identical name and order\n"
	  <<"\t-r: reverse complement the read2 sequences\n"
      ;
  cout<<"\t-v and -h: version and help\n";
  cout<<"Note: 1)when we do the concatenation, we reverse complement the r2 and \n"
	  <<"reverse the quality string if they are the fastq files.\n";
	
  cout<<"2)The output file are either fasta or fastq. \n"
	  <<"3)The output file name is either determined from user defined by \"-o\"\n"
	  <<"\t or it is determined using the R1 input file name. \n"
	  <<"\twhen user input the file name, we assume it is the prefix without the file type.\n"
	  <<"\tthe file type is determined by the input file type.\n"
	  ;
	  
  cout<<"Example:\n"
		<<"\t"<<argv[0]<<"-s file1.fasta(.gz) -t file2.fasta(.gz) -r\n"
		;
  cout<<"**Note**:\n"
	<<"\tthis program could be used to add the umi barcodes to the beginning of R1 file \n"
	<<"\tin case of nextseq data, which has R2 as the umi and R1 and R3 are the regular \n"
	<<"\tR1 and R2 input file, respectively.\n"
	;
  cout   <<"\t\t\t *************@9/2/2019 by Feng\n";
  cout<<"\t\t\t\t *****updated 12/29/2019 by Feng: adding output file name -o\n";
     
}

void printCallingCommand(int argc, char* argv[])
{
	cout<<"Calling:"<<endl;
	for(int i =0;i<argc;i++)
	{
		cout<<argv[i]<<" ";
	}
	cout<<endl;
}

