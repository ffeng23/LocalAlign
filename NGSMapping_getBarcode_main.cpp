#include <iostream>
#include <string.h>
#include <vector>
#include <zlib.h>
#include "Accessory/string_ext.hpp"
#include "Accessory/FileHandler.hpp"

//from now on you need to specify the c libraries explicitly
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "FASTQ.hpp"
#include "FastqHandler.hpp"
#include "SequenceHandlerBarcode.hpp"
#include "FastaHandler.hpp"
using namespace std;

//global function for help and input
static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);

//define global varible mainly for commandline input

static string sequenceFile_name;//input file for sequence data

//static unsigned int MinimumOverlapLength=10;//not too short
static bool dualIndex=false;
static string indexFileR1_name;
static string indexFileR2_name;
static bool indexFromSeq=false;
static unsigned lenOfBarcode=8;

//note for match  matrix
//the match matrix mmf has been declare and initialized in 
//score.hpp and score.cpp, respectively. we can use it 
//directly
//   MatchMatrix mmf; 

int main(int argc, char* argv[])
{
	//for parsing commandline arguement
  const char *opts = "hvs:di:j:nl:";
  //d: dual index, indicating the inpute is dual index
  //s: sequence data file, contains the real sequence data, having the same order and name as index files.
  //		in this programe, the sequence as input only because we want to read index/barcodes 
  //		from the name of the sequences. 
  //n: get index from the barcode 
  //	s and n should be set together to work.
  //i,j: index data file, optional, could be specified or not. for dual index mode, we will reach index 
  //		and get barcodes. NOT IMPLEMENTED YET!!!!
  //l: length of barcode, 8 nts
  // v and h: version and help

  cout<<"Start parsing the input commands......."<<endl;
  parseArguments(argc, argv, opts);
  
  cout<<"***Input parameters Summary:\n";
  //start print out the input and check for validaty of inputs
  if(indexFromSeq)
  {
	  cout<<"\tRead index from sequence names:TRUE\n";
	if(sequenceFile_name.length()==0)
	{
	  cout<<"ERROR:no sequence input file specified!!!"<<endl;
	  exit(-1);
	}
	cout<<"\tSequence file name:"<<sequenceFile_name<<endl;
  }
  else //now we need to set up the 
  {
	  cout<<"\tRead index from sequence names:FALSE\n";
	if(indexFileR1_name.length()==0)
	{
	  cout<<"ERROR:no index (1) input file specified!!!"<<endl;
	  exit(-1);
	}
	cout<<"\tIndex R1 file name:"<<indexFileR1_name<<endl;
	if(dualIndex)
	{
		if(indexFileR2_name.length()==0)
		{
		  cout<<"ERROR:no index (2) input file specified on a dual index run!!!"<<endl;
		  exit(-1);
		}
		cout<<"\tIndex R2 file name:"<<indexFileR2_name<<endl;
	}
	cout<<"NOTE::::::This option has NOT been implemented.\n Please run read index from sequence mode\n "<<endl;
	return 0;
  }
  cout<<"\tDual index:"<<dualIndex<<endl;
  cout<<"\tLength of barcode:"<<lenOfBarcode<<endl;
  
  //check file format, so far we only read either fastq or gzip'ed fastq
  FileType ft=getFileType(sequenceFile_name);
  if(ft!=FASTQ&&ft!=GZ_FASTQ)
  {
	cout<<"ERROR: so far we only process the fastq or gzip'ed fastq files"<<endl;
	return 0;
  }
  //Start testing the readfastq function
  cout<<"---------start testing read fastq (either gz'ed or regular fastq) ......."<<endl;
  vector<Fastq> vec;
  
  unsigned int num=ReadFastq(sequenceFile_name, vec, false);
  cout<<"total number of sequences read: "<<num<<endl;
  
  vector<SequenceString> v_seq;
  //now we need to print out for debugging.
  for(unsigned i=0;i<vec.size();i++)
  {
	  //cout<<"i:"<<i<<"--"<<vec.at(i).toString()<<endl;
	  v_seq.push_back(vec.at(i).GetSequenceString());
  }
  
  cout<<"DONE"<<endl;
  
  cout<<"=====Start proessing indexes ..........."<<endl;
  
  //first get the indexes from the header of the sequences.
  vector<SequenceString> v_index1;
  vector<SequenceString> v_index2;
  
  unsigned lenOfBarcode=8;
  bool dualIndex=true;
  
  cout<<"Reading indexes from Sequence names......."<<endl;
  unsigned numOfSeq=ReadIndexFromSequenceName(v_seq, v_index1, v_index2, lenOfBarcode, dualIndex);
  
  cout<<"\tnumber of sequences processed:"<<numOfSeq<<endl;
  
  /*for(unsigned i=0;i<v_index1.size();i++)
  {
	  cout<<v_index1.at(i).toString()<<endl;
	  cout<<v_index2.at(i).toString()<<endl;
  }*/
  
  //now we have the indexes, pick the unique ones as barcodes from them
  vector<SequenceString> v_BarSeq1;
  vector<SequenceString> v_BarSeq2;
  vector<unsigned> v_count;
  cout<<"Counting unique barcodes from indexes......."<<endl;
  unsigned numOfBar=GetBarcodes(v_index1, v_index2, dualIndex, v_BarSeq1, v_BarSeq2, v_count);
  cout<<"\n\tnumber of bar:"<<numOfBar<<endl;
  
  cout<<"Start writting output files..........."<<endl;
  //WriteFasta("index1.fasta", v_index1);
  //WriteFasta("index2.fasta", v_index2);
  WriteFasta(sequenceFile_name+"bar1.fasta", v_BarSeq1);
  WriteFasta(sequenceFile_name+"bar2.fasta", v_BarSeq2);
  WriteTextFile(sequenceFile_name+"count.txt", v_count);

  //make it ready as a table
  vector<vector<string> > stat_vec;
  vector<string> bar1;
  vector<string> bar2;
  vector<string> count_str;
  for(unsigned i=0;i<v_BarSeq1.size();i++)
  {
	  bar1.push_back(v_BarSeq1.at(i).GetSequence());
	  if(dualIndex)
	  {
		  bar2.push_back(v_BarSeq2.at(i).GetSequence());
	  }
	  count_str.push_back(to_string((v_count.at(i))));
  }
  stat_vec.push_back(bar1);
  if(dualIndex)
	  stat_vec.push_back(bar2);
  stat_vec.push_back(count_str);
  vector<string> header;
  header.push_back("R1_Barcode");
  if(dualIndex)
  {
	  header.push_back("R2_Barcode");
  }
  header.push_back("count");
  
  //now calling to get it written down to the disk
WriteTextTableFile(sequenceFile_name+"_stat.txt",stat_vec, '\t',true, ios_base::out,header); 
  
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
	case 'i':
	  indexFileR1_name=optarg;
	  break;
	
	case 'j':
	  indexFileR2_name=optarg;
	  break;
	   
	case 's':
	  sequenceFile_name=optarg;
	  break;
	  
	case 'n':
	  indexFromSeq=true;
	  break;
	  //case 'i':
	  //isotype_flag=true;
	  //break;

	case 'd':
	  dualIndex=true;
	  break;
	
	case 'l':
	  lenOfBarcode=atoi(optarg);
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
  cout<<argv[0]<<", the program used to go through the index or sequence (name for index)\n\n"
	  <<"to collect the unique barcodes and its states (counts). It would read in\n"
      <<"indexes and insert into a vector and count them uniquely. It could read either gz'ed\n"
      <<"or regulare fastq. So far it only support fastq. WILL DO FASTA laster!!!\n"
	  <<"We could do index from sequence name (-s and -n together ) or from index files\n"
      <<"(index 1 and index 2, -i and -j and -d). For now, only read from sequence name\n"
	  <<". WILL DO INDEX FILE LATER!!!\n"
      <<"To find the unique barcodes, we simply assuming the length is fixed and we will \n"
	  <<"compare the strings (string.compare), not doing alignment or match. we do insertion\n"
	  <<"sort kind of algorithem. the obtained barcodes are sorted and counted"
	  <<"the length of barcode input is only used by reading indexes from, but not for getting\n"
	  <<"the unique barcodes.\n";

  //d: dual index, indicating the inpute is dual index
  //s: sequence data file, contains the real sequence data, having the same order and name as index files.
  //		in this programe, the sequence as input only because we want to read index/barcodes 
  //		from the name of the sequences. 
  //n: get index from the barcode 
  //	s and n should be set together to work.
  //i,j: index data file, optional, could be specified or not. for dual index mode, we will reach index 
  //		and get barcodes. NOT IMPLEMENTED YET!!!!
  //l: length of barcode, 8 nts
  // v and h: version and help
  cout<<"\tOption string: \"hvs:di:j:nl:\" \n";
  cout<<"Usage:\n";
  cout<<"\t"<<argv[0]<<" [-s sequence file] [-n]  [-d] [-l length of barcode]\n"
      <<"\t\t [-i  index R1 file] [-j index R2 file]\n"
       <<"\t\t [-h] [-v]\n\n"
    ;
	cout<<"\t-s: sequece data file, contains the real sequence data, \n"
      <<"\t\thaving the same order and name. This input is only good \n"
	  <<"\t\tbecause we need to read indexes from the name of the sequences\n"
	  <<"\t\tit has to be set together with -n reading index from name\n"
	  <<"\t\tNOTE:: CURRENTLY ONLY SUPORT FASTQ OR GZ'ed FASTQ file!!!!!\n"
		;
  cout<<"\t-i: the input file name read 1 index\n"
      <<"\t-j: the input file name read 2 index\n" 
      <<"\t\ti and j contain the indexes from reads. These are require, if\n"
	  <<"\t\tno -n is set. NOTE::::THIS HAS NOT BEEN IMPLEMENTED!!!\n";
  //cout<<"\t-p: paired end read, boolean, false by default\n" ;
  cout<<"\t-n: read indexs from the names of each sequence and assuming \n"
      <<"\t\tthe illumina style:\n"
      <<"\t\t\t  xxxxxx:xxxxx:xxxx:actgkd+acdfad\n" 
      <<"\t\t\t\tseparated by ':' and the last field contains the barcode (pair)\n"
      <<"\t\t\t\tfalse by default\n";
  
  cout<<"\t-d: dual index (true or false), true by default\n"
      <<" \t\tor simple output the stats";
  cout<<"\t-l: length of barcodes, used to read indexes form sequence names\n"
		<<"\t 8 by default\n"
      ; 

  cout<<"\t-v and -h: version and help\n";
  cout<<"Note: 1)when do string comparison for uniqueness of index/barcode"
      <<"not anything fancy like matching or alignment\n"
      <<"\t\t\t *************@1/3/2019 by Feng\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********updated by Feng @ BU 2019\n"; 

}
