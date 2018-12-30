#include <iostream>
#include <fstream>
#include <vector>
#include <zlib.h>
#include "Accessory/string_ext.hpp"

//from now on you need to specify the c libraries explicitly
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "score.hpp"
#include "SequenceString.hpp"
//#include "OverlapAlignment.hpp"
#include "FastaHandler.hpp"
//#include "SequenceHandlerIsotype.hpp"
#include "SequenceHandlerCommon.hpp"
#include "SequenceHandlerBarcode.hpp"

using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);
//static int lookUpScoreMatrix(const string* _scoreMatrixNameArray,const int& len, const string& scoreMatrixName);

//all file are in fasta format
static string barFileR1_name;//default input file for single index, R1
static string barFileR2_name;
static string indexFileR1_name;
static string indexFileR2_name;
//static string outFile_name;//for stats;

static string sequenceFileR1_name;//input file for sequence data
static string sequenceFileR2_name;//read 2

static unsigned  misMatchNum=1; //not too many mismatch 

//static unsigned int MinimumOverlapLength=10;//not too short
//enum mapType {FivePrime, ThreePrime}; defined in SequenceHandlerBarcode.hpp
static mapType mapEnd=FivePrime; //only work for R2, to reverse complement the sequence index 2,
                                 //we never reverse complement the barcode. 
static bool pairEnd=false;
static bool indexFromSeq=false;
static bool demux=false;//write stats only, no demux
static bool dualIndex=false;

//note for match  matrix
//the match matrix mmf has been declare and initialized in 
//score.hpp and score.cpp, respectively. we can use it 
//directly
//   MatchMatrix mmf; 


int main(int argc, char* argv[])
{
  //for parsing commandline arguement
  const char *opts = "hva:b:f:e:ptdxrs:q:m:";
  //a: the input file barcode name (r2 for dual index input), assuming has identical name and order
  //b: the input barcode file name (r1 or for single index input
  //       a and b list the expected barcode to demux the input
  //f: the input file name read 1 index
  //e: the input file name read 2 index, 
  //       f and r contains the read index from sequencing to be demux'ed
  //p: paired end read, boolean 
  //t: read indexs from the names of each sequence and assuming the illumina style
  //         xxxxxx:xxxxx:xxxx:actgkd+acdfad, separated by ':' and the last field 
  //         contains the barcode (pair)
  //x: demux, meaning to write out the sequences or simple output the stats
  //d: dual index
  //r: reverse complement the read2 barcode or not to do demux

  //s: sequence data file, contains the real sequence data, having the same order and name.
  //q: sequence data file, optional, could be specified or not. if dumux, then needed to check for this.
  //m: mismatch, allow degenerate barcode, meaning the barcode (expected) could be degenerated, but not the other way around
  // barcode (expected)-->row, index (sequence read) ---->column
  //      A C G T 
  //note: 1)when do comparison, always use barcode(expected) to align index (data), in terms 
  //          of alignment, barcode is sujbect, index is pattern
  //             index (pattern)
  //             |||||
  //             barcode (subject)
  //      2)allow degeneracy on barcode, but less so on index
  //      3)it is possible to have index longer than barcode, but not vice versa
  //      4)barcode and index must start on 5' position 1 (position 0 in c)
  //      5)we don't do alignment in comparison, but string comparison, see comparison matrix
  //                  match.matrix.feng
  // v and h: version and help

  cout<<"Start parsing the input commands......."<<endl;
  parseArguments(argc, argv, opts);
  
  cout<<"***Input parameters Summary:\n";
  
  if(barFileR1_name.size()==0)
    {
      cout<<"please specify the barcode  data input fasta file name.......\n";
      //printUsage(argc, argv);
      cout<<"type "<<argv[0]<<" -h \" for usage help\n";
      exit(-1);
    }
  cout<<"barcode file (Read 1):\""<<barFileR1_name<<"\".\n";
  //dual index
  if(dualIndex)
    {
      cout<<"Index type: dual index"<<endl;
      if(barFileR2_name.size()==0)
	{
	  cout<<"please specify the barcode (R2) data input fasta file name.......\n";
	  //printUsage(argc, argv);
	  cout<<"type "<<argv[0]<<" -h \" for usage help\n";
	  exit(-1);
	}
      cout<<"barcode file (Read 2):\""<<barFileR2_name<<"\".\n";
    }
  else
    {
      cout<<"Read type: single index"<<endl;
    }

  //pair end read 
  if(demux||indexFromSeq)
    {//so we have to have sequence file
      if(sequenceFileR1_name.size()==0)
	{
	  cout<<"please specify the sequence data input fasta file name.......\n";
	  //printUsage(argc, argv);
	  cout<<"type "<<argv[0]<<" -h \" for usage help\n";
	  exit(-1);
	}
      cout<<"Sequence data file:"<<sequenceFileR1_name<<endl;
      //pairend yes or no
      if(pairEnd)
	{
	  cout<<"Pair End read: TRUE\n"<<endl;
	  if(sequenceFileR2_name.size()!=0)
	    {
	      cout<<"Sequence data R2 file:"<<sequenceFileR2_name<<endl;
	    }
	}
      else
	{
	  cout<<"Pair End read: FALSE\n"<<endl;
	  if(sequenceFileR2_name.size()!=0)
	    {
	      cout<<"Sequence data R2 file:"<<sequenceFileR2_name<<endl;
	      cout<<"\tspecifed, but will be ignore"<<endl;
	      sequenceFileR2_name="";
	    }
	}
    }
  
  cout<<"Demux:"<< demux<<endl;
  cout<<"index from Sequence:"<<indexFromSeq<<endl;
  
  if(!indexFromSeq)
    {
      if(indexFileR1_name.size()==0)
		{
		  cout<<"***please specify the index data (R1) input fasta file name.......\n";
		  //printUsage(argc, argv);
		  cout<<"type "<<argv[0]<<" -h \" for usage help\n";
		  exit(-1);
		}
      cout<<"index read file (Read 1):\""<<indexFileR1_name<<"\".\n";
      if(dualIndex)
		{
		  if(indexFileR2_name.size()==0)
			{
			  cout<<"****please specify the index data (R2) input fasta file name.......\n";
			  //printUsage(argc, argv);
			  cout<<"type "<<argv[0]<<" -h \" for usage help\n";
			  exit(-1);
			}
		  cout<<"index read file (Read 2):\""<<indexFileR2_name<<"\".\n";
		}
    }
  cout<<"Map type:"<<mapEnd<<endl;
  cout<<"# mismatch:"<<misMatchNum<<endl;
  
  //*********get output files names
  string outFileR1_name;//=(sequenceFile_name+".mapped.fasta"); //the output file for mapped both files
  string outFileR2_name;//=(sequenceFile_name+".mapNone.fasta");//map none
  string outSeqFileR1_name;
  string outSeqFileR2_name;
  size_t pt=0;
  if(indexFromSeq)
    {
      outFileR1_name=sequenceFileR1_name;
	  pt=outFileR1_name.find_last_of('.');
	  
	  if(demux)//now we need to distinguish the index and seq output file
	  {
		outFileR1_name+="_index";
		outSeqFileR1_name=sequenceFileR1_name.substr(0,pt);
		if(pairEnd)
		{
			//work on the second read 
			outSeqFileR2_name=sequenceFileR2_name;
			unsigned pt2=outSeqFileR2_name.find_last_of('.');
			outSeqFileR2_name=outSeqFileR2_name.substr(0,pt2);
		}
	  }  //otherwise, we are fine. let's don't work on outseqfile name, simply leave it empty
	  
	  //get rid of the file suffix
	  if(pt!=string::npos)
	  {
		outFileR1_name=outFileR1_name.substr(0, pt);
	  }
	  
	  if(dualIndex)
	  {
		outFileR2_name=outFileR1_name;
		
		pt=outFileR2_name.rfind("R1");
		if(pt!=string::npos)
		{
			cout<<"To doing replacing, pt:"<<pt<<";outFileR1_name:"<<outFileR1_name<<endl;
			outFileR2_name=outFileR2_name.replace(pt,2,"R2");
			cout<<"after replacing:outFileR1_name:"<<outFileR1_name<<";outFileR2_name:"<<outFileR2_name<<endl;
		}
		else
		{
			outFileR2_name=outFileR1_name+"_R2";
		}
	  }
    }
  else  //from index file;
    {
      outFileR1_name=indexFileR1_name;
	  pt=outFileR1_name.find_last_of('.');
	  if(pt!=string::npos)
		{
			outFileR1_name=outFileR1_name.substr(0, pt);
		}
	  if(dualIndex)
	  {
		outFileR2_name=indexFileR2_name;
		pt=outFileR2_name.find_last_of('.');
		if(pt!=string::npos)
		{
			outFileR2_name=outFileR2_name.substr(0, pt);
		}
	  }
	  
	  if(demux)//now we need to distinguish the index and seq output file
	  {
		//outFileR1_name+="_index";
		pt=sequenceFileR1_name.find_last_of('.');
		outSeqFileR1_name=sequenceFileR1_name.substr(0,pt);
		if(pairEnd)
		{
			//work on the second read 
			outSeqFileR2_name=sequenceFileR2_name;
			unsigned pt2=outSeqFileR2_name.find_last_of('.');
			outSeqFileR2_name=outSeqFileR2_name.substr(0,pt2);
		}
	  }  //otherwise, we are fine. let's don't work on outseqfile name, simply leave it empty
	  
    }

  cout<<"\tThe output file name index R1: \""<<outFileR1_name<<"\"\n"
	  <<"\toutputFile index R2:\""<<outFileR2_name<<"\";\n"
	  <<"\toutputFile seq R1:\""<<outSeqFileR1_name<<"\";\n"
	  <<"\toutputFile seq R2:\""<<outSeqFileR2_name<<"\";\n"
    //<<"\tmapEnd:"<<mapEnd<<"\n"
      <<"\n";
  
  MatchMatrix* mm= &mmf;//ScoreMatrixArr[scoreMatrixIndex];
      //mmf has been done declaration and initialization in score.hpp/cpp
  cout<<"The matching matrix for mm(A,T):"<<mm->GetScore('A','T')<<endl;

  //start doing the reading....assuming all fasta, will do fastq later
  //fasta handler:reading fasta files depending on the input
  cout<<"--------------------\n";
  cout<<"Starting reading the data from input file.....\n";
  vector<SequenceString> vec_index1;
  vector<SequenceString> vec_index2;
  vector<SequenceString> vec_bar_seq1;//this will hold processed sequences on the forward end
  vector<SequenceString> vec_bar_seq2;
  vector<SequenceString> vec_seq1;
  vector<SequenceString> vec_seq2;

  //------all file are in fasta format
  if(indexFromSeq)
    {
		if(sequenceFileR1_name.length()==0)
		{
			cout<<"*****ERROR: request to read index from the sequence files, "
				<<"\tbut the sequence R1 file has not been specified !!!"<<endl;
			exit(-1);
		}
      //reading R1
      cout<<"reading sequence data file (Read 1): "<<ReadFasta(sequenceFileR1_name, vec_seq1)<<"sequences read"<<endl;
      if(demux&&pairEnd)  //in this case, we don't have to get index from it, assuming indexes are in the name xxxxx+xxxxxx
		{
			if(sequenceFileR2_name.length()==0)
			{
				cout<<"*****ERROR: request to do demux for pairEnd data, "
					<<"\tbut the sequence R2 file has not been specified !!!"<<endl;
				exit(-1);
			}
		  //reading R2
		  cout<<"reading sequence data file (Read 2): "<<ReadFasta(sequenceFileR2_name, vec_seq2)<<"sequences read"<<endl;
		}
    }
  else //reading
    {
	  if(indexFileR1_name.length()==0)
		{
			cout<<"*****ERROR: request to read index from the index files, "
				<<"\tbut the index R1 file has not been specified !!!"<<endl;
			exit(-1);
		}
      cout<<"reading index data file (Read 1): "<<ReadFasta(indexFileR1_name, vec_index1)<<"indexes read"<<endl;
      if(dualIndex)
		{
			if(indexFileR2_name.length()==0)
			{
				cout<<"*****ERROR: request to read dual indexes from the index files, "
					<<"\tbut the index R2 file has not been specified !!!"<<endl;
				exit(-1);
			}
			cout<<"reading index data file (Read 2): "<<ReadFasta(indexFileR2_name, vec_index2)<<"indexes read"<<endl;
		}
      if(demux)
		{
		  if(sequenceFileR1_name.length()==0)
			{
				cout<<"*****ERROR: request to demux from the sequence files, "
					<<"\tbut the sequence R1 file has not been specified !!!"<<endl;
				exit(-1);
			}
		  cout<<"reading sequence data file (Read 1): "<<ReadFasta(sequenceFileR1_name, vec_seq1)<<"sequences read"<<endl;
		  if(pairEnd)
			{
				if(sequenceFileR2_name.length()==0)
				{
					cout<<"*****ERROR: request to demux for pairEnd sequence files, "
						<<"\tbut the sequence R2 file has not been specified !!!"<<endl;
					exit(-1);
				}
			  cout<<"reading index data file (Read 2): "<<ReadFasta(sequenceFileR2_name, vec_seq2)<<"sequences read"<<endl;
			}
		}      
    }
	
	
	if(barFileR1_name.length()==0)
	{
		cout<<"*****ERROR: no bar code file R1 has not been specified !!!"<<endl;
		exit(-1);
	}
  cout<<"reading barcode data file  (Read 1): "<<ReadFasta(barFileR1_name, vec_bar_seq1)<<"barcodes read"<<endl;
  if(dualIndex)
    {
		if(barFileR2_name.length()==0)
		{
			cout<<"*****ERROR: request to do dual index, "
				<<"\tbut the bar code R2 file has not been specified !!!"<<endl;
			exit(-1);
		}
      cout<<"reading barcode data file  (Read 2): "<<ReadFasta(barFileR2_name, vec_bar_seq2)<<"barcodes read"<<endl;
    }
  
  //cout<<"reading sequence data file: "<<ReadFasta(sequenceFileR1_name, vec_seq1)<<endl;
  //cout<<"1/1000:"<<vec_seq1.at(0).toString()<<endl;
  cout<<"----------------------\n";
  cout<<"reading and processing sequences for barcodes/demux......"<<endl;

  //now start testing the read indexes from seq
  cout<<"vec_Seq1.size():"<<vec_seq1.size()<<endl;
  if(indexFromSeq){
	  unsigned num=ReadIndexFromSequenceName(vec_seq1, vec_index1, vec_index2, 8, dualIndex);
	  cout<<"vec_Seq1.size():"<<vec_seq1.size()<<endl;
	  cout<<"Total number of sequences and indexes read:"<<num<<endl;
	  unsigned x=7;
	  cout<<"vec_bar_seq1:"<<vec_bar_seq1.at(x).toString()<<endl;
	  if(dualIndex){
		cout<<"vec_bar_seq2:"<<vec_bar_seq2.at(x).toString()<<endl;
	  }
	  cout<<"ve_seq1:"<<vec_seq1.at(x).toString()<<endl;
  }
  else{
	cout<<"*******reading index from index files;"<<endl; 
  }
  //testing the mathcing barcoe scoring system.
	cout<<"TEsting matching score..............."<<endl;
  SequenceString seq1("seq1", "ACTN");
	SequenceString barcode("barcode", "ACTGAN");
  double score=MatchBarcodes(seq1, barcode, mm);
  cout<<"seq1:"<<seq1.toString()<<"\nbarcode:"<<barcode.toString()<<endl;
  cout<<"Matching score is "<<score<<endl;
//now we have everything, we just need to do the job, I mean mapping, here.
  MappingBarcodes(vec_seq1, vec_seq2, vec_index1, vec_index2,vec_bar_seq1, vec_bar_seq2, 
		  misMatchNum, pairEnd, dualIndex, indexFromSeq, mapEnd, demux, mm,
		  outFileR1_name /*index file name*/, outFileR2_name /*index file R2*/,
		  outSeqFileR1_name /*seq file name*/, outSeqFileR2_name /*seq file*/
		  ); 
  /*
  //test the writting the text table
  vector<string> header;
  header.push_back("gene");
  //header.push_back("stats");
  
  vector<double> gene_vec;
  gene_vec.push_back(1.0);
  gene_vec.push_back(2.0);
  
  vector<double> stats_vec;
  stats_vec.push_back(3.0);
  stats_vec.push_back(5.0);
  stats_vec.push_back(7.0);
  
  vector<vector<double> > v_v;
  v_v.push_back(gene_vec);
  v_v.push_back(stats_vec);
  WriteTextTableFile("test.txt", v_v, '\t',
		     true, ofstream::trunc, header);
  */

  cout<<"Done!!!"<<endl;
  //ofs.close();

  //cout<<"Total "<<gene_info.size()<<" genes are processed"<<endl;
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
	case 'f':
	  indexFileR1_name=optarg;
	  break;
	
	case 'e':
	  indexFileR2_name=optarg;
	  break;
	
	case 'b':
	  barFileR1_name=optarg;
	  break;
	
	case 'a':
	  barFileR2_name=optarg;
	  break;
       
	case 's':
	  sequenceFileR1_name=optarg;
	  break;
	case 'q':
	  sequenceFileR2_name=optarg;
	  break;
	case 'm':
	  misMatchNum=atoi(optarg);
	  break;
	  
	case 't':
	  indexFromSeq=true;
	  break;
	  //case 'i':
	  //isotype_flag=true;
	  //break;

	case 'x':
	  demux=true;
	  break;
	case 'd':
	  dualIndex=true;
	  break;
	case 'r':
	  mapEnd=ThreePrime;
	  break;
	case 'p':
	  pairEnd=true;
	  break;
	  /*
	case 'd':
	  temp=atoi(optarg);
	  if(temp==1)
	    {
	      mapEnd=FivePrime;
	    }
	  if(temp==2)
	    {
	      mapEnd=ThreePrime;
	    }
	  if(temp!=1&&temp!=2)
	    {
	      cout<<"ERROR: unknown mapping type, can only be 1 or 2 (for 5prime or 3 prime mapping, respecitively);\n"<<endl;
	      printUsage(argc, argv);
	      exit(-1);
	    }
	  break;

	case 'l':
	  MinimumOverlapLength=atoi(optarg);
	  break;
	case 'e':
	  gapextension=atoi(optarg);
	  if(gapextension>0)
	    gapextension*=-1;
	  gapextensionFlag=true;
	  break;
	case 'g':
	  gapopen=atoi(optarg);
	  if(gapopen>=0)
	    gapopen*=-1;
	 
	  if(!gapextensionFlag)
	    gapextension=gapopen;
	  break;
	case 'x':
	  demux=true;
	  break;
	  */
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
  cout<<argv[0]<<", the program used to demux or runs stats of indexes. It would read in\n"
      <<"barcodes and sequence indexes and then match(instead of align) to find them.\n"
      <<"We assume that the index (read) and barcode (expected) are having identical length.\n"
      <<"Well, actually, the index could be long, but the barcodes are having fixed length.\n"
      <<"we always match/align barcodes again index(as template) and we simply match (compare)\n"
      <<"instead of aligning them. We also allow degeneracy on barcode, but not on indexes.\n"
      <<"Therefore, the match score, mmf, is asymmetric. We also allow to read sequence \n"
      <<"index from the name of sequences instead of sequence read file. When sequences\n"
      <<"are specified, we can do demux. Otherwise, we will only output stats about\n"
      <<"barcoded sequences. We allow single or dual indexes. It is also possible to \n"
      <<"do 5' or 3' aligned on the second read depending on how the input specified.\n";

  cout<<"\tOption string: \"hva:b:f:e:ptdrxs:n:\" \n";
  cout<<"Usage:\n";
  cout<<"\t"<<argv[0]<<" [-s sequence file] -a adaptor file (R1) [-b barcode file (R2)]\n"
      <<"\t\t [-f  sequence index R1 file] [-e sequence index R2 file]\n"
      <<"\t\t [-p (pair end read,true or false)] [-t (read indexes from the name of sequences] \n"
      <<"\t\t [-x (demux or not)] [-r (reverse implement r2 sequence index)]\n"
      <<"\t\t [-q sequence file r2, optional]\n"
      <<"\t\t [-m maximal num of mismatch allowed] [-d (dual index)]\n"
       <<"\t\t [-h] [-v]\n\n"
    ;
  cout<<"\t-a: the input file barcode name (r2 for dual index input), \n"
      <<"\t\tassuming has identical name and order\n";
  cout<<"\t-b: the input barcode file name (r1 or for single index input\n"
      <<"\t a and b list the expected barcode to demux the input\n";
  cout<<"\t-f: the input file name read 1 index\n"
      <<"\t-e: the input file name read 2 index\n" 
      <<"\tf and e contain the read index from sequencing to be demux'ed\n";
  cout<<"\t-p: paired end read, boolean, false by default\n" ;
  cout<<"\t-t: read indexs from the names of each sequence and assuming \n"
      <<"\t\tthe illumina style:\n"
      <<"\t\t\t  xxxxxx:xxxxx:xxxx:actgkd+acdfad\n" 
      <<"\t\t\t\tseparated by ':' and the last field contains the barcode (pair)\n"
      <<"\t\t\t\tfalse by default\n";
  cout<<"\t-x: demux (true or false), meaning to write out the sequences\n"
      <<" \t\tor simple output the stats. False by default";
  cout<<"\t-d: dual index (true or false), true by default\n"
      <<" \t\tor simple output the stats";
  cout<<"\t-r: reverse complement the read2 barcode or not to do demux\n"
      <<"\t\tonly work for Read2, to reverse complement the sequence index 2,\n"
      <<"\t\twe never reverse complement the barcode.False(FivePrime) by default\n"; 
  cout<<"\t-q: sequence data file (R2), optional, could be\n"
      <<"\t\tspecified or not. if dumux, then needed to check for this.\n";
  cout<<"\t-s: sequece data file, contains the real sequence data, \n"
      <<"\t\thaving the same order and name.\n";

  cout<<"\t-m: maximal number of mismatches allowed, \n"
      <<"\t\tallow degenerate barcode, meaning the barcode (expected)\n"
      <<"\t\t could be degenerated, but not the other way around.\n"
      <<"\t\t barcode (expected)-->row, index (sequence read) ---->column\n"
      <<"\t\t also this number is the total number of mismatches allow \n"
      <<"\t\t for both barcodes if there are dual indexeds. 1 by default\n";

  cout<<"\t-v and -h: version and help\n";
  cout<<"Note: 1)when do comparison, always use barcode(expected) to \n"
      <<"\talign index (data), in terms of alignment, barcode is sujbect,\n"
      <<"\t index is pattern\n"
      <<"\t             index (pattern)\n"
      <<"\t             |||||\n"
      <<"\t             barcode (subject)\n"
      <<"\t2)allow degeneracy on barcode, but less so on index\n"
      <<"\t3)it is possible to have index longer than barcode, but not vice versa\n"
      <<"\t4)barcode and index must start on 5' position 1 (position 0 in c)\n"
      <<"\t5)we don't do alignment in comparison, but string comparison, \n"
      <<"\t\tsee comparison matrix, match.matrix.feng (mmf)\n"
      <<"\t6)the match score could be fraction number, because of degeneracy\n"
      <<"\t7)all the files should have the same name or order in order to reference them\n"
    
      <<"\t\t\t *************@4/2/2014 by Feng\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********updated by Feng @ BU 2018\n"; 

}

