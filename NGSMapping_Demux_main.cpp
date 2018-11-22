#include <iostream>
#include <fstream>
#include <vector>
#include "string_ext.hpp"

//from now on you need to specify the c libraries explicitly
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "score.hpp"
#include "SequenceString.hpp"
//#include "OverlapAlignment.hpp"
#include "FastaHandler.hpp"
#include "SequenceHandlerIsotype.hpp"
#include "SequenceHandlerCommon.hpp"
using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);
//static int lookUpScoreMatrix(const string* _scoreMatrixNameArray,const int& len, const string& scoreMatrixName);

//all file are in fasta format
static string barcodeFile_name("HumanIGConstant_Lib.fas");//default input file for single index, R1

static string sequenceFile_name;//input file for sequence data

static unsigned  misMatchNum=1; //not too many mismatch 
//static unsigned int MinimumOverlapLength=10;//not too short

static mapType mapEnd=FivePrime;

static bool demux=false;//write stats only, no demux

//note for match  matrix
//the match matrix mmf has been declare and initialized in 
//score.hpp and score.cpp, respectively. we can use it 
//directly
//   MatchMatrix mmf; 


int main(int argc, char* argv[])
{
  //for parsing commandline arguement
  const char *opts = "hva:b:f:e:ptdr:s:n:";
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
  //d: demux, meaning to write out the sequences or simple output the stats
  //r: reverse complement the read2 barcode or not to do demux

  //s: sequece data file, contains the real sequence data, having the same order and name.

  //n: mismatch, allow degenerate barcode, meaning the barcode (expected) could be degenerated, but not the other way around
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
  parseArguments(argc, argv, opts);
    
  if(sequenceFile_name.size()==0)
    {
      cout<<"please specify the sequece data input fasta file name.......\n";
      //printUsage(argc, argv);
      cout<<"type \"./ngsmapping -h \" for usage help\n";
      exit(-1);
    }
  
  //*********get output files
  string outputFileMap_name=(sequenceFile_name+".mapped.fasta"); //the output file for mapped both files
  string outputFileNone_name=(sequenceFile_name+".mapNone.fasta");//map none

  cout<<"***Input parameters Summary:\n";
  cout<<"\tIsotype file name:\""<<isotypeFile_name<<"\".\n";
  cout<<"\tsequence data file name:\""<<sequenceFile_name<<"\".\n";
  cout<<"\tThe output file name : \""<<outputFileMap_name<<"\",\""
      <<outputFileNone_name<<"\";\n"
      <<"\tmapEnd:"<<mapEnd<<"\n"
      <<"\n";
    
  cout<<"\tgap open penalty:"<<gapopen<<"\n"
      <<"\tgap extension penalty:"<<gapextension<<"\n"
      
      <<"\toffset on forward end:"<<Offset<<"\n"

      <<"\tmatch rate threshold:"<<matchRateThreshold<<"\n"
      <<"\tminimum overlap length:"<<MinimumOverlapLength<<"\n"
      <<"\tdemux output:"<<demux<<"\n";
  cout<<"  ****************\n";

  MatchMatrix* mm= &mmf;//ScoreMatrixArr[scoreMatrixIndex];
      //mmf has been done declaration and initialization in score.hpp/cpp


  //fasta handler:reading fasta
  vector<SequenceString> vec_seq;
  vector<SequenceString> vec_Isotype_seq;//this will hold processed sequences on the forward end
  
  cout<<"reading sequence data file: "<<ReadFasta(sequenceFile_name, vec_seq)<<endl;
  cout<<"1/1000:"<<vec_seq.at(0).toString()<<endl;
  
  cout<<"reading and processing isotype sequences for mapping......"<<endl;

  //SetUpByIsotypeOutputFlag(true);
  cout<<"Reading constant library file:"<<ReadFasta(isotypeFile_name, vec_Isotype_seq)<<endl;


  //reverse the Ig constant region sequence
  /*for(unsigned int i=0;i<vec_forward_seq.size();i++)
    {
      vec_forward_seq[i]=ReverseComplement(vec_forward_seq.at(i));
    }
  */
 
//now we have everything, we just need to do the job, I mean mapping, here.
  MappingIsotypes(vec_seq, vec_Isotype_seq, 
		  mapEnd,sm, gapopen, gapextension,
		  matchRateThreshold, MinimumOverlapLength, 
		  Offset, 
		  outputFileMap_name, outputFileNone_name,demux
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
  int temp;
  while ((optchar = getopt(argc, argv, opts)) != -1)
    {
      switch(char(optchar))
	{	  
	case 'f':
	  isotypeFile_name=optarg;
	break;
	
	//case 'r':
	//  reverseFile_name=optarg;
	//break;
	
	case 's':
	  sequenceFile_name=optarg;
	break;
	
	case 'm':
	  scoreMatrixName=optarg;
	  break;
	  
	  //case 't':
	  //trim=atoi(optarg);
	  //break;
	  //case 'i':
	  //isotype_flag=true;
	  //break;

	case 'k':
	  scale=atof(optarg);
	  break;
	case 'n':
	  matchRateThreshold=atof(optarg);
	  break;
	case 'p':
	  Offset=atoi(optarg);
	  break;

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
  cout<<"argv[0], the program used to identify the isotypes of the sequences and also\n"
      <<"\tdemutiplex them. This will work mainly on the new design of IgSeq data, which\n"
      <<"\tare pair-end and have We try to follow the style of cutadapt. Basically, we will\n"
      <<"\tonly ask for the isotype sequence and then map them to one side \n"
      <<"\tof the sequences, which could be either 5' or 3'\n";
  cout<<"\tOption string: \"hva:b:f:e:ptdr:s:n:\" \n";
  cout<<"Usage:\n";
  cout<<"\t"<<argv[0]<<" [-s sequence file] -a adaptor file [-b barcode file]  \n"
      <<" [-f  sequence file name] [-e "
      <<" [-d mapping type]\n"
      <<"\t  [-m score matrix] [-l MinimumOverlapLength]\n"
      <<"\t [-k scale] [-g gapopen panelty] [-e gap extension]\n"
    //<<"\t [-i]\n"
      <<"\t [-n match rate threshold] [-p offset on forward end] \n"     
      <<"\t [-x] [-h] [-v]\n\n"
    ;

  cout<<"\t\t-s filename -- the sequence fasta data filename \n"
      <<"\n";

  cout<<"\t\t-f filename -- the isotype sequences fasta filename \n"
      <<"\t\t\tif 3' is chosen, the isotype sequences will be reverse\n"
      <<"\t\t\tcomplemented before mapping.\n"
      <<"\n";

  cout<<"\t\t-d # -- the mapping type, 5' (by default) or 3' mapping.\n"
      <<"\t\t\t 1 for 5' and 2 for 3' mapping\n"
      <<"\n";
  
  cout<<"\t\t-m scorematrix -- the socre matrix name used for the alignment, \n"
      <<"\t\t\tonly support nuc44. nuc44 by default for nucleotide\n"
      <<"\n";

  //cout<<"\t\t-t # -- the number 0 indicate no trimmed data to be saved; \n"
  //    <<"\t\t\t  All other number means to save the trimmed data to file\n\n";

  cout<<"\t\t-k scale -- the scale factor used to set the returned score to the correct unit,\n"
      <<"\t\t\t 1 by default. The programe first uses the scale factor coming with matrix\n"
      <<"\t\t\t  to the return score and the the scale set by this option\n\n";

  cout<<"\t\t-g gapopen -- the gap open value. will be turned into negative if not\n"
      <<"\t\t\t 15 by default\n"
      <<"\n";

  //cout<<"\t\t-i -- set to output data by isotypes\n"
  //    <<"\n";

  cout<<"\t\t-e gapextension -- the gap extension value. will be turned into negative if not\n"
      <<"\t\t\t10 by default\n"
      <<"\n";
  cout<<"\t\t-l minimum overlap length, 10 by default\n "
      <<"\t\t\tsee note below in -n match rate threshold!!!\n"
      <<"\n"; 
  cout<<"\t\t-n match rate threshold,0.75 by default. This is not\n"
      <<"\t\t\tthe error rate, but the match rate. It is empirically 0.3\n"
      <<"\t\t\tminimum overlap length and match rate are not conflicting with each other,\n"
      <<"\t\t\tsince partial matching is allowed. We specify the minimum overlap in case\n"
      <<"\t\t\tit happens that a very small overlap with higher matching rate. For example\n"
      <<"\t\t\tit could purely by chance be that 4 matches partial show in the matching.\n"
      <<"\n"; 
  cout<<"\t\t-p offset on the sequence ends, 15 by default\n"
      <<"\n"; 
  //cout<<"\t\t-q offset on the reverse end\n"
  //    <<"\n"; 
  cout<<"\t\t-x this is one to indicating whether to write demultiplexed output\n"
      <<"\t\t\t note: when this is set, the output will write to have only\n"
      <<"\t\t\t sequences arranged into different isotyped, but not mapped\n"
      <<"\n";
  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to map the adaptor to the sequence\n"
      <<"\t\t\tthere are 3 different outputs, two of which will always be\n"
      <<"\t\t\twritten. The second one will only be written as requested by -x.\n"
      <<"\t\t\tFirst one is the mapped file (and unmapped too), containing\n"
      <<"\t\t\tthe sequence and aligned isotype sequence. The second is the\n"
      <<"\t\t\tstats file, which contains the information about the alignment\n"
      <<"\t\t\tincluding match rate and overlaping length, starting and ending\n"
      <<"\t\t\tposition of the sequence and isotype as well as the length of\n"
      <<"\t\t\tthe sequence.This stats file is included to make it easy to do\n"
      <<"\t\t\tsummary\n"
      <<"\t\t\t @4/2/2014 by Feng\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********updated by Feng @ BU 2018\n"; 


  //exit(-1);
}

