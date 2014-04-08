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
#include "OverlapAlignment.hpp"
#include "FastaHandler.hpp"
#include "SequenceHandler.hpp"
using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);
static int lookUpScoreMatrix(const string* _scoreMatrixNameArray,const int& len, const string& scoreMatrixName);

//all file are in fasta format
static string adaptorFile_name("454_adaptors.fa"); //the input file for adaptor
static string barcodeFile_name("454_barcodes.fa");//input ifle for barcode
static string forwardFile_name("primer_forward.fa");//input ifle for forward constant 
static string reverseFile_name("primer_reverseUPM.fa");//input ifle for reverse UPM
static string sequenceFile_name;//input file for sequence data

static double MismatchRateThreshold=0.8; //not too many mismatch 
static unsigned int MinimumOverlapLength=10;//not too short

//how far we allow the alignment to be away from the ends. can not be too far, since they are supposed to be aligned on the ends.
unsigned int OffsetForward=10;//###10 might too big????
unsigned int OffsetReverse=10;//

static string scoreMatrixName="nuc44"; //name as an input for specifing the name of the score matrix. for example nuc44

/*static enum SeqType
  {
    AminoAcid,
    Nucleotide
  } ; //by default
*/
static string supportedScoreMatrixNameArr[]={"nuc44","blosum50", "tsm1", "tsm2"};

static ScoreMatrix* ScoreMatrixArr[]={&nuc44, &blosum50, &tsm1, &tsm2};

static double scale=1; //this is the one on top of matrix, the programe will run score and use the 
//matrix specified scale first and then apply the scale set by this one.

static double gapopen=-8;
static double gapextension=-8;
static bool gapextensionFlag=false;

static int trim=1;

int main(int argc, char* argv[])
{
  //for parsing commandline arguement

  const char *opts = "hva:b:f:r:g:e:t:m:s:k:n:p:q:l:";
  //a: the input file name adaptor
  //b: the input file name barcode
  //f: the input file name forward primer
  //r: the input file name reverse primer

  //g: gap open
  //e: gap extension
  //m: score matrix
  //s: sequece data file
  //t: get trimmed data file (1) or no trimmed data (0)
  //k: scale factor

  //n: mismatch ratio threshold 
  //p: offset in for the forward
  //q: offset for the reverse end
  //l: minimum overlap length
  parseArguments(argc, argv, opts);
  
  
  if(sequenceFile_name.size()==0)
    {
      cout<<"please specify the sequece data input fasta file name.......\n";
      //printUsage(argc, argv);
      cout<<"type \"./ngsmapping -h \" for usage help\n";
      exit(-1);
    }
  
  //*********get output files
  string outputFileBoth_name=(sequenceFile_name+".mapBoth.fasta"); //the output file for mapped both files
  string outputFileForward_name=(sequenceFile_name+".mapForward.fasta");//map forward
  string outputFileReverse_name=(sequenceFile_name+".mapReverse.fasta");//map reverse
  string outputFileNone_name=(sequenceFile_name+".mapNone.fasta");//map none

  string outputStatFileBoth_name=(sequenceFile_name+".mapBothStat.fasta"); //the output file for mapped both files
  string outputStatFileForward_name=(sequenceFile_name+".mapForwardStat.fasta");//map forward
  string outputStatFileReverse_name(sequenceFile_name+".mapReverseStat.fasta");//map reverse
  string outputStatFileNone_name(sequenceFile_name+".mapNoneStat.fasta");//map none

  string outputTrimmedFile_name(sequenceFile_name+".trim.fasta");

  cout<<"***Input parameters Summary:\n";
  cout<<"\tAdaptor file name:\""<<adaptorFile_name<<"\".\n";
  cout<<"\tBarcode file name:\""<<barcodeFile_name<<"\".\n";
  cout<<"\tforward primer file name:\""<<forwardFile_name<<"\".\n";
  cout<<"\treverse primer file name:\""<<reverseFile_name<<"\".\n";
  cout<<"\tsequence data file name:\""<<sequenceFile_name<<"\".\n";
  
  cout<<"\tThe output file name : \""<<outputFileBoth_name<<"\",\""<< outputFileForward_name<<"\",\""<< outputFileReverse_name<<"\",\""
      <<outputFileNone_name<<"\";\n"
      <<"\t\t"<<outputStatFileBoth_name<<"\",\""<< outputStatFileForward_name<<"\",\""<< outputStatFileReverse_name<<"\",\""
      <<outputStatFileNone_name<<"\";\n\n";
  cout<<"\t\t\""<<outputTrimmedFile_name<<"\"\n";
  
  //look up the score matix
  int scoreMatrixIndex=lookUpScoreMatrix(supportedScoreMatrixNameArr, sizeof(supportedScoreMatrixNameArr)/sizeof(supportedScoreMatrixNameArr[0]),scoreMatrixName);
  if(scoreMatrixIndex==-1)
    {
      cout<<"\tscore matrix specified by input was not found. Using the default scorematrix.\n";
      scoreMatrixName="nuc44";
      
      scoreMatrixIndex=lookUpScoreMatrix(supportedScoreMatrixNameArr, sizeof(supportedScoreMatrixNameArr)/sizeof(supportedScoreMatrixNameArr[0]), scoreMatrixName);
    }
  cout<<"\tscore matrix:"<<scoreMatrixName<<"\n"
  
      <<"\tscale matrix scale:"<<scale<<"\n"
  
      <<"\tgap open penalty:"<<gapopen<<"\n"
      <<"\tgap extension penalty:"<<gapextension<<"\n"
      
      <<"\toffset on forward end:"<<OffsetForward<<"\n"
      <<"\toffset on reverse end:"<<OffsetReverse<<"\n"
      <<"\tmismatch rate threshold:"<<MismatchRateThreshold<<"\n"
      <<"\tminimum overlap length:"<<MinimumOverlapLength<<"\n";
  cout<<"  ****************\n";

  cout<<"Testing ScoreMatrix:\n";
  cout<<"\tthe index: "<<scoreMatrixIndex<<endl;
  ScoreMatrix* sm= ScoreMatrixArr[scoreMatrixIndex];
  char c1='C', c2='C';

  cout<<"\t("<<c1<<","<<c2<<")="<<sm->GetScore(c1,c2)<<endl;
  
  
  //testing fasta handler
  vector<SequenceString> vec_seq;
  vector<SequenceString> vec_forward_seq;//this will hold processed sequences on the forward end
  vector<SequenceString> vec_reverse_seq;//this will hold processed sequences on the reverse end
  
  cout<<"reading sequence data file: "<<ReadFasta(sequenceFile_name, vec_seq)<<endl;
  
  
  cout<<"1/1000:"<<vec_seq.at(0).toString()<<endl;
  //cout<<"1000/1000;"<<vec_seq.at(999).toString()<<endl;

  cout<<"reading and processing adaptor sequences for mapping......"<<endl;
  ProcessAdaptorSequences(adaptorFile_name, barcodeFile_name, forwardFile_name, reverseFile_name, vec_forward_seq, vec_reverse_seq);

  
  
  WriteFasta("forward_adaptor_sequnces_set.fa",vec_forward_seq);
  WriteFasta("reverse_adaptor_sequnces_set.fa",vec_reverse_seq);  

//now we have everything, we just need to do the job, I mean mapping, here.
  MappingAdaptors(vec_forward_seq, vec_reverse_seq, vec_seq, 
		  sm, gapopen, gapextension, trim,
		  MismatchRateThreshold, MinimumOverlapLength, OffsetForward, OffsetReverse, 
		  outputFileBoth_name, outputFileForward_name, outputFileReverse_name,outputFileNone_name
		  ); 
  

  cout<<"Done!!!"<<endl;
  //ofs.close();

  //cout<<"Total "<<gene_info.size()<<" genes are processed"<<endl;
  cout<<"Thanks for using our program and have a nice day!!"<<endl;
  return 0;
}


static void parseArguments(int argc, char **argv, const char *opts)
{
  int optchar;

  while ((optchar = getopt(argc, argv, opts)) != -1)
    {
      switch(char(optchar))
	{
	case 'a':
	  adaptorFile_name=optarg;
	  break;

	case 'b':
	  barcodeFile_name=optarg;
	  break;  
	  
	case 'f':
	  forwardFile_name=optarg;
	break;
	
	case 'r':
	  reverseFile_name=optarg;
	break;
	
	case 's':
	  sequenceFile_name=optarg;
	break;
	
	case 'm':
	  scoreMatrixName=optarg;
	  break;
	  
	case 't':
	  trim=atoi(optarg);
	  break;
	case 'k':
	  scale=atof(optarg);
	  break;
	case 'n':
	  MismatchRateThreshold=atof(optarg);
	  break;
	case 'p':
	  OffsetForward=atoi(optarg);
	  break;
	case 'q':
	  OffsetReverse=atoi(optarg);
	  break;
	case 'l':
	  MinimumOverlapLength=atoi(optarg);
	  break;
	case 'e':
	  gapextension=atoi(optarg);
	  gapextensionFlag=true;
	  break;
	case 'g':
	  gapopen=atoi(optarg);
	  if(!gapextensionFlag)
	    gapextension=gapopen;
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
  //"hva:b:f:r:g:e:t:m:s:k:"
  cout<<"Usage:\n";
  cout<<"\t"<<argv[0]<<" -s second sequence file [-a adaptor file] [-b barcode file]  \n"
      <<"\t [-t reverse primer filename] [-f forward primer filename]\n"
      <<"\t  [-m score matrix] [-t trim ] [-l MinimumOverlapLength]\n"
      <<"\t [-k scale] [-g gapopen panelty] [-e gap extension]\n"
      <<"\t [-n Mismatch rate threshold] [-p offset on forward end] [-q offset on reverse end]\n"     
      <<"\t [-h] [-v]\n\n"
    ;

  cout<<"\t\t-s filename -- the sequence fasta data filename \n"
      <<"\n";

  cout<<"\t\t-a filename -- the adaptor fasta filenames\n"
      <<"\n";

  cout<<"\t\t-b filename -- the barcode fasta filenames \n"
      <<"\n";

  cout<<"\t\t-f filename -- the forward primer fasta filenames \n"
      <<"\n";

  cout<<"\t\t-r filename -- the reverse primer fasta filenames \n"
      <<"\n";
  
  cout<<"\t\t-m scorematrix -- the socre matrix name used for the alignment, \n"
      <<"\t\t\tonly support nuc44. nuc44 by default for nucleotide\n"
      <<"\n";

  cout<<"\t\t-t # -- the number 0 indicate no trimmed data to be saved; \n"
      <<"\t\t\t  All other number means to save the trimmed data to file\n\n";

  cout<<"\t\t-k scale -- the scale factor used to set the returned score to the correct unit,\n"
      <<"\t\t\t 1 by default. The programe first uses the scale factor coming with matrix\n"
      <<"\t\t\t  to the return score and the the scale set by this option\n\n";

  cout<<"\t\t-g gapopen -- the gap open value\n"
      <<"\n";

  cout<<"\t\t-e gapextension -- the gap extension value\n"
      <<"\n";
  cout<<"\t\t-l minimum overlap length\n"
      <<"\n"; 
  cout<<"\t\t-n mismatch rate threshold\n"
      <<"\n"; 
  cout<<"\t\t-p offset on the forward end\n"
      <<"\n"; 
  cout<<"\t\t-q offset on the reverse end\n"
      <<"\n"; 
  
  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to map the adaptor to the sequence\n"
    <<"\t\t\t @4/2/2014\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********by Feng @ BU 2013\n"; 


  exit(-1);
}
//return the index that could be used to pick up the ScoreMatrix object from ScoreMatrixArr
//-1 for can not find
static int lookUpScoreMatrix(const string* _scoreMatrixNameArray,const int& len, const string& scoreMatrixName)
{
  //int index=-1;
  //int size=sizeof(_scoreMatrixNameArray)/sizeof(_scoreMatrixNameArray[0]);
  for(int i=0;i<len;i++)
    {
      if(scoreMatrixName.compare(_scoreMatrixNameArray[i])==0)
	{
	  return i;
	}
    }
  return -1;
}

