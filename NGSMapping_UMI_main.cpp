#include <iostream>
#include <fstream>
#include <vector>
#include "Accessory/string_ext.hpp"

//from now on you need to specify the c libraries explicitly
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "score.hpp"
#include "Accessory/SequenceString.hpp"
#include "OverlapAlignment.hpp"
#include "Accessory/FastaHandler.hpp"
#include "Accessory/FastqHandler.hpp"
#include "Accessory/FileHandler.hpp"
#include "SequenceHandler.hpp"
#include "SequenceHandler_umi.hpp"
using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);
static int lookUpScoreMatrix(const string* _scoreMatrixNameArray,const int& len, const string& scoreMatrixName);

//all file are in fasta format
static string R1File_name; //the input file for R1
static string R2File_name;//input ifle for R2
static string OR1File_name;//output ifle for R1 
static string OR2File_name;//output ifle for R2

static string OR1FileNU_name;//output file for sequence data, R1 containing no UMI
static string OR2FileNU_name;//output file for sequence data, R2 containing no UMI

static double Mismatches=0; //not too many mismatch 
//static unsigned int MinimumOverlapLength=10;//not too short, this only important for "anchors", if there are one.
static bool trim=false;
static bool extract=false;

static string umi_pattern;
//how far we allow the alignment to be away from the ends. can not be too far, since they are supposed to be aligned on the ends.
unsigned int OffsetForward=0;//###10 might too big????
//unsigned int OffsetReverse=10;//

static string scoreMatrixName="nuc44DM1"; //name as an input for specifing the name of the score matrix. for example nuc44
					//this is hardcoded to be nuc44DM1, since we use this matrix to calculate the error/mismathes in the handler code;
static string supportedScoreMatrixNameArr[]={/*"nuc44","blosum50", "tsm1", "tsm2", */"nuc44DM1"};
static ScoreMatrix* ScoreMatrixArr[]={/*&nuc44, &blosum50, &tsm1, &tsm2, */&nuc44DM1};

static double scale=1; //this is the one on top of matrix, the programe will run score and use the 
//matrix specified scale first and then apply the scale set by this one.

static double gapopen=-8;
static double gapextension=-8;
static bool gapextensionFlag=false;

int main(int argc, char* argv[])
{
  //for parsing commandline arguement

  const char *opts = "hvp:s:r:o:q:b:d:g:e:m:txk:n:f:";
  //p: pattern for UMI
  //s: sequece data file R1
  //r: sequence data file R2
  
  //o: ouput data file R1 with umi found and processed 
  //q: ouput data file R2 with umi found and processed
  //
  //b: output data file R1 with seqeunce no umi found  
  //d: output data file R2 with sequence no umi found 
  
  //g: gap open
  //e: gap extension
  //m: score matrix
   
  //t: trim the sequence keep the barcode in the sequence. by default we trim it off and set the umi to the sequence name.
  //x: extract the umi and added to the sequence name following the umi_tools style
  //k: scale factor

  //n: # of mismatch allowd 
  //f: biggest offset/shift allowed for the umi (reasonable only when there is an anchors

  parseArguments(argc, argv, opts);
  //printUsage(argc, argv);
  if(!trim&&!extract)
  {
	  cout<<"Please specify at least either to trim or to extract\n"
		  <<"\tOtherwise no output at all\n"<<endl;
	  cout<<"type \""<<argv[0]<<" -h \" for usage help\n";
      exit(-1);
  }	  
  if(R1File_name.size()==0)
    {
      cout<<"please specify the sequece data input fasta file name.......\n";
      //printUsage(argc, argv);
      cout<<"type \""<<argv[0]<<" -h \" for usage help\n";
      exit(-1);
    }
  if(umi_pattern.size()==0)
    {
      cout<<"please specify the umi pattern.......\n";
      //printUsage(argc, argv);
      cout<<"type \""<<argv[0]<<" -h \" for usage help\n";
	  printUsage(argc, argv);
      exit(-1);
    }
	
  // ++++ check the file names
  FileType ft1, ft2;
  ft2=UNKNOWN;
  
  if(exist(R1File_name.c_str())&&is_file(R1File_name.c_str()))
  {
	  ft1=getFileType(R1File_name, false);
	  //cout<<"the input R1 is .......\n";
      //printUsage(argc, argv);
  }
  else
  {
	  cout<<"The R1 input file either is not a file or doesn't exist, please check.\n";
	  exit(-1);
  }
  
  if(OR1File_name.size()==0)
  {
	  OR1File_name.append(R1File_name);
	  switch(ft1)
	  {
		  case GZ_FASTQ:
		  case FASTQ:
			OR1File_name.append("out.fastq");
			break;
		  case GZ_FASTA:
		  case FASTA:
			OR1File_name.append("out.fasta");
			break;
		  default:
			cout<<"the R1 input file type is not supported, please input as either fastq/fasta (gz'ed)\n";
			exit(-1);
			break;
	  }
	  
  }
  if(OR1FileNU_name.size()==0)
  {
	  OR1FileNU_name.append(R1File_name);
	  switch(ft1)
	  {
		  case GZ_FASTQ:
		  case FASTQ:
			OR1FileNU_name.append("outNU.fastq");
			break;
		  case GZ_FASTA:
		  case FASTA:
			OR1FileNU_name.append("outNU.fasta");
			break;
		  default:
			cout<<"the R1 input file type is not supported, please input as either fastq/fasta (gz'ed)\n";
			exit(-1);
			break;
	  }
  }
  //now we also need to check for the  R2 files
  if(R2File_name.size()>0) //there is R2 paired read
  {
	  if(R1File_name.compare(R2File_name)==0)
	  {
		  cerr<<"ERROR: the input R1 and R2 files are the same. please double check\n"<<endl;
		  cerr<<"please type \""<<argv[0]<<" -h \" for help "<<endl; 
		  exit(-1);
	  }
	  //check for the type 
	  ft2=getFileType(R2File_name.c_str(), false);
			
	  if(OR2File_name.size()==0)
	  {
		  OR2File_name.append(R2File_name);
		  switch(ft2)
		  {
			  case GZ_FASTQ:
			  case FASTQ:
				OR2File_name.append(".out.fastq");
				break;
			  case GZ_FASTA:
			  case FASTA:
				OR2File_name.append(".out.fasta");
				break;
			  default:
				cout<<"the R2 input file type is not supported, please input as either fastq/fasta (gz'ed)\n";
				exit(-1);
				break;
		  }
		  
	  }
	  if(OR2FileNU_name.size()==0)
	  {
		  OR2FileNU_name.append(R2File_name);
		  switch(ft2)
		  {
			  case GZ_FASTQ:
			  case FASTQ:
				OR2FileNU_name.append("outNU.fastq");
				break;
			  case GZ_FASTA:
			  case FASTA:
				OR2FileNU_name.append("outNU.fasta");
				break;
			  default:
				cout<<"the R2 input file type is not supported, please input as either fastq/fasta (gz'ed)\n";
				exit(-1);
				break;
		  }
	  }
  }

//now start display parameters
	cout<<"Calling command:\n"
		<<"\t\t ";
	for(int i=0;i<argc;i++)
		cout<<argv[i]<<" ";
	cout<<endl;
  cout<<"***Input parameters Summary:\n";
  cout<<"\tSequence file name (R1):\""<<R1File_name<<"\".\n";
  cout<<"\tSequence output file name (R1):\""<<OR1File_name<<"\".\n";
  cout<<"\tSequence output file name (R1, no UMI):\""<<OR1FileNU_name<<"\".\n";
  
  if(R2File_name.size()>0)
  {
	cout<<"\tSequence file name (R2):\""<<R2File_name<<"\".\n";
	cout<<"\tSequence output file name (R2):\""<<OR2File_name<<"\".\n";
	cout<<"\tSequence output file name (R2, no UMI):\""<<OR2FileNU_name<<"\".\n";
  }
  cout<<"\tUMI pattern:\""<<umi_pattern<<"\".\n";
  //cout<<"\ttrim the sequence:\""<<trim<<"\".\n";
  cout<<"\textract umi:"<<extract<<".\n";
  cout<<"\toffset umi:"<<OffsetForward<<".\n";
  cout<<"\tmismatch allowd:"<<Mismatches<<"\n";
  if(trim)
    cout<<"\ttrimmed data to be written:TRUE"<<endl;
  else
    cout<<"\ttrimmed data to be written:FALSE"<<endl;

  //look up the score matix  */
  int scoreMatrixIndex=lookUpScoreMatrix(supportedScoreMatrixNameArr, sizeof(supportedScoreMatrixNameArr)/sizeof(supportedScoreMatrixNameArr[0]),scoreMatrixName);
  if(scoreMatrixIndex==-1)
    {
      cout<<"\tscore matrix specified by input was not found. Using the default scorematrix.\n";
      scoreMatrixName=supportedScoreMatrixNameArr[0];
      
      scoreMatrixIndex=lookUpScoreMatrix(supportedScoreMatrixNameArr, sizeof(supportedScoreMatrixNameArr)/sizeof(supportedScoreMatrixNameArr[0]), scoreMatrixName);
    }
  cout<<"\tscore matrix:"<<scoreMatrixName<<"\n"
      <<"\tscale matrix scale:"<<scale<<"\n"
      <<"\tgap open penalty:"<<gapopen<<"\n"
      <<"\tgap extension penalty:"<<gapextension<<"\n";
    
  cout<<"  ****************\n";

  cout<<"Testing ScoreMatrix:\n";
  cout<<"\tthe index: "<<scoreMatrixIndex<<endl;
  ScoreMatrix* sm= ScoreMatrixArr[scoreMatrixIndex];
  char c1='C', c2='N';
  cout<<"\t("<<c1<<","<<c2<<")="<<sm->GetScore(c1,c2)<<endl;
  unsigned umi_length=0;
  unsigned* umi_positions=new unsigned[umi_pattern.size()]; 
  unsigned anchor_length=0;
  unsigned* anchor_positions=new unsigned[umi_pattern.size()]; 
  bool valid_umi=parseUmiPattern(umi_pattern, umi_positions, umi_length, anchor_positions, anchor_length);
  if(!valid_umi)
  {
	 cout<<"Error: the umi pattern string is not valid"<<endl;
	 printUsage(argc, argv);
	 exit(-1);
  }
  cout<<"parsing results:umi_length,"<<umi_length<<"; anchor length:"<<anchor_length<<endl;
  cout<<"\tumi pos:";
  for(unsigned i=0;i<umi_length;i++)
  {
	cout<<umi_positions[i]<<",";
  }
  cout<<endl;
  cout<<"\tumi cha:";
  for(unsigned i=0;i<umi_length;i++)
  {
	cout<<umi_pattern.at(umi_positions[i]);
  }
  cout<<endl;
  cout<<"\tanchor pos:";
  for(unsigned i=0;i<anchor_length;i++)
  {
	cout<<anchor_positions[i]<<",";
  }
  cout<<endl;
  cout<<"\tanchor cha:";
  for(unsigned i=0;i<anchor_length;i++)
  {
	cout<<umi_pattern.at(anchor_positions[i]);
  }
  cout<<endl;
  
  //
  cout<<"\t**********************"<<endl;
  if(anchor_length>0)
	cout<<"anchor string for alignment:"<<getAnchor(umi_pattern, anchor_positions, anchor_length)<<endl;
  else
	cout<<"anchor string is empty"<<endl;
  //cout<<"umi string for insertion:"<<getUMI(umi_pattern, umi_positions, umi_length)<<endl;
  
  
  //testing fasta handler
  vector<SequenceString> vec_seq;
  vector<string> quality;
  if(ft1==GZ_FASTA||ft1==FASTA)
	cout<<"======>reading sequence data file: "<<ReadFasta(R1File_name, vec_seq)<<endl;
  if(ft1==GZ_FASTQ||ft1==FASTQ)
	cout<<"======>reading sequence data file: "<<ReadFastq(R1File_name, vec_seq, quality)<<endl;
  if(ft2==GZ_FASTA||ft2==FASTA) 
	cout<<"======>reading sequence data file: "<<ReadFasta(R2File_name, vec_seq)<<endl;
  if(ft2==GZ_FASTQ||ft2==FASTQ)
	cout<<"======>reading sequence data file: "<<ReadFastq(R2File_name, vec_seq, quality)<<endl;

  //cout<<"1/"<<vec_seq.size()<<":"<<vec_seq.at(0).toString()<<endl;
  //cout<<"1000/1000;"<<vec_seq.at(999).toString()<<endl;

  SequenceString ss1(vec_seq.at(0));
  cout<<"ss1:"<<ss1.toString()<<endl;
  /*cout<<"Name:"<<ss1.GetName()<<endl;
  cout<<"++++++++++++++++++++++++"<<endl;
  cout<<"**extracted and inserted umi:"<<insertUmi2SName_umitools(ss1.GetName(), "ACATGGG")<<endl;
  */
  unsigned seq_str_start=0, seq_str_end=0;
  cout<<"--align sequence:"<<prepAlignSequence(ss1.GetSequence(), OffsetForward, 
		anchor_positions, anchor_length, seq_str_start, seq_str_end)<<endl;
  cout		<<"\n\tSeq start: "<<seq_str_start<<"; end: "<<seq_str_end
		<<endl;
  /*
  //+++++++++testing alignment
  cout<<"====================="<<endl;
  
  SequenceString p("p","TACACCTCGGGTT");
  SequenceString s("s","ACATNNGGGNT");
  cout<<"the pattern string: "<<p.toString()<<endl;
  cout<<"the subject string: "<<s.toString()<<endl;
  OverlapAlignment ol(&p, &s, sm, gapopen, gapextension
						,scale,0);
	AlignmentString aso=ol.GetAlignment(); //pattern is on top; and subject is on bottom.
	unsigned p_start=aso.GetPatternIndexStart();
	unsigned p_end=aso.GetPatternIndexEnd();
	unsigned s_start=aso.GetSubjectIndexStart();
	unsigned s_end=aso.GetSubjectIndexEnd();
	cout<<"p_start:"<<p_start<<"--p_end:"<<p_end<<endl;
	cout<<"pattern:"<<aso.GetPattern()<<endl;
	cout<<"s_start:"<<s_start<<"--s_end:"<<s_end<<endl;
	cout<<"pattern:"<<aso.GetPattern()<<endl;
	cout<<"subject:"<<aso.GetSubject()<<endl;
	cout<<"alignment:\n"<<aso.toString()<<endl;

	// ******************
	cout<<"tesiting extracting umi........."<<endl;
	cout<<"\t Get GetUmiFromAlignmentString:"<<GetUmiFromAlignmentString(aso)<<endl;
	cout<<"testing trim sequence ....."<<endl;
	cout<<"\trimming the sequence at "<<11<<endl;
	cout<<"\t "<<trimSequenceUmiAnchorFromStart(p.GetSequence(), 11)<<endl;
*/	
	//*********testing extracting the umi
	cout<<"*****************testing extraction umi "<<endl;
	string ex_umi;
	bool t=extractUmi(ss1, umi_pattern, 
				anchor_positions, anchor_length,
				umi_positions, umi_length,
				trim, extract, 
				sm, gapopen, gapextension,
				scale,
				Mismatches, OffsetForward,
				ex_umi
				);
	if(t)
		cout<<"extraction success!!"<<endl;
	else
		cout<<"extractoin is not successful"<<endl;
	cout<<"the extracted umi:"<<ex_umi<<endl;
/*  cout<<"===>reading and processing adaptor sequences for mapping......"<<endl;
  if(trim==1)
    {
      cout<<"set up the trim flag"<<endl;
      SetUpTrimFlag(true);
    }
//  SetUpByIsotypeOutputFlag(isotype_flag);
//  ProcessAdaptorSequences(adaptorFile_name, barcodeFile_name, forwardFile_name, reverseFile_name, vec_forward_seq, vec_reverse_seq);

  
  
  WriteFasta("forward_adaptor_sequnces_set.fa",vec_forward_seq);
  WriteFasta("reverse_adaptor_sequnces_set.fa",vec_reverse_seq);  

  //SetupConstantPrimerInfo(1,forwardFile_name);
//now we have everything, we just need to do the job, I mean mapping, here.
  cout<<"Start doing the mapping..........."<<endl;
//  MappingAdaptors(vec_forward_seq, vec_reverse_seq, vec_seq, 
//		  sm, gapopen, gapextension, trim,
//		  MismatchRateThreshold, MinimumOverlapLength, OffsetForward, OffsetReverse, 
//		  outputFileBoth_name, outputFileForward_name, outputFileReverse_name,outputFileNone_name
//		  ); 
  
*/
  //-------------clean up the memory
  delete [] umi_positions;
  delete [] anchor_positions;
  cout<<"Done!!!"<<endl;
  //ofs.close();

  //cout<<"Total "<<gene_info.size()<<" genes are processed"<<endl;
  cout<<"Thanks for using our program and have a nice day!!"<<endl;
  return 0;
}


static void parseArguments(int argc, char **argv, const char *opts)
{
	//"hvp:s:r:o:q:b:d:g:e:m:txk:n:f:"
  int optchar;

  while ((optchar = getopt(argc, argv, opts)) != -1)
    {
      switch(char(optchar))
	{
	case 's':
	  R1File_name=optarg;
	  break;
	case 'r':
	  R2File_name=optarg;
	  break;
	case 'p':
	  umi_pattern=optarg;
	  break;
	case 'o':
	  OR1File_name=optarg;
	  break;
    case 'q':
	  OR2File_name=optarg;
	  break;
	  
	case 'b':
	  OR1FileNU_name=optarg;
	  break;
	case 'd':
	  OR2FileNU_name=optarg;
	  break;  

	case 'm':
		cout<<"NOTE: this has been hard coded to be \"nuc44DM1\"!!!!"<<endl;
		cout<<"\t\t ********can not to be changed"<<endl;
	  //scoreMatrixName=optarg;
	  scoreMatrixName="nuc44DM1";
	  break;
	  
	case 't':
	  trim=true;
	  break;
	case 'x':
	  extract=true;
	  break;
	case 'k':
	  scale=atof(optarg);
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
	  
	case 'n':
	  Mismatches=atof(optarg);
	  break;
	case 'f':
	  OffsetForward=atoi(optarg);
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
 //	"hvp:s:r:o:q:b:d:g:e:m:txk:n:f:"
	
  //o: ouput data file R1 with umi found and processed 
  //q: ouput data file R2 with umi found and processed
  //
  //b: output data file R1 with seqeunce no umi found  
  //d: output data file R2 with sequence no umi found 
  
  //g: gap open
  //e: gap extension
  //m: score matrix
   
  //t: trim the sequence keep the barcode in the sequence. by default we trim it off and set the umi to the sequence name.
  //x: extract the umi and added to the sequence name following the umi_tools style
  //k: scale factor

  //n: # of mismatch allowd 
  //f: biggest offset/shift allowed for the umi (reasonable only when there is an anchors
  cout<<"This is a program used to mimicking the umi_tools extract to extract the umi and put it \n"
	  <<"\tto the sequence header. The reason we implement this piece of code is that in our project\n"
	  <<"\twe have mixed umi vs no-umi sequences. We need a software can distinguish them and separate\n"
	  <<"\tthem. We need a program that can search for an \"anchor\", eg ACATGGG, and then based on the \n"
	  <<"\tanchor to separate sequences and extract the code. Another issue for umi_tools is that \n"
	  <<"\tit has no ability to detect shift of the umi code. This is reasonable if there are only random\n"
	  <<"\tumi bits. But here in our project, we have an anchor, so we can detect shift and find the \n"
	  <<"\tcorrect umi barcodes.\n"
	  <<"The algorithm is like this: we first read/parse the umi string pattern; if the umi contains \n"
	  <<"\tno anchor, things are easy. We can simply turn off offset =0, output file name (no UMI) and \n"
	  <<"\twe only need to take the umi codes at the specific locations/positions; If otherwise, we "
	  <<"\twe need to set up the offset, output file (no UMI) and take the umi codes positions. then the\n"
	  <<"\tnext step to do an alignment and get the best alignment; we take the aligmnet strings and check\n"
	  <<"\tgainst minimum length, # of mismatch. After this we go ahead to extract umi and then we are done\n"
	  <<"This program works for pair end reads, but for now we only extract the umi on R1 (will need to\n"
	  <<"\texpend to add umi on R2). This program can also be used to cut off adaptors without extracting\n"
	  <<"\tumi. "
	  ;
  cout<<"Usage:\n";
  cout<<"\t"<<argv[0]<<" -s sequence file(R1) [-r sequence file (R2)] -p UMI pattern string  \n"
      <<"\t [-o output Sequence file (R1)] [-q output Sequence file (R2)] \n"
	  <<"\t [-b output Sequence file (R1, no UMI found)] [-d output Sequence file (R2, no UMI found)]\n"
      <<"\t  [-m score matrix] [-t (trimming sequences)] [-x (extract the umi and add to the sequence name)] \n"
      <<"\t [-k scale] [-g gapopen panelty] [-e gap extension]\n"
      <<"\t [-n # of Mismatch allowed] [-f bigged offset/shift allowed (where there is an anchor)]\n"     
      <<"\t [-h] [-v]\n\n"
    ;

  cout<<"\t\t-s filename -- the sequence fasta/fastq(.gz) (R1)data filename \n"
      <<"\n";
  cout<<"\t\t-r filename -- the sequence fasta/fastq(.gz) (R2)data filename \n"
      <<"\n";
  cout<<"\t\t-p umi pattern -- the umi pattern string. A character \'N\' is a random barcode \n"
	  <<"\t\t\t; other IUPAC codes indicating anchor bits. An example is \"NNNNNNNNNACATGGG\", which \n"
	  <<"\t\t\t is the current 9 random UMI on Bio-nextera TSO\n"
      <<"\n";

  cout<<"\t\t-o filename -- output sequence file (R1) when either we have sequence trimmed\n"
	  <<"\t\t\tor umi added to the sequence name or both.NOTE: in this case, the sequences \n"
	  <<"\t\t\tincluded are the ones that were found to contain the umi barcode (with anchor)\n"
	  <<"\t\t\t where there is no anchor in the umi, all sequences will be included in this\n"
	  <<"\t\t\t file"
      <<"\n";

  cout<<"\t\t-q filename -- output sequence file (R2) when either we have sequence trimmed\n"
	  <<"\t\t\tor umi added to the sequence name or both. NOTE: in this case, the sequences \n"
	  <<"\t\t\tincluded are the ones that were found to contain the umi barcode (with anchor)\n"
	  <<"\t\t\t where there is no anchor in the umi, all sequences will be included in this\n"
	  <<"\t\t\t file"
      <<"\n";
  cout<<"\t\t-b filename -- output sequence file (R1) contains only sequences that contains\n"
	  <<"\t\t\tno umi-anchor. This could be empty if the sequences are all containing the umi\n"
      <<"\n";

  cout<<"\t\t-d filename -- output sequence file (R2), simliar to -b filename for the file  \n"
	  <<"\t\t\t with no umi"
      <<"\n";
  
  cout<<"\t\t-m scorematrix -- the socre matrix name used for the alignment, \n"
      <<"\t\t\tonly support nuc44DM1. nuc44DM1 by default for nucleotide\n"
      <<"\n";

  cout<<"\t\t-t  -- trim the sequence and save data to output file; \n"
      <<"\n\n";
  cout<<"\t\t-x  -- extract umi and added to the seqeunce name following the umi_tools style; \n"
      <<"\n";
	  
  cout<<"\t\t-k scale -- the scale factor used to set the returned score to the correct unit,\n"
      <<"\t\t\t 1 by default. The programe first uses the scale factor coming with matrix\n"
      <<"\t\t\t  to the return score and the the scale set by this option\n\n";

  cout<<"\t\t-g gapopen -- the gap open value.will be turned into negative if not\n"
      <<"\n";
  cout<<"\t\t-e gapextension -- the gap extension value. will be turned into negative if not\n"
      <<"\n";
  
  cout<<"\t\t-n # -- number of mismatch allowed\n"
      <<"\n"; 
  cout<<"\t\t-p # -- offset/shift allowed.\n"
      <<"\n"; 
  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to cut off the adaptor and extract umi \n"
    <<"\t\t\t @1/11/2020\n"
     ;
  cout<<"\t\t\t*********by Feng @ BU 2020\n"; 


  exit(-1);
}
//return the index that could be used to pick up the ScoreMatrix object from ScoreMatrixArr
//-1 for can not find
static int lookUpScoreMatrix(const string* _scoreMatrixNameArray,const int& len, const string& scoreMatrixName)
{
	//cout<<"|||||testing score matrix index look up......"<<endl;
	//cout<<"\tscoreMatrixName:"<<scoreMatrixName<<endl;
  //int index=-1;
  //int size=sizeof(_scoreMatrixNameArray)/sizeof(_scoreMatrixNameArray[0]);
  for(int i=0;i<len;i++)
    {
		//cout<<"\ti:"<<i<<endl;
      if(scoreMatrixName.compare(_scoreMatrixNameArray[i])==0)
	{
		//cout<<"\tfound a match done."<<endl;
	  return i;
	}
    }
	//cout<<"\tno match ,done"<<endl;
  return -1;
}

