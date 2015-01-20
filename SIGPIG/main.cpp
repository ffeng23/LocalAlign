#include <iostream>
#include <fstream>
#include <vector>
#include "../string_ext.hpp"

//#include "read_gene_info_func.hpp"

//from on you need to specify the c libraries explicitly, 11/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#include "../score.hpp"
#include "../SequenceString.hpp"
#include "../LocalAlignment.hpp"
#include "../GlobalAlignment.hpp"
#include "../OverlapAlignment.hpp"
#include "../FastaHandler.hpp"
#include "../SequenceHandler.hpp"
#include "../AlignmentString.hpp"

#include "GenomicJ.hpp"
#include "GenomicD.hpp"
#include "GenomicV.hpp"
#include "genomicSegments.hpp"
#include "LoadData.hpp"
#include "AlignmentSettings.hpp"
#include "Alignment.hpp"
//#include 

using namespace std;

static string scoreMatrixName="nuc44"; //name as an input for specifing the name of the score matrix. for example nuc44

static string supportedScoreMatrixNameArr[]={"nuc44","blosum50", "tsm1", "tsm2", "nuc44HP"};
//tsm2: score matrix from Geoffrey Barton reference.

static ScoreMatrix* ScoreMatrixArr[]={&nuc44, &blosum50, &tsm1, &tsm2, &nuc44HP};
double gapopen=-8;
double gapextension=-8;
//static double scale=1; //this is the one on top of matrix, the programe will run score and use the 
//matrix specified scale first and then apply the scale set by this one.

double errorCost=-4;
string seqFileName;
string vSegmentFileName("genomicVs_alleles.fasta");;
string dSegmentFileName("genomicDs.fasta");
string jSegmentFileName("genomicJs_all_curated.fasta");

string outputFileNameBase;
unsigned numberOfSeqProcessedEach;

int numberOfThread=2;//the number of working thread, doesn't included 
pthread_mutex_t progressMutex;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);

int main(int argc, char* argv[])
{
  //for parsing commandline arguement

  const char *opts = "hiv:d:j:c:s:o:n:t:";
  //f: the input file name containing the fasta format gene sequence
  //c: flag to indicate to count each gene length.
  //o: output file name

  numberOfSeqProcessedEach=AlignmentSettings::N_per_file;

  parseArguments(argc, argv, opts);
  if(numberOfThread<1)
    {
      cout<<"The number of working threads requested is not valid.......\n";
      cout<<"type \""<<argv[0]<<"\" -h\" for usage help\n";
      exit(-1);
      
    }
  
  if(seqFileName.size()==0)
    {
      cout<<"please specify the seq input fasta file name.......\n";
      cout<<"type \"./testalign -h\" for usage help\n";
      exit(-1);
      //printUsage(argc, argv);
    }

  if(outputFileNameBase.size()==0)
    {
      outputFileNameBase=seqFileName;
    }
  cout<<"***Input parameters Summary:\n";
  cout<<"\tInput file name:\""<<seqFileName<<"\".\n";
  cout<<"\tOutput file name (base):\""<<outputFileNameBase<<"\".\n";
  cout<<"\tV gene segment file:\""<<vSegmentFileName<<"\"\n";
  cout<<"\tD gene segment file:\""<<dSegmentFileName<<"\"\n";
  cout<<"\tJ gene segment file:\""<<jSegmentFileName<<"\"\n";
  cout<<"\tThe error cost: "<<errorCost<<"\n";
  cout<<"\tThe number of sequences processed for each file: "<<numberOfSeqProcessedEach<<"\n";
  
  AlignmentSettings::N_per_file=numberOfSeqProcessedEach;

  //testing the alignment string
  ScoreMatrix* sm= ScoreMatrixArr[0];
  char c1='A', c2='T';

  cout<<"\t("<<c1<<","<<c2<<")="<<sm->GetScore(c1,c2)<<endl;
  
  //testing local alignment
  //SequenceString Seq2("seq1","AGCTAGAGACCCCAGTCTGAGGTAGA");
  //SequenceString Seq1("seq2", "AGCTAGAGACCAGCTATCTAGAGGTAGA");
  //SequenceString Seq1("seq1","ACCCCAG");
  //SequenceString Seq2("seq2", "ACCAG");


  SequenceString Seq1("seq1","ATGAGGTAGA");
  SequenceString Seq2 ("seq2", "CTAGAGGTAGA");

  cout<<"showing sequence string\n"<<Seq1.toString()<<Seq2.toString()<<endl;

  //SequenceString tempSStr=ReverseComplement(Seq2);
  
  //now testing alignment
  cout<<"Testing alignment:"<<endl;

  cout<<"\t("<<c1<<","<<c2<<")="<<sm->GetScore('A','T')<<endl;
  
  //LocalAlignment la(&Seq1,&Seq2,sm, gapopen, gapextension,1, 5);
  //LocalAlignment la(&Seq1,&tempSStr,sm, gapopen, gapextension,1, 100);
  /*cout<<"\tdone and the score is "<<la.GetScore()<<endl;
  cout<<"\t"<<la.GetAlignment().toString()<<endl;

  //testing multiple alignment
  cout<<"Total number of local aligments:"<<la.GetNumberOfAlignments()<<endl;
  for(unsigned int i=0;i<la.GetNumberOfAlignments();i++)
    {
      cout<<i<<"/"<<la.GetNumberOfAlignments()<<":"<<endl;
      cout<<la.GetAlignmentArr()[i].toString()<<endl;
    }
  */
  
  //testing globalAlignment
  cout<<"Testing global alignment:"<<endl;
  //GlobalAlignment gla(&Seq1,&Seq2,sm, gapopen, gapextension,1);
  OverlapAlignment ola(&Seq1,&Seq2,sm, gapopen, gapextension,1);
  cout<<"\tdone and the score is "<<ola.GetScore()<<endl;
  cout<<"\t********"<<ola.GetAlignment().toString()<<endl;
  
  
  //testing overlapAlignment
  /*
  cout<<"Testing overlap alignment:"<<endl;
  OverlapAlignment ola(&Seq1,&Seq2,sm, gapopen, gapextension,1);
  
  //OverlapAlignment ola(&Seq1,&tempSStr,sm, gapopen, gapextension,1);
  cout<<"\tdone and the score is "<<ola.GetScore()<<endl;
  cout<<"\t"<<ola.GetAlignment().toString()<<endl;
  */
  cout<<"done"<<endl;
  

  //testing fasta handler
  /*vector<SequenceString> vec_seq;
  cout<<"reading in "<<ReadFasta(inputFile1_name, vec_seq)<<endl;

  cout<<"1/1000:"<<vec_seq.at(0).toString()<<endl;
  cout<<"999/1000;"<<vec_seq.at(999).toString()<<endl;

  WriteFasta("feng.fa",vec_seq);

  //testing reverseComp
  cout<<"Testing reverse comp:"<<endl;
  SequenceString tempSString=ReverseComplement(vec_seq.at(0));

  cout<<"rev comp 0/1000"<<tempSString.toString()<<endl;

  cout<<"Test substr()"<<endl;
  string tempString("feng");
  cout<<"0,1 : "<<tempString.substr(0,1)<<endl;
  cout<<"0,-100: "<<tempString.substr(0,100)<<endl;
  cout<<"length is "<<tempString.length()<<endl;
  cout<<"4,1 :" <<tempString.substr(4)<<endl;
  */
  cout<<"writing output........."<<endl;

  //start doing the reading of genomic segments
  //first we need to declare the genomic segments
  GenomicV* genV=NULL;
  GenomicJ* genJ=NULL;
  GenomicD* genD=NULL;
  
  unsigned totalNumJ=  ReadGenomicJ("genomicJs_all_curated.fasta",&genJ);
  
  cout<<totalNumJ<<" J genomic segments are read in."<<endl; 
  cout<<"\tshowing the J seq 1:"<<genJ[1].Get_Sequence()<<endl;
  
  unsigned totalNumD=  ReadGenomicD("genomicDs.fasta",&genD);
  
  cout<<totalNumD<<" D genomic segments are read in."<<endl; 

  unsigned totalNumV=  ReadGenomicV("genomicVs_alleles.fasta",&genV);
  
  cout<<totalNumV<<" V genomic segments are read in."<<endl; 


  //now testing load sequence data
  vector<SequenceString> data_vec;
  vector<string> header_vec;
  vector<unsigned> count_vec;
  unsigned totalNumSequences=LoadData("sample.data", header_vec, data_vec, count_vec);
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

  //testing parse field function
  cout<<"Testing parsefield function:\n";
  unsigned N_read=data_vec.size();
  vector<SequenceString>all_Sequences=data_vec;
  unsigned Read_Length=ParseField(header_vec, "Read_Length");
  if((signed)Read_Length!=-1)
    {
      cout<<"\tread length:"<<Read_Length<<endl;
    }
  else
    {
      cout<<"\tread length:"<<-1<<endl;
    }
  cout<<"\tnumber of read:"<<N_read<<endl;

  //determine the output file names.
  vector<string> outFileNames;
  outFileNames=DetermineOutputFileNames(outputFileNameBase, AlignmentSettings::N_per_file, N_read);

  //start testing the alignment, first mathJ
  cout<<"%%%%%%%%%%%testing matchJ()>>>>>>"<<endl;
  SequenceString test_seq=all_Sequences.at(0);
  cout<<"\ttesting sequence 0:"<<endl;
  cout<<test_seq.toString()<<endl;
  
  double error_cost=5;
  Alignment_Object J_obj(totalNumJ);
  //now calling it
  bool seq_j_ok=match_J(test_seq, genJ, totalNumJ, AlignmentSettings::J_minimum_alignment_length, AlignmentSettings::J_maximum_deletion, AlignmentSettings::negative_excess_deletions_max, AlignmentSettings::J_allowed_errors, error_cost, J_obj);
  //now we check the output
  if(J_obj.numOfAligned>0)
    {
      cout<<"\tthe number of aligned:"<<J_obj.numOfAligned<<endl;
      cout<<"\tthe aligned length:"<<J_obj.align_length[0]<<","
	  <<J_obj.align_length[1]<<","<<J_obj.align_length[2]<<endl;
      cout<<"\tthe min_deletions:"<<J_obj.min_deletions[0]<<"."
	  <<J_obj.min_deletions[1]<<","<<J_obj.min_deletions[2]<<"."<<endl;
      cout<<"\tthe n_errors:"<<J_obj.n_errors[0]<<","
	  <<J_obj.n_errors[1]<<","<<J_obj.n_errors[2]<<","<<endl;
    }
  else
    {
      cout<<"failed alignment"<<endl;
    }
  /*
  //now we are ready to do the alignment??
  //first need to figure out number of output files
  unsigned numOfOutFiles=outFileNames.size();

  //declare/initialize the working threads
  pthread_t* workThreads;
  int ret;
  ret = pthread_mutex_init(&progressMutex, NULL);
      if(ret != 0) {
	perror("mutex init failed\n");
	exit(EXIT_FAILURE);
      }


  for(unsigned i=0;i<numOfOutFiles;i++) //outer for loops for the alignment
    {
      //determine the number of sequences
      unsigned numOfThisFile=(N_read-i*AlignmentSettings::N_per_file);
      
      //allocate the threads
      workThreads=new pthread_t[numberOfThread];
      
      //prepare the output
      for(int j=0;j<numberOfThread;j++)
	{
	  ret =  pthread_create(workThreads[j], NULL, thread_func, NULL);
	  if(ret != 0) {
	    perror("pthread_create failed\n");
	    exit(EXIT_FAILURE);
	  }
	}
      
      
      if(numOfThisFile>AlignmentSettings::N_per_file)
	{
	  numOfThisFile=AlignmentSettings::N_per_file;
	}
      //now determine the sequences to be used in this files
      unsigned startIndexOfSequence=i*AlignmentSettings::N_per_file;
      startIndexOfSequence=startIndexOfSequence+0;
      //now here is where we need to figure out the thread thing

      // *******cleaning up the threads
      delete[] workThreads;
    }
  */
  //********clean up
  cout<<"Clean up the memories....."<<endl;
  if(genV!=NULL)
    delete [] genV;
  if(genD!=NULL)
    delete [] genD;
  if(genJ!=NULL)
    delete [] genJ;

  cout<<"Done!!!"<<endl;

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

	case 'v':
	  vSegmentFileName=optarg;
	  break;  

	case 'd':
	  dSegmentFileName=optarg;
	  break;

	case 'j':
	  jSegmentFileName=optarg;
	  break;
	  
	case 'o':
	  outputFileNameBase=optarg;
	break;

	case 'n':
	  numberOfSeqProcessedEach=atoi(optarg);
	  break;
	case 'c':
	  errorCost=atoi(optarg);
	  errorCost*=-1;
	  break;
	case 't':
	  numberOfThread=atoi(optarg);
	  break;

	  /*case 'e':
	  gapextension=atoi(optarg);
	  if (gapextension>0)
	    gapextension*=-1;
	  gapextensionFlag=true;
	  break;
	case 'g':
	  gapopen=atoi(optarg);
	  if(gapopen>0)
	    gapopen*=-1;
	  if(!gapextensionFlag)
	    gapextension=gapopen;
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
	
	case 'i':  
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
  cout<<"\t"<<argv[0]<<"-s sequence inputfile [-v gene segment file] \n" 
      <<"\t [-d D gene segment inputfile] [-j J gene segment file] \n"
      <<"\t [-o base output filename] [-n number of sequences processed for each file] [-c error cost] \n"
    //<<"\t [-k scale] [-g gapopen panelty] [-e gap extension]\n"
    // <<"\t [-n number of local alignment returned]"
    ;

  cout<<"\t\t-s filename -- the sequence file as the input \n"
      <<"\n";

  cout<<"\t\t-v/d/j filename -- the V/D/J gene segment files, by default, \n"
      <<"\t\t \"genomicVs_alleles.fasta\", \"geneomicJs_all_curated.fasta\" and "
      <<"\t\t \"genomicDs.fasta\""
      <<"\n";

  cout<<"\t\t-o filename -- the files contains output alignments. this is \n"
      <<"\t\t the base name. By default, it is the name of the input file.\n"
      <<"\t\t will also be suffixed by the numbers based on the total number of sequences\n"
      <<"\t\t as well as the number of processed for each file (-n option)\n"
      <<"\n";
  
  cout<<"\t\t-n number -- the number of sequences processed by each output file\n"
      <<"\t\t\t this number is also specified in the aligment setting file.\n"
      <<"\t\t\t so only if we want to overwrite this number, specify it here\n"
      <<"\n";

  cout<<"\t\t-t number -- the number of working threads for doing the alignment\n"
      <<"\t\t\t Be careful not to setting a number too big. Not bigger than\n"
      <<"\t\t\t the actual physical CPUs/processors. Please run script to \n"
      <<"\t\t\t determine the number of the physical CPUs and logical CPUs\n"
      <<"\t\t\t (they are not the same when the hyper-threading is enabled\n"
      <<"\t\t\t Also, you need to specify at least 1, \n"
      <<"\t\t\t (total number of thread = main+one working thread)\n"
      <<"\t\t\t main thread is managing, reporting, writing output, etc.\n"
      <<"\t\t\t only the working thread(s) doing the alignments."
      <<"\n";

  /*
  cout<<"\t\t-a  -- the type of alignment a for protein alignment, aa; \n"
      <<"\t\t\t without this one, it is by default doing nucleotide alignment.\n"
      <<"\t\t\t  This type decides which default matrix will be used\n\n";

  cout<<"\t\t-k scale -- the scale factor used to set the returned score to the correct unit,\n"
      <<"\t\t\t 1 by default. The programe first uses the scale factor coming with matrix\n"
      <<"\t\t\t  to the return score and the the scale set by this option\n\n";
  */

  cout<<"\t\t-c error cost -- the cost for a miss match. will be turned into negative if not\n"
      <<"\n";
  /*  cout<<"\t\t-e gapextension -- the cost to open a gap.\n"
      <<"\t\t\t -8 by default. The gap value has to negative.\n"
      <<"\t\t\t if a non-negative value is specified, it will be\n"
      <<"\t\t\t turned into negative (by being multiplied by -1)\n"
      <<"\n";
  */
  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to do the alignment\n"
    <<"\t\t\t @10/20/2014\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********by Feng @ BU 2014\n"; 


  exit(-1);
}

