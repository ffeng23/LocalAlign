#include <iostream>
//#include <fstream>
//#include <vector>
//#include "string_ext.hpp"
//#include "read_gene_info_func.hpp"

//from on you need to specify the c libraries explicitly, 11/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "VDJ_cuts_insertion_dinuc_ntbias_model.hpp"
#include "Assigns.hpp"
#include "Counter.hpp"

#include <fstream>
#include "../string_ext.hpp"

#include <pthread.h>

//#include "../score.hpp"
#include "../SequenceString.hpp"
//#include "../LocalAlignment.hpp"
//#include "../GlobalAlignment.hpp"
//#include "../OverlapAlignment.hpp"
#include "../FastaHandler.hpp"
//#include "../SequenceHandler.hpp"
//#include "../AlignmentString.hpp"

#include "../SIGPIG/GenomicJ.hpp"
#include "../SIGPIG/GenomicD.hpp"
#include "../SIGPIG/GenomicV.hpp"
#include "../SIGPIG/genomicSegments.hpp"
//#include "LoadData.hpp"VD
//#include "../SIGPIG/AlignmentSettings.hpp"
#include "../SIGPIG/Alignment.hpp"
#include "../SIGPIG/Alignment_V.hpp"
#include "../SIGPIG/Alignment_D.hpp"
//#include "do_VDJ_alignment.hpp"
#include "../SIGPIG/DataIO.hpp"
#include "do_probabilistic_model_fitting.hpp"

using namespace std;

string vSegmentFileName("genomicVs_alleles.fasta");
string dSegmentFileName("genomicDs.fasta");
string jSegmentFileName("genomicJs_all_curated.fasta");

string alignmentFileName("sample.data_0.align");
string outputFileNameBase;
unsigned numOfThread=2;

bool start_from_flat_prior=true;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);

//in this version, I am doing single threading processing
int main(int argc, char* argv[])
{
  
  //for parsing commandline arguement

  const char *opts = "hv:j:d:s:o:t:";
  //f: the input file name containing the fasta format gene sequence
  //c: flag to indicate to count each gene length.
  //o: output file name

  parseArguments(argc, argv, opts);

  //for input parameters
if(alignmentFileName.size()==0)
    {
      cout<<"please specify the seq input fasta file name.......\n";
      cout<<"type \"./"<<argv[0]<<" -h\" for usage help\n";
      exit(-1);
      //printUsage(argc, argv);
    }
// unsigned dim
//Maxtrix<double> x(2,
  if(outputFileNameBase.size()==0)
    {
      outputFileNameBase=alignmentFileName;
    }
  cout<<"***Input parameters Summary:\n";
  cout<<"\tInput alignment file name:\""<<alignmentFileName<<"\".\n";
  cout<<"\tOutput file name (base):\""<<outputFileNameBase<<"\".\n";
  cout<<"\tV gene segment file:\""<<vSegmentFileName<<"\"\n";
  cout<<"\tD gene segment file:\""<<dSegmentFileName<<"\"\n";
  cout<<"\tJ gene segment file:\""<<jSegmentFileName<<"\"\n";
  //cout<<"\tThe error cost: "<<errorCost<<"\n";
  //cout<<"\tThe number of sequences processed for each file: "<<numberOfSeqProcessedEach<<"\n";
  //start testing
  cout<<"printing...........-1"<<(unsigned)-1<<endl;
  //first open up the files============================
  //first we need to declare the genomic segments
  GenomicV* genV=NULL;
  GenomicJ* genJ=NULL;
  GenomicD* genD=NULL;
  unsigned totalNumV, totalNumD, totalNumJ;
  if(!DoGenomicTemplateReading(vSegmentFileName, dSegmentFileName, jSegmentFileName,
			      &genV, &genD, &genJ, totalNumV, totalNumD, totalNumJ))
    {
      cout<<"Error in reading genomic template files, quit"<<endl;
      exit(-1);
    }
  for(unsigned i=0;i<totalNumD;i++)
    {
      cout<<genD[i].toString()<<endl;//Get_n_alleles()<<",";
    }
  cout<<endl;
  //==================================
  //now deserialize the alignment file
//testing the deserialization
  Alignment_Object* v_align=NULL;
  Alignment_D* d_align=NULL;
  Alignment_Object* j_align=NULL;
  unsigned total_alignment=0;
  SequenceString* seq=NULL;
  cout<<"doing deserialization...."<<endl;
  if(!DoDeserialization(alignmentFileName, total_alignment,  &v_align, 
			&d_align, &j_align, &seq)
     )
    {
      cout<<"Error in deserializing the alignment"<<endl;
      exit(-1);
    } 
  cout<<"Succefully deserialize the alignments"<<endl;
  cout<<"Deserialization summary:"<<endl;
  cout<<"\tnumber of aligments read in:"<<total_alignment<<endl;

  //now we start building the model and get counter and assigns
  bool start_from_flat_prior=true;
  unsigned numOfIters=1;
  bool mfit=do_probabilistic_model_fitting
    (seq, v_align, d_align, j_align, total_alignment, genV, totalNumV,
     genD, totalNumD, genJ, totalNumJ, start_from_flat_prior,numOfIters);
  if(mfit)
    {
      cout<<"Successfully done the fitting!"<<endl;
    }
  else
    {
      cout<<"model fitting failed!"<<endl;
    }
     
  
  cout<<"Thanks for using our program and have a nice day!!"<<endl;

  //clean up the memory
  delete[] genV;
  delete[] genD;
  delete[] genJ;
  
  delete[] v_align;
  delete[] d_align;
  delete[] j_align;
  delete[] seq;
  
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
	  alignmentFileName=optarg;
	  break;

	case 't':
	  numOfThread=atoi(optarg);
	  break;  
	  
	case 'o':
	  outputFileNameBase=optarg;
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
	
	  //case 'v':  
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
  cout<<"\t"<<argv[0]<<" [-s first sequence file] [-t thread #] \n"
      <<"\t [-o output filename] [-v vGenSegment file] [-d dGeneSegment file] \n"
      <<"\t [-j jGeneSegment file] \n"
    //<<"\t [-n number of local alignment returned]"
    ;

  cout<<"\t\t-s filename -- the binary alignment filenames \n"
      <<"\n";

  cout<<"\t\t-t # -- the number of thread to use\n"
      <<"\n";

  cout<<"\t\t-o filename -- the output filenames \n"
      <<"\n";
  
  cout<<"\t\t-v/d/j -- v/d/j gene segment data file. has to been certain format\n"
    //<<"\t\t\tonly support nuc44 and blosum50. nuc44 by default for nucleotide\n"
    //<<"\t\t\tand pam50 by default for protein alignment  \n"
      <<"\n";
  /*
  cout<<"\t\t-a  -- the type of alignment a for protein alignment, aa; \n"
      <<"\t\t\t without this one, it is by default doing nucleotide alignment.\n"
      <<"\t\t\t  This type decides which default matrix will be used\n\n";

  cout<<"\t\t-k scale -- the scale factor used to set the returned score to the correct unit,\n"
      <<"\t\t\t 1 by default. The programe first uses the scale factor coming with matrix\n"
      <<"\t\t\t  to the return score and the the scale set by this option\n\n";

  cout<<"\t\t-g gapopen -- the cost to open a gap.\n"
      <<"\t\t\t -8 by default. The gap value has to e negative.\n"
      <<"\t\t\t if a non-negative value is specified, it will be\n"
      <<"\t\t\t turned into negative (by being multiplied by -1)\n"
      <<"\n";
  cout<<"\t\t-e gapextension -- the cost to open a gap.\n"
      <<"\t\t\t -8 by default. The gap value has to negative.\n"
      <<"\t\t\t if a non-negative value is specified, it will be\n"
      <<"\t\t\t turned into negative (by being multiplied by -1)\n"
      <<"\n";
  */
  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to do model fitting using alignment file\n"
    <<"\t\t\t @3/12/2015\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********by Feng @ BU 2015\n"; 


  exit(-1);
}


