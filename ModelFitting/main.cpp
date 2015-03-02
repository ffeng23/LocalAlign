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

using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);

int main(int argc, char* argv[])
{
  //for parsing commandline arguement

  const char *opts = "hvg:e:o:s:t:k:am:";
  //f: the input file name containing the fasta format gene sequence
  //c: flag to indicate to count each gene length.
  //o: output file name

  parseArguments(argc, argv, opts);

  //start testing

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
	  /*	case 's':
	  inputFile1_name=optarg;
	  break;

	case 't':
	  inputFile2_name=optarg;
	  break;  
	  
	case 'o':
	  outputFile_name=optarg;
	break;

	case 'm':
	  scoreMatrixName=optarg;
	  break;
	  
	case 'a':
	  seqtype=AminoAcid;
	  break;
	case 'k':
	  scale=atoi(optarg);
	  break;
	case 'e':
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
  cout<<"\t"<<argv[0]<<" -s first sequence file -t second sequence file \n"
      <<"\t [-o output filename] [-m score matrix] [-a (protein alignment)] \n"
      <<"\t [-k scale] [-g gapopen panelty] [-e gap extension]\n"
      <<"\t [-n number of local alignment returned]"
    ;

  cout<<"\t\t-s filename -- the first sequence fasta filenames \n"
      <<"\n";

  cout<<"\t\t-t filename -- the second sequence fasta filenames\n"
      <<"\n";

  cout<<"\t\t-o filename -- the filenames contains the local alignments \n"
      <<"\n";
  
  cout<<"\t\t-m scorematrix -- the socre matrix name used for the alignment, \n"
      <<"\t\t\tonly support nuc44 and blosum50. nuc44 by default for nucleotide\n"
      <<"\t\t\tand pam50 by default for protein alignment  \n"
      <<"\n";

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
    
  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to do the alignment in the fasta file\n"
    <<"\t\t\t @11/12/2013\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********by Feng @ BU 2013\n"; 


  exit(-1);
}


