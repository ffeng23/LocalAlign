#include <iostream>
#include <fstream>
#include <vector>
#include "../string_ext.hpp"

//from on you need to specify the c libraries explicitly, 11/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#include "../SequenceString.hpp"
#include "../FastaHandler.hpp"

using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);

int main(int argc, char* argv[])
{
  //for parsing commandline arguement

  const char *opts = "hiv:d:j:c:s:o:n:t:r";
  //f: the input file name containing the fasta format gene sequence
  //c: flag to indicate to count each gene length.
  //o: output file name
  parseArguments(argc, argv, opts);
  

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
	  //seqFileName=optarg;
	  break;

	case 'n':
	  //numberOfSeqProcessedEach=atoi(optarg);
	  break;
	case 'c':
	  /*errorCost=atoi(optarg);
	  if(errorCost<0)
	    errorCost*=-1;
	  */
	  break;

	case 'r':
	  //flip_flag=true;
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


  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to do the alignment\n"
    <<"\t\t\t @10/20/2014\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********by Feng @ BU 2014\n"; 


  exit(-1);
}

