#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "../string_ext.hpp"
//#include "read_gene_info_func.hpp"

//from on you need to specify the c libraries explicitly, 11/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "../SequenceHandler.hpp"
#include "../FastaHandler.hpp"
#include "../SequenceHandlerCommon.hpp"
using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);

static string inputFile_name;//input filename

static string outputFile_name; //the output filename 

int main(int argc, char* argv[])
{
   const char *opts = "hvf:o:";
  //f: the input file name containing the fasta format gene sequence
  //o: output file name

  parseArguments(argc, argv, opts);
  
  if(inputFile_name.size()==0)
    {
      cout<<"please specify the input fasta file name.......\n";
      cout<<"type \"./"<<argv[0]<<"-h\" for usage help\n";
      exit(-1);
      //printUsage(argc, argv);
    }
  if(outputFile_name.size()==0)
    {
	outputFile_name=inputFile_name;
	outputFile_name.append(".out");
      cout<<"no output file name specified, use \""<<inputFile_name<<".out\".......\n";
      //printUsage(argc, argv);
      
    }

	//testing sequence string operator '<'
	SequenceString seq1("a","a");
	SequenceString seq2("b", "b");
	SequenceString seq3("c", "c");
	vector<SequenceString> test_vec;
	test_vec.push_back(seq2);
	test_vec.push_back(seq1);
	test_vec.push_back(seq1);
	test_vec.push_back(seq3);
	sort(test_vec.begin(), test_vec.end());
for(unsigned int i=0;i<test_vec.size();i++)
{
	cout<<test_vec.at(i).GetSequence()<<endl;
}

cout<<"seq1:"<<seq1.GetSequence()<<";Seq2:"<<seq2.GetSequence()<<";Compare <:"<<(seq1<seq2)<<endl;
//now we call read in the fasta file
	vector<SequenceString> vec_seq;
  cout<<"reading in "<<ReadFasta(inputFile_name, vec_seq)<<" sequences from the input file"<<endl;
cout<<"1st sequence:"<<vec_seq.at(0).toString()<<endl;

//for(unsigned int i=0;i<vec_seq.size();i++)
//{
//	cout<<vec_seq.at(i).toString();
//}
//reverse the sequence
for(unsigned int i=0;i<vec_seq.size();i++)
{
	vec_seq[i]=Reverse(vec_seq.at(i));
}

//cout<<"*************after reverse"<<endl;
/*for(unsigned int i=0;i<vec_seq.size();i++)
{
	cout<<vec_seq.at(i).toString();
}*/
cout<<"done..."<<endl;
cout<<"now start sorting........"<<endl;
	sort(vec_seq.begin(), vec_seq.end());
cout<<"done with sorting.........."<<endl;

//reverse the sequence
/*for(unsigned int i=0;i<vec_seq.size();i++)
{
	cout<<vec_seq.at(i).toString();
}*/

//open up the output file to write output
ofstream ofs((outputFile_name).c_str(), ios::out);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<outputFile_name<<"\" can not be opened, quit....\n";
      exit(-1);
    }

cout<<"run through to summarize the abundance........"<<endl;


SequenceString curr_seq("",vec_seq.at(0).GetSequence());
unsigned curr_count=0;

for(unsigned int i=0;i<vec_seq.size();i++)
{
	//cout<<"round "<<i<<endl;
	SequenceString runningSeq=vec_seq.at(i);
	if(curr_seq.GetSequence().compare(runningSeq.GetSequence())==0)
	{
		//cout<<"found a match, update the count"<<endl;
		curr_count++;
	}
	else
	{
		//cout<<"not a match, write the sequence and udpate the counter"<<endl;
		ofs<<curr_seq.GetSequence()<<"\t"<<curr_count<<"\n";
		curr_seq.SetSequence(runningSeq.GetSequence());
		curr_count=1;
	}
	//cout<<test_vec.at(i).GetSequence()<<endl;
}
//last, we need to write down whatever in the counter
ofs<<curr_seq.GetSequence()<<"\t"<<curr_count<<"\n";

//testing the reverse of the string
//SequenceString seqR=Reverse(curr_seq);
//cout<<"before reverse:"<<curr_seq.toString()<<endl;
//cout<<"after reverse:"<<seqR.toString()<<endl;

//close the output file 
ofs.close();
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
	case 'f':
	  inputFile_name=optarg;
	  break;
	  
	case 'o':
	  outputFile_name=optarg;
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
  cout<<"Usage:\n";
  cout<<"\t"<<argv[0]<<" -f input file name [-t output file name] \n"
	<<"\t\t[-h][-v]\n"
    ;

  cout<<"\t\t-f filename -- the input sequence fasta filenames \n"
      <<"\n";

  cout<<"\t\t-o filename -- the output filename, if not provided, use \"inputfilename.out\" \n"
      <<"\n";
  
    
  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to identify the unique sequences and also count their multiplicity\n"
    <<"\t\t\t @09/16/2014\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********by Feng @ BU 2014\n"; 


  exit(-1);
}

