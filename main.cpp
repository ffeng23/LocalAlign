#include <iostream>
#include <fstream>
#include <vector>
#include "string_ext.hpp"
//#include "read_gene_info_func.hpp"

//from on you need to specify the c libraries explicitly, 11/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>


using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);

static string inputFile1_name; //the input file for seq1
static string inputFile2_name;//input ifle for seq2  
//both above files are in fasta file and should contains the same sets of genes.

static string outputFile_name("localalignment.txt"); //the output file for the local alignment results with scores

static string scoreMatrixName; //name as an input for specifing the name of the score matrix. for example nuc44

//static int read_gene_sequence(const string& temp_seq, vector<string>& promoter_sequence);



int main(int argc, char* argv[])
{
  //for parsing commandline arguement
  const char *opts = "hvf:co:";
  //f: the input file name containing the fasta format gene sequence
  //c: flag to indicate to count each gene length.
  //o: output file name

  parseArguments(argc, argv, opts);
  
  /*if(inputFile_name.size()==0)
    {
      cout<<"please specify the input fasta file name.......\n";
      printUsage(argc, argv);
    }
  if(!flag_countLength)
    {
      cout<<">>>>>>WARNING:no action has been specified, quit...........\n";
      exit(-1);
    }

  cout<<"The output file name is \""<<outputFile_name<<".\n";
  ifstream ifs1_p(inputFile1_name.c_str());
  
  ofstream ofs((outputFile2_name).c_str());
  
  
  if(!ifs_p.is_open())
    {
      cout<<">>>>>>ERROR:the input file \""<<inputFile_name<<"\" can not be opened, quit....\n";
      exit(-1);
    }

  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<outputFile_name<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  */
  /*
  //so far, we have open files, now we need to read the file first.
  string line; //out_line;

  vector<string> gene_info;
  vector<string> gene_sequence;
  
  int gene_number=0;
  
  int line_count=0;
  string temp_seq;

  cout<<"read the fasta file........."<<endl;
  cout<<"line#: ";
  //  double d_storage, d_current;
  string temp_string;
  while(!ifs_p.eof())
    {
      line_count++;
      if(line_count/5000*5000==line_count)
	{
	  cout<<"..... "<<line_count;
	  cout.flush();
	}
      //found_one=false;
      getline(ifs_p, line);
      if(line.compare("\n")==0||line.length()==0||line.compare("\t")==0||line.compare("\t\t")==0)
	{
	  cout<<"...read an empty row (skip!)"<<endl;
	  continue;
	}
            //cout<<"here......"<<endl;
      if(line[0]=='>')
	{
	  //this is a new refSeq header line,
	  temp_string=line;
	  gene_info.push_back(temp_string);
	  
	  
	  if(gene_number>0) //we already read something, so we have to store it to the vector
	    {
	      if(read_gene_sequence(temp_seq,gene_sequence)==-1)
		{
		  cout<<">>>>>>>>>>>>ERROR:in reading gene sequences, quit..........."<<endl;
		  exit(-1);
		}
	      
	    }//otherwise, we are doing the first one, so do bother

	  temp_seq.clear();
	  gene_number++;
	  
	}
      else
	{//append the seqence to it
	  temp_seq.append(line);
	  //temp_seq.append("\n");
	}

    }
  
  //here we have to update/add the last gene
  if(read_gene_sequence(temp_seq, gene_sequence)==-1)
    exit(-1);
  

  //we are done with reading the promoter sequence file.
  cout<<"\nfinish reading the promoter file........\n"
      <<"\tsummary: total "<<line_count<<" line read in and \n"
      <<"\t\t"<<gene_number<<" sequences store in the vector...."<<endl;
  */
  //ifs_p.close();
  cout<<"writing output........."<<endl;

  //write_output(ofs, gene_info, gene_sequence);

  cout<<"Done!!!"<<endl;
  //ofs.close();

  //cout<<"Total "<<gene_info.size()<<" genes are processed"<<endl;
  cout<<"Thanks for using our program and have a nice day!!"<<endl;
  return 0;
}

/*void write_output(ostream& ofs, vector<string>& gene_info, vector<string>& gene_sequence)
{
  
  for(unsigned int i=0;i<gene_info.size();i++)
    {
      ofs<<">"<<gene_info[i]<<"\t";
      ofs<<gene_sequence[i].size()<<endl;
    }
    }*/

//here in this function, we append the sequece into the gene sequece vector
//-1 for error,
//0 for good.
/*static int read_gene_sequence(const string& temp_seq, vector<string>& gene_sequence)
{
  gene_sequence.push_back(temp_seq);
  return 0;
  }*/

static void parseArguments(int argc, char **argv, const char *opts)
{
  int optchar;

  while ((optchar = getopt(argc, argv, opts)) != -1)
    {
      switch(char(optchar))
	{
	  
	case 's':
	  inputFile1_name=optarg;
	  break;

	case 't':
	  inputFile2_name=optarg;
	  break;  
	  
	case 'o':
	  outputFile_name=optarg;
	break;


	case 'm':
	 ScoreMatrixName=optarg;
	break;
	  
	
	  /*case '':
	  refSeq_size= atoi(optarg);
	  break;
	  	case 'e':
	  to_id = atoi(optarg);
	  break;
	
	  case 'p':
	  delimitor=optarg[1];
	  break;
	    
	case 'x':
	  if(parseNumberString(optarg, unknown_fields)!=0)
	    {
	      printUsage(argc,argv);
	      exit(-1);
	    }
	  //commond_field_flag=1;	  
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
  cout<<"\t"<<argv[0]<<"-s first sequence file -t second sequence file \n"
      <<"[-o output filename] [-m score matrix] [-a (protein alignment)] \n"
      <<"[-k scale] [-g gapopen panelty] [-e gap extension] [-n number of local alignment returned]"
    ;

  cout<<"\t\t-s filename -- the first sequence fasta filenames \n"
      <<"\n";

  cout<<"\t\t-t filename -- the second sequence fasta filenames\n"
      <<"\n";

  cout<<"\t\t-o filename -- the filenames contains the local alignments \n"
      <<"\n";
  
  cout<<"\t\t-m scorematrix -- the socre matrix name used for the alignment, only support nuc44 and \n"
      <<"\t\t\t nuc44 by default for nucleotide alignment and pam50 by default for protein alignment  \n"
      <<"\n";

  cout<<"\t\t-a  -- the type of alignment a for protein alignment, aa; without setting this one,\n"
      <<"\t\t\t it is by default doing nucleotide alignment. This type decides which default matrix will be used\n"
      <<"\n";

  cout<<"\t\t-k scale -- the scale factor used to set the returned score to be of correct unit, by default 1\n"
      <<"\t\t\t the programe first used the scale factor coming with matrix to the return score and the the scale set by this option\n"
      <<"\n";

  cout<<"\t\t-g gapopen -- the second sequence fasta filenames\n"
      <<"\n";
    
  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to count the gene lenghth in the fasta file @11/12/2013\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********by Feng @ BU 2013\n"; 


  exit(-1);
}

