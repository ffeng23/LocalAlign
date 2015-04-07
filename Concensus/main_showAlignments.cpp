#include <iostream>
#include <fstream>
#include <vector>
#include "../string_ext.hpp"

//from on you need to specify the c libraries explicitly, 11/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
//#include <pthread.h>

#include "../SequenceString.hpp"
#include "../FastaHandler.hpp"

//#include "../SIGPIG/LoadData.hpp"
#include "../string_ext.hpp"

//#include <pthread.h>

//#include "../score.hpp"
#include "../SequenceString.hpp"
//#include "../LocalAlignment.hpp"
//#include "../GlobalAlignment.hpp"
//#include "../OverlapAlignment.hpp"
//#include "../FastaHandler.hpp"
//#include "../SequenceHandler.hpp"
//#include "../AlignmentString.hpp"

#include "../SIGPIG/GenomicJ.hpp"
#include "../SIGPIG/GenomicD.hpp"
#include "../SIGPIG/GenomicV.hpp"
#include "../SIGPIG/genomicSegments.hpp"
//#include "LoadData.hpp"VD
#include "../SIGPIG/AlignmentSettings.hpp"
#include "../SIGPIG/Alignment.hpp"
#include "../SIGPIG/Alignment_V.hpp"
#include "../SIGPIG/Alignment_D.hpp"
//#include "do_VDJ_alignment.hpp"
#include "../SIGPIG/DataIO.hpp"
//#include "do_probabilistic_model_fitting.hpp"

using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);

string alignmentFileName("");
string outFileName("");

string vSegmentFileName("genomicVs_alleles.fasta");
string dSegmentFileName("genomicDs.fasta");
string jSegmentFileName("genomicJs_all_curated.fasta");
unsigned numOfSeqSave=0;
int main(int argc, char* argv[])
{
  //for parsing commandline arguement

  const char *opts = "hs:o:v:d:j:n:";
  parseArguments(argc, argv, opts);

  if(alignmentFileName.length()==0)
    {
      cout<<"please specify the input file........quit."<<endl;
      exit(-1);
    }
  if(outFileName.length()==0)
    {
      outFileName=alignmentFileName+"_summary.txt";
    }

  //now check for validity
  ifstream ifs_ss(alignmentFileName.c_str(), std::ios::in);
  if(!ifs_ss.is_open())
    {
      cout<<"******ERROR: can not open file, quit..."<<endl;
      exit(-1);
    }
  ifs_ss.close();
  ofstream ofs_ss(outFileName.c_str(), std::ios::out);
  if(!ofs_ss.is_open())
    {
      cout<<"********Error: can not open out file, quit......"<<endl;
      exit(-1);
    }
  //ofs_ss.close();
  //***********************


  //first read the alignment file
  cout<<"***Input parameters Summary:\n";
  cout<<"\tInput alignment file name:\""<<alignmentFileName<<"\".\n";
  cout<<"\tOutput file name (base):\""<<outFileName<<"\".\n";
  cout<<"\tV gene segment file:\""<<vSegmentFileName<<"\"\n";
  cout<<"\tD gene segment file:\""<<dSegmentFileName<<"\"\n";
  cout<<"\tJ gene segment file:\""<<jSegmentFileName<<"\"\n";
  if((unsigned)numOfSeqSave<=0)
    cout<<"\tfirst nth sequence to rewrite to file:\""<<numOfSeqSave<<"\"";

  //cout<<"\tThe error cost: "<<errorCost<<"\n";
  //cout<<"\tThe number of sequences processed for each file: "<<numberOfSeqProcessedEach<<"\n";
  //start testing
  
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

  //now saving into readable format
  //save the header first
  //

  ofs_ss<<"ID\t"<<"seqName\t"<<"V_gen_1\t"<<"V_gen_1_index\t"<<"V_gen_2\t"<<"V_gen_2_index\t"<<"V_gen_3\t"<<"V_gen_2_index\t"
	                     <<"D_gen_1\t"<<"D_gen_1_index\t"<<"D_gen_2\t"<<"D_gen_2_index\t"<<"D_gen_3\t"<<"D_gen_3_index\t"
	                     <<"J_gen_1\t"<<"J_gen_1_index\t"<<"J_gen_2\t"<<"J_gen_2_index\t"<<"J_gen_3\t"<<"J_gen_3_index\t"
	<<endl;
  cout<<"doing rewriting......."<<endl;
  for(unsigned i=0;i<total_alignment;i++)
    {
      cout<<"."<<i<<":"<<seq[i].GetName()<<endl;
      ofs_ss<<i<<"\t"
	    <<seq[i].GetName()<<"\t";
      //v alignment
      cout<<"v genes:"<<endl;
      if((unsigned(v_align[i].numOfAligned))>0)
	{
	  ofs_ss<<genV[v_align[i].alleles_all[0]].Get_Name()<<"\t";
	  ofs_ss<<genV[v_align[i].alleles_all[0]].Get_GeneIndex()<<"\t";
	}
      else
	ofs_ss<<"\t\t";
      
      if((unsigned(v_align[i].numOfAligned))>1)
	{
	  ofs_ss<<genV[v_align[i].alleles_all[1]].Get_Name()<<"\t";
	  ofs_ss<<genV[v_align[i].alleles_all[1]].Get_GeneIndex()<<"\t";
	}
      else
	ofs_ss<<"\t\t";
      
      if((unsigned(v_align[i].numOfAligned))>2)
	{
	  ofs_ss<<genV[v_align[i].alleles_all[2]].Get_Name()<<"\t";
	  ofs_ss<<genV[v_align[i].alleles_all[2]].Get_GeneIndex()<<"\t";
	}
      else
	ofs_ss<<"\t\t";

      //d alignment
      cout<<"d genes:"<<endl;
      //if((unsigned(d_align[i].allele_order[))>0)
      ofs_ss<<genD[d_align[i].allele_order[0]].Get_Name()<<"\t";
      ofs_ss<<genD[d_align[i].allele_order[0]].Get_GeneIndex()<<"\t";
      ofs_ss<<genD[d_align[i].allele_order[1]].Get_Name()<<"\t";
      ofs_ss<<genD[d_align[i].allele_order[1]].Get_GeneIndex()<<"\t";
      ofs_ss<<genD[d_align[i].allele_order[2]].Get_Name()<<"\t";
      ofs_ss<<genD[d_align[i].allele_order[2]].Get_GeneIndex()<<"\t";
      //else
      //ofs_ss<<"\t";

      //j_alignment
      cout<<"j genes:"<<endl;
      if((unsigned(j_align[i].numOfAligned))>0)
	{
	  ofs_ss<<genJ[j_align[i].alleles_all[0]].Get_Name()<<"\t";
	  ofs_ss<<genJ[j_align[i].alleles_all[0]].Get_GeneIndex()<<"\t";
	  cout<<"\t"<<genJ[j_align[i].alleles_all[0]].Get_Name()<<"\n";
	}
      else
	{
	  ofs_ss<<"\t\t";
	  cout<<"empty\t"<<endl;
	}

      if((unsigned(j_align[i].numOfAligned))>1)
	{
	  ofs_ss<<genJ[j_align[i].alleles_all[1]].Get_Name()<<"\t";
	  ofs_ss<<genJ[j_align[i].alleles_all[1]].Get_GeneIndex()<<"\t";
	  cout<<"\t"<<genJ[j_align[i].alleles_all[1]].Get_Name()<<"\n";
	}
      else
	{
	  cout<<"\tempty"<<endl;
	  ofs_ss<<"\t\t";
	}

      if((unsigned(j_align[i].numOfAligned))>2)
	{
	  ofs_ss<<genJ[j_align[i].alleles_all[2]].Get_Name()<<"\t";
	  ofs_ss<<genJ[j_align[i].alleles_all[2]].Get_GeneIndex();
	}
      else
	ofs_ss<<"\t";

      ofs_ss<<"\n";
      cout<<";"<<endl;
    }

  
  //ofstream ofs_rw;

  //now save the first nth aligned seqs to a different file
  if((unsigned)numOfSeqSave>0|| (unsigned)numOfSeqSave<total_alignment)
    {
      cout<<"Saving the first nth alignment seqs to file......"<<endl;
      //if we are requesting to do a reasonable number of seqs.
      /*
      ofs_rw.open(outFileNam+"_firstSeq.align", std::ios::out);
      if(!ofs_rw.is_open())
	{
	  cout<<"********Error: can not open out file, quit......"<<endl;
	  exit(-1);
	  }*/
      //now start doing the rwriteing
      if(!DoSerialization(v_align, d_align, j_align, seq, total_alignment, alignmentFileName+"_firstSeqs.align", numOfSeqSave))
	{
	  cout<<"ERror in selecting first nth sequence for saving. quit......"<<endl;
	  cerr<<"ERror in selecting first nth sequence for saving. quit......"<<endl;
	  exit(-1);
	}
	
      //done
    }
  
  cout<<"Done..!!"<<endl;

  //clean up
  delete[] genV;
  delete[] genD;
  delete[] genJ;
  delete[] v_align;
  delete[] d_align;
  delete[] j_align;
  delete[] seq;
  ofs_ss.close();

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
	  
	case 'o':
	  outFileName=optarg;
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

	case 'n':
	  numOfSeqSave=atoi(optarg);
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
	
	  //case '':  
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
  cout<<"\t"<<argv[0]<<"-s sequence inputfile \n"
      <<"\t\t [-o output filename]\n"
      <<"\t\t [-v/d/j v/d/j segment file]\n"
      <<"\t\t [-n #]"
    ;

  cout<<"\t\t-s filename -- the sequence file as the input \n"
      <<"\n";

  cout<<"\t\t-o filename -- the output file name\n"
      <<"\n";

  cout<<"\t\t-v/d/j filename -- the v/d/j gene segment/template file name\n";
  
  cout<<"\t\t-n # -- the number of first aligned sequences to save to another file\n";

  cout<<"\t\t-h -- help\n";
  cout<<"\n\t\tv0.1 this code is used to rewrite the alignment file to text format\n"
      <<"\t\t\t so that the input sequence file can be read/used by human,\n"
      <<"\t\t\t it also try to save the first n aligned seqs to another file when -n is specified.\n"
      <<"\t\t\t .\n" 
    <<"\t\t\t @3/20/2015\n"
     ;

  //cout<<"\n\t**************************"<<endl;
  cout<<"\t\t\t*********by Feng @ BU 2015\n"; 

  exit(-1);
}

