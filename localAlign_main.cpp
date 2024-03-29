#include <iostream>
#include <fstream>
#include <vector>
#include "Accessory/string_ext.hpp"
//#include "read_gene_info_func.hpp"

//from on you need to specify the c libraries explicitly, 11/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "score.hpp"
#include "Accessory/SequenceString.hpp"
#include "LocalAlignment.hpp"
#include "GlobalAlignment.hpp"
#include "OverlapAlignment.hpp"
#include "Accessory/FastaHandler.hpp"
#include "SequenceHandler.hpp"
#include "LocalAlignment_CT.hpp"

using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);
static int lookUpScoreMatrix(const string* _scoreMatrixNameArray,const int& len, const string& scoreMatrixName);

static string inputFile1_name; //the input file for seq1
static string inputFile2_name;//input ifle for seq2  
//both above files are in fasta file and should contains the same sets of genes.

static string outputFile_name("localalignment.txt"); //the output file for the local alignment results with scores

static string scoreMatrixName="nuc44"; //name as an input for specifing the name of the score matrix. for example nuc44
static enum SeqType
  {
    AminoAcid,
    Nucleotide
  } seqtype=Nucleotide; //by default
static string supportedScoreMatrixNameArr[]={"nuc44","blosum50", "tsm1", "tsm2", "nuc44HP", "nuc44DM1"};
//tsm2: score matrix from Geoffrey Barton reference.

static ScoreMatrix* ScoreMatrixArr[]={&nuc44, &blosum50, &tsm1, &tsm2, &nuc44HP, &nuc44DM1};

static double scale=1; //this is the one on top of matrix, the programe will run score and use the 
//matrix specified scale first and then apply the scale set by this one.

double gapopen=-8;
double gapextension=-8;
bool gapextensionFlag=false;
//static int read_gene_sequence(const string& temp_seq, vector<string>& promoter_sequence);


int main(int argc, char* argv[])
{
  //for parsing commandline arguement

  const char *opts = "hvg:e:o:s:t:k:am:";
  //f: the input file name containing the fasta format gene sequence
  //c: flag to indicate to count each gene length.
  //o: output file name

  parseArguments(argc, argv, opts);
  
  if(inputFile1_name.size()==0)
    {
      cout<<"please specify the seq1 input fasta file name.......\n";
      cout<<"type \"./"<<argv[0]<<" -h\" for usage help\n";
      exit(-1);
      //printUsage(argc, argv);
    }
  if(inputFile2_name.size()==0)
    {
      cout<<"please specify the seq2 input fasta file name.......\n";
      //printUsage(argc, argv);
      cout<<"type \"./"<<argv[0]<<" -h \" for usage help\n";
      exit(-1);
    }

  cout<<"***Input parameters Summary:\n";
  cout<<"\tSeq1 file name:\""<<inputFile1_name<<"\".\n";
  cout<<"\tSeq2 file name:\""<<inputFile2_name<<"\".\n";
  cout<<"\tThe output file name : \""<<outputFile_name<<"\".\n";
  if(seqtype==AminoAcid)
    cout<<"\tSequency Type:Amino Acids\n";
  else
    cout<<"\tSequency Type:Nucleotides\n";

  //look up the score matix
  int scoreMatrixIndex=lookUpScoreMatrix(supportedScoreMatrixNameArr, sizeof(supportedScoreMatrixNameArr)/sizeof(supportedScoreMatrixNameArr[0]),scoreMatrixName);
  if(scoreMatrixIndex==-1)
    {
      cout<<"\tscore matrix specified by input was not found(\""<<scoreMatrixName<<"\") . "
	  <<"\n\tUsing the default scorematrix.\n";
      scoreMatrixName="nuc44";
      if(seqtype==AminoAcid)
	{
	  scoreMatrixName="blosum50";
	}
      scoreMatrixIndex=lookUpScoreMatrix(supportedScoreMatrixNameArr, sizeof(supportedScoreMatrixNameArr)/sizeof(supportedScoreMatrixNameArr[0]), scoreMatrixName);
    }
  cout<<"\tscore matrix:"<<scoreMatrixName<<"\n";
  
  cout<<"\tscale matrix scale:"<<scale<<"\n";
  
  cout<<"\tgap open penalty:"<<gapopen<<"\n";
  cout<<"\tgap extension penalty:"<<gapextension<<"\n";
  cout<<"  ****************\n";

  cout<<"Testing ScoreMatrix:\n";
  cout<<"\tthe index: "<<scoreMatrixIndex<<endl;
  ScoreMatrix* sm= ScoreMatrixArr[scoreMatrixIndex];
  char c1='C', c2='C';

  cout<<"\t("<<c1<<","<<c2<<")="<<sm->GetScore(c1,c2)<<endl;
  
  //fasta handler
  cout<<"Reading input file 1:";
  vector<SequenceString> vec_seq1;
  int num1=ReadFasta(inputFile1_name, vec_seq1);
  if(num1<0)
       exit(-1);
  cout<<num1<<" sequence(s) read in"<<endl;
	if(num1>1)
	  {
		  cout<<"\tmore than one sequences read in from the input file,"
			  <<"\t\tonly the first one will be used"<<endl; 		  
	  }
	if(num1==0)
	{
		cout<<"ERROR***:no sequence read from the input, please check the file...."<<endl;
		exit(-1);
	}
	SequenceString Seq1=vec_seq1.at(0);//("seq1", "ATGACGGTGATGATGGTAGTGGCCATTTCTTTGCCTACCCACTGTGGCCTTGTTGTGAAGAACTTTGATGCCCTGTTGTATGTGATCGTTATGACAATCCTGTGAAGTTTTTCAATGATGAAGGAAACTGAGGCTTAGAGAGTTGGAGTAATCTGTGAAGGCTCAGGATGGGCAAGAGGTCCCGCCCAGGTTTGAGCCCAGATGCGAGGTTACCACGCTTCCTGGTGAGGTGTTTTACAACTAAGGCCAAGCCAGGCAAAACCCATTGTTCTGCAGCTTCTGGCTTGGATTGGGTGTCTTGTTGAGTATGTGGGCAGTGGATCTGATGTTTTCCACTTCCACCAAGGGCCCATCGGTCTTCCCCCTGT");
  
  //fasta handler
  cout<<"Reading input file 2:";
  vector<SequenceString> vec_seq2;
  int num2=ReadFasta(inputFile2_name, vec_seq2);
  cout<<num2<<" sequences(s) read in"<<endl;
  if(num2<0)
      exit(-1);
	if(num2>1)
	  {
		  cout<<"\tmore than one sequences read in from the input file,"
			  <<"\t\tonly the first one will be used"<<endl; 		  
	  }
	 if(num2==0)
	{
		cout<<"ERROR***:no sequence read from the input, please check the file...."<<endl;
		exit(-1);
	}
  SequenceString Seq2=vec_seq2.at(0);//("seq2","GTATTATGATTACGTTTGGGGGAGTTATCGTTATACC");
  cout<<"showing sequence string\n"<<Seq1.toString()<<Seq2.toString()<<endl;

  SequenceString tempSStr=ReverseComplement(Seq2);
  
  //now testing alignment
  cout<<"Running local alignment:"<<endl;

cout<<"\t("<<c1<<","<<c2<<")="<<sm->GetScore('A','A')<<endl;
  
 LocalAlignment la(&Seq1,&Seq2,sm, gapopen, gapextension,scale, 1,0);//,1);
  //LocalAlignment la(&Seq1,&tempSStr,sm, gapopen, gapextension,1, 100);
  cout<<"\tdone and the score is "<<la.GetScore()<<endl;
  cout<<"\t"<<la.GetAlignment().toString()<<endl;
  cout<<"\toffset:"<<la.GetAlignment().GetPatternIndexStart()<<endl;

  //testing multiple alignment
  cout<<"Total number of local aligments:"<<la.GetNumberOfAlignments()<<endl;
  for(unsigned int i=0;i<la.GetNumberOfAlignments();i++)
    {
      cout<<i<<"/"<<la.GetNumberOfAlignments()<<":"<<endl;
      cout<<la.GetAlignmentArr()[i].toString()<<endl;
    }
  
  /*
  if(vec_seq.size()>=0)
	cout<<"1/1000:"<<vec_seq.at(0).toString()<<endl;
  if(vec_seq.size()>=999)
	cout<<"999/1000;"<<vec_seq.at(999).toString()<<endl;

  
  //testing reverseComp
  cout<<"Testing reverse comp:"<<endl;
  SequenceString tempSString=ReverseComplement(vec_seq.at(0));

  cout<<"rev comp 0/1000"<<tempSString.toString()<<endl;
*/
  cout<<"Test substr(\"feng\")"<<endl;
  string tempString("feng");
  cout<<"0,1 : "<<tempString.substr(0,1)<<endl;
  cout<<"0,-100: "<<tempString.substr(0,100)<<endl;
  cout<<"length is "<<tempString.length()<<endl;
  cout<<"4,1 :" <<tempString.substr(4)<<endl;


//
cout<<"---------------Start testing special cases............."<<endl;
    cout<<"see notes in ../UMI_barcode/"<<endl;
     SequenceString p("p","GGGATATATA");
  SequenceString s("s","GGGTATATAT");
  
  cout<<"the pattern string: "<<p.toString()<<endl;
  cout<<"the subject string: "<<s.toString()<<endl;
  LocalAlignment_CT ol(&p, &s, sm, gapopen, gapextension
  						,scale,4,0);//last 2 parameter, # of best local alignment found and and type gapmodel, 0 for affine 
   // LocalAlignment_CT ol(&p, &s, sm, gapopen, gapextension
	//					,scale,4,0);//last 2 parameter, # of best local alignment found and and type gapmodel, 0 for affine                
	AlignmentString* aso_arr=ol.GetAlignmentArr(); //pattern is on top; and subject is on bottom.
	cout<<"first alignment"<<endl;
	cout<<aso_arr[0].toString()<<endl;
	AlignmentString aso=aso_arr[1];
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
  /*ifstream ifs1_p(inputFile1_name.c_str());
  
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

