#include <iostream>
//#include <fstream>
//#include <vector>
//#include "string_ext.hpp"
//#include "read_gene_info_func.hpp"

//from on you need to specify the c libraries explicitly, 11/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
//#include "score.hpp"
//#include "SequenceString.hpp"
//#include "LocalAlignment.hpp"
//#include "GlobalAlignment.hpp"
//#include "OverlapAlignment.hpp"
//#include "FastaHandler.hpp"
//#include "SequenceHandler.hpp"
#include "Matrix.hpp"
using namespace std;

static void printUsage(int argc, char* argv[]);
static void parseArguments(int argc, char **argv, const char *opts);
//static int lookUpScoreMatrix(const string* _scoreMatrixNameArray,const int& len, const string& scoreMatrixName);

//static string inputFile1_name; //the input file for seq1
//static string inputFile2_name;//input ifle for seq2  
//both above files are in fasta file and should contains the same sets of genes.

//static string outputFile_name("localalignment.txt"); //the output file for the local alignment results with scores

//static string scoreMatrixName="nuc44"; //name as an input for specifing the name of the score matrix. for example nuc44
/*static enum SeqType
  {
    AminoAcid,
    Nucleotide
  } seqtype=Nucleotide; //by default
  static string supportedScoreMatrixNameArr[]={"nuc44","blosum50", "tsm1", "tsm2", "nuc44HP"};*/
//tsm2: score matrix from Geoffrey Barton reference.
/*
static ScoreMatrix* ScoreMatrixArr[]={&nuc44, &blosum50, &tsm1, &tsm2, &nuc44HP};

static double scale=1; //this is the one on top of matrix, the programe will run score and use the 
//matrix specified scale first and then apply the scale set by this one.

double gapopen=-8;
double gapextension=-8;
bool gapextensionFlag=false;
//static int read_gene_sequence(const string& temp_seq, vector<string>& promoter_sequence);

*/
int main(int argc, char* argv[])
{
  //for parsing commandline arguement

  const char *opts = "hvg:e:o:s:t:k:am:";
  //f: the input file name containing the fasta format gene sequence
  //c: flag to indicate to count each gene length.
  //o: output file name

  parseArguments(argc, argv, opts);
  
  /*if(inputFile1_name.size()==0)
    {
      cout<<"please specify the seq1 input fasta file name.......\n";
      cout<<"type \"./testalign -h\" for usage help\n";
      exit(-1);
      //printUsage(argc, argv);
    }
  if(inputFile2_name.size()==0)
    {
      cout<<"please specify the seq2 input fasta file name.......\n";
      //printUsage(argc, argv);
      cout<<"type \"./testalign -h \" for usage help\n";
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
  */
  
  cout<<"start testing matrix function........."<<endl;
  cout<<"first building the matrix"<<endl;
  Matrix<unsigned> m1;
  cout<<"m1:"<<m1.toString()<<endl;

  //now we try to initialize it
  unsigned dim_size[]={2,3,2};
  unsigned data[]={1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12};
  m1.initialize(3, dim_size, data);
  m1=m1-1;
  cout<<"after initializeation"<<endl;
  cout<<"m1:\n"<<m1.toString()<<endl;
  m1(0,0,0)+=100;
  cout<<"m1:\n"<<m1.toString()<<endl;

  //create another one
  double data1[]={1,2,3,4,5,6,7,8,9,10,11,12,1,2,3,4,5,6,7,8,9,10,11,12};
  Matrix<double> m2(3, dim_size, data1);
  m2=m2+2;
  cout<<"m2:\n"<<m2.toString()<<endl;


  Matrix<double> m3=m2-12;
  cout<<"m3:\n"<<m3.toString()<<endl;
  
  //size
  unsigned s=m3.size(1);
  cout<<"m3 size of dim 1:"<<s<<endl;
  cout<<"m3 size :"<<m3.size().toString()<<endl;

  //testing submatrix
  cout<<"testing submatrix:"<<endl;
  int sub[]={1,-1,-1};
  cout<<m3.SubMatrix(3, sub).toString()<<endl;
  
  int sub2[]={-1,-1,0};
  cout<<m3.SubMatrix(3,sub2).toString()<<endl;

  cout<<"===>testing a 4D matrix:"<<endl;
  //start a new test with 4D
  double data2[]={ 1, 2, 3, 4, 4.7,
		   5, 6, 7, 8, 8.7,
		   9,10,11,12, 12.7,
		   9.8,10.8, 11.8,12.8,12.85,
		  
		   -1,-2,-3,-4, -4.7,
		   -5,-6,-7,-8,-8.7,
		   -9,-10,-11,-12, -12.7,
		   -9.5,-19.8, -11.8,-12.8, -12.85,
		   /*4D-0*/
		   1.1, 2.1, 3.1, 4.1, 4.17,
		   5.1, 6.1, 7.1, 8.1, 8.17,
		   9.1,10.1,11.1,12.1, 12.17,
		   9.18, 10.18, 11.18, 12.18, 12.175,
		  
		   -1.1,-2.1,-3.1,-4.1, -4.17,
		   -5.1,-6.1,-7.1,-8.1, -8.17,
		   -9.1,-10.1,-11.1,-12.1, -12.17,
		   -9.18, -10.18, -11.18, -12.18, -12.175,

		   /*4D-1*/
		   1.2, 2.2, 3.2, 4.2, 4.27,
		   5.2, 6.2, 7.2, 8.2, 8.27,
		   9.2,10.2,11.2,12.2, 12.27,
		   9.28, 10.28, 11.28, 12.28, 12.275,
		  
		   -1.2,-2.2,-3.2,-4.2, -4.27,
		   -5.2,-6.2,-7.2,-8.2, -8.27,
		   -9.2,-10.2,-11.2,-12.2, -12.27,
		   -9.28, -10.28, -11.28, -12.28, -12.275
		   /*4D-2*/
		     };
  unsigned dim_size4[]={3,2,4,5};
  Matrix<double> m4(4, dim_size4, data2);
  cout<<m4.toString()<<endl;

  //start testing the submatrix on 4d matrix
  int sub4[]={1,-1,-1,3};
  cout<<m4.SubMatrix(4,sub4).toString()<<endl;
  cout<<"Done!!!"<<endl;

  cout<<m4.toString()<<endl;
  //start doing the su
  cout<<sum(m4).toString()<<endl;
  //ofs.close();

  //cout<<"Total "<<gene_info.size()<<" genes are processed"<<endl;

  //start testing the divide by dimension
  unsigned dim_size_vec[]={2};
  Matrix<double> vec(1, dim_size_vec, 10);
  vec(1)=100;
  //vec(2)=1000;
  //vec(3)=1;
  cout<<vec.toString()<<endl;
  cout<<m2.toString()<<endl;
  m2.divide_by_dimension(vec, 0);
  cout<<m2.toString()<<endl;

  //testing addition
  cout<<"testing dot additon of two matrices"<<endl;
  cout<<"===>m1:\n"<<m1.toString()<<endl;
  cout<<"===>m2:\n"<<m2.toString()<<endl;
  m2=(m2+m2);
  cout<<m2.toString()<<endl;

  //testing accessor by array
  unsigned ins[]={0,0,0,0};
  cout<<"==================="<<endl;
  cout<<m4.toString()<<endl;
  cout<<m4(ins, 4)<<endl;

  cout<<"changing the element "<<endl;
  m4(ins,4)=1000;
  cout<<m4.toString()<<endl;
  //
  cout<<"testing scalar matrix, dim=0"<<endl;
  Matrix<unsigned> m0(0, NULL, 15);
  cout<<m0.toString()<<endl;
  cout<<m0(ins, 0)<<endl;

  //testing a nan
  cout<<"====>testing nan"<<endl;
  double x=0;
  double y=0/x;
  cout<<y<<endl;
  if(isnan(y))
    {
      cout<<"y is a nan"<<endl;
    }
  //testing max()
  cout<<"=============>testing max()";
  cout<<m4.toString()<<endl;
  cout<<max(m4,3).toString()<<endl;

  //testing m2vec()
  cout<<"------testing m2vec----"<<endl;
  Matrix<double> mvec=m4.m2vec();
  cout<<mvec.toString()<<endl;

  cout<<"-====-----m4"<<endl;
  cout<<m4.toString()<<endl;

  //testing comparison and bitwise operation
  cout<<"666666666666testing bitwise operation"<<endl;
  Matrix<bool> m4Big0(m4>0);
  Matrix<bool> m4Big10(m4>10);
  cout<<"m4>0"<<m4Big0.toString()<<endl;
  cout<<"m4>10"<<m4Big10.toString()<<endl;
  cout<<"m4>10&m4>0"<<(m4Big10&m4Big0).toString()<<endl;

  cout<<"m4 >0 elements"<<m4.GetElements(m4Big0&(m4<0)).toString()<<endl;
  cout<<"m4:\n"<<m4.toString()<<endl;
  cout<<"Thanks for using our program and have a nice day!!"<<endl;

  int dimSub[]={0,-1,-1,-1};
  Matrix<double> m4Sub=m4.SubMatrix(4, dimSub);
  cout<<"m4Sub (0):"<<m4Sub.toString()<<endl;

  m4.SetSubMatrix(1, m4Sub);
  cout<<"m4:\n"<<m4.toString()<<endl;
  dimSub[2]=1;
  Matrix<double> m4Sub2D=m4.SubMatrix(4, dimSub);
  cout<<"m4Sub2D (0,1):"<<m4Sub2D.toString()<<endl;

  double inputArr[]={0,0,1,1.1,-1000};
  m4Sub2D.SetSubMatrix(1,5,inputArr);
  cout<<"m4Sub2D (0,1):"<<m4Sub2D.toString()<<endl;
  
  //1d multiplication
  dim_size_vec[0]=2;
  Matrix<double> multi_m1(1,dim_size_vec, 1);
  Matrix<unsigned> multi_m2(1, dim_size_vec, 2);
  cout<<"1d multiplication:"<<matrix_multiply_1D<unsigned>(multi_m1, multi_m2)<<endl;
  
  Matrix<double> m4SubSum=sum(m4Sub2D,0);
  cout<<"sum:"<<m4SubSum.toString()<<endl;

  cout<<"m4Sub: this is 3D\n"<<m4Sub.toString()<<endl;
  double zeroArr[]={0,0,0,0,0};
  m4Sub.SetSubMatrix(1,1,5,zeroArr);
  cout<<"m4Sub: after setting zeros:"<<m4Sub.toString()<<endl;

  cout<<"....=====>comparing signed vs. unsigned integer"<<endl;
  signed sx=5; unsigned sy=-1;
  cout<<"comparing 5 unsigned vs -1 signed :"<< ((signed)(sx-1) > (signed)(sy-1)) <<endl;

  cout<<"unsigned int 5 subtracting negative number -1:"<<sx-sy<<endl;

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


