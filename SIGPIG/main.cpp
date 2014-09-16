#include <iostream>
#include <fstream>
#include <vector>
#include "string_ext.hpp"
//#include "read_gene_info_func.hpp"

//from on you need to specify the c libraries explicitly, 11/
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

#include "../score.hpp"
#include "../SequenceString.hpp"
#include "../LocalAlignment.hpp"
#include "../GlobalAlignment.hpp"
#include "../OverlapAlignment.hpp"
#include "../FastaHandler.hpp"
#include "../SequenceHandler.hpp"
#include "../AlignmentString.hpp"

using namespace std;

static string scoreMatrixName="nuc44"; //name as an input for specifing the name of the score matrix. for example nuc44
/*static enum SeqType
  {
    AminoAcid,
    Nucleotide
  } seqtype=Nucleotide; //by default
*/
static string supportedScoreMatrixNameArr[]={"nuc44","blosum50", "tsm1", "tsm2", "nuc44HP"};
//tsm2: score matrix from Geoffrey Barton reference.

static ScoreMatrix* ScoreMatrixArr[]={&nuc44, &blosum50, &tsm1, &tsm2, &nuc44HP};

//static double scale=1; //this is the one on top of matrix, the programe will run score and use the 
//matrix specified scale first and then apply the scale set by this one.

double gapopen=-8;
double gapextension=-8;
bool gapextensionFlag=false;
//static int read_gene_sequence(const string& temp_seq, vector<string>& promoter_sequence);


int main(int argc, char* argv[])
{
  //testing the alignment string
  ScoreMatrix* sm= ScoreMatrixArr[0];
  char c1='A', c2='T';

  cout<<"\t("<<c1<<","<<c2<<")="<<sm->GetScore(c1,c2)<<endl;
  
  //aminoacid
  
  //SequenceString Seq1("seq1","VSPAGMASGYDPGKA");
  //SequenceString Seq2("seq2", "IPGKATREYDVSPAG");

  //SequenceString Seq1("seq1","ASGYDPGKA");
  //SequenceString Seq2("seq2", "ATREYDVSPAG");


  //testing local alignment
  //SequenceString Seq2("seq1","AGCTAGAGACCCCAGTCTGAGGTAGA");
  //SequenceString Seq1("seq2", "AGCTAGAGACCAGCTATCTAGAGGTAGA");
  //SequenceString Seq1("seq1","ACCCCAG");
  //SequenceString Seq2("seq2", "ACCAG");


  SequenceString Seq1("seq1","ATGAGGTAGA");
  SequenceString Seq2 ("seq2", "CTAGAGGTAGA");

  //****the sequences from the "Geoffrey Barton reference"
  //SequenceString Seq1("seq1","CCAATCTACTACTGCTTGCAGTACTTGT");
  //SequenceString Seq2 ("seq2", "AGTCCGAGGGCTACTCTACTGAAC");

  //SequenceString Seq1("seq1","AGACGCACTCGTTCGGGAAGTAGTCCTTGACCAGGCAGCCCACGGCGCTGTC");
  //SequenceString Seq2 ("seq2", "CGTATCGCCTCCCTCGCGCCATCAGACGAGTGCGTGTTCGGGGAAGTAGTCCTTGAC");
  //SequenceString Seq2("seq2", "CGTATCGCCTCCCTCGCGCCATCAGACGAGTGCGTCAGGAGACGAGGGGGAA");
  //SequenceString Seq2("seq2","CGTATCGCCTCCCTCGCGCCATCAGACGAGTGCGTACGTTGGGTGGTACCCAGTTAT");
  //SequenceString Seq2("seq2", "CGTATCGCCTCCCTCGCGCCATCAGACGAGTGCGTACGACTCACTATAGGGCAAGCAGTGGTAACAACGCAGAGTACGCGGG");
  //SequenceString Seq1("seq1","ATCGA");
  //SequenceString Seq2 ("seq2", "GATTGA");

  //SequenceString Seq2("seq1","CGGGA");
  //SequenceString Seq1 ("seq2", "CGGA");

  //SequenceString Seq1("seq1","ATAT");
  //SequenceString Seq2 ("seq2", "ACHKAT");
  //SequenceString Seq1("seq","AGACGCACTGTTCGGGAAGTAGTCCTTGACCAGGCAGCCACCCATGTACTCTGCGTTGATACCACTGCTTGCCCTATAGTGAGTCGTGAGTGCGTCTCTGACGGGCTGGCAAGGCGCATAG");
  //SequenceString Seq2("Constant","cctccaccaagggcccatcggtcttccccctggcgccctgctccaggagcacctccgagagcacagcggccctgggctgcctggtcaaggactacttccccgaaccggtgacggtgtcgtggaactcaggcgctctgaccagcggcgtgcacaccttcccggctgtcctacagtcctcaggactctactccctcagcagcgtggtgaccgtgacctccagcaacttcggcacccagacctacacctgcaacgtagatcacaagcccagcaacaccaaggtggacaagacagttg");
  
  //SequenceString Seq1("seq1","CCTATAGTGAGTCGTACGCACTCGTCTGAGCGGGCTGGCAAGGCGCATAG");
  //SequenceString Seq2("seq2", "CCCGCGTACTCTGCGTTGTTACCACTGCTTGCCCTATAGTGAGTCGT");

  //SequenceString Seq1("seq1", "AGTGCGTCAGGGACGAGGGGGTGGATTCACCCATGTACTCTGCGTTGATACCACTGCTTGCCCTATAGTGAGTCGTACGCACTCGTCTGANCGGGCTGGCAAGGCGCATAG");
  //SequenceString Seq2("IgMCH1", "CTGGAAGAGGCACGTTCTTTTCTTTGTTGCCGTTGGGGTGCTGGACTTTGCACACCACGTGTTCGTCTGTGCCCTGCATGACGTCCTTGGAAGGCAGCAGCACCTGTGAGGTGGCTGCGTACTTGCCCCCTCTCAGGACTGATGGGAAGCCCCGGGTACTGCTGATGTCAGAGTTGTTCTTGTATTTCCAGGACAAAGTGATGGAGTCGGGAAGGAAGTCCTGTGCGAGGCAGCCAACGGCCACGCTGCTCGTATCCGACGGGGAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCC");
  //SequenceString Seq2("IgGCH1", "CAACTCTCTTGTCCACCTTGGTGTTGCTGGGCTTGTGATTCACGTTGCAGGTGTAGGTCTGGGTGCCCAAGCTGCTGGAGGGCACGGTCACCACGCTGCTGAGGGAGTAGAGTCCTGAGTACTGTAGGACAGCCGGGAAGGTGTGCACGCCGCTGGTCAGGGCGCCTGAGTTCCACGACACCGTCACCGGTTCGGGGAAGTAGTCCTTGACCAGGCAGCCCAGGGCCGCTGTGCCCCCAGAGGTGCTCCTGGAGCAGGGCGCCAGGGGGAAGACCGATGGGCCCTTGGTGGAAG");
  //SequenceString Seq2("IgA1", "CTGGGCAGGGCACAGTCACATCCTGGCTGGGATTCGTGTAGTGCTTCACGTGGCATGTCACGGACTTGCCGGCTAGGCACTGTGTGGCCGGCAGGGTCAGCTGGCTGCTCGTGGTGTACAGGTCCCCGGAGGCATCCTGGCTGGGTGGGAAGTTTCTGGCGGTCACGCCCTGTCCGCTTTCGCTCCAGGTCACACTGAGTGGCTCCTGGGGGAAGAAGCCCTGGACCAGGCAGGCGATGACCACGTTCCCATCTGGCTGGGTGCTGCAGAGGCTCAGCGGGAAGACCTTGGGGCTGGTCGGGGATG");
  //SequenceString Seq2("IgD1", "CTGGCCAGCGGAAGATCTCCTTCTTACTCTTGCTGGCGGTGTGCTGGACCACGCATTTGTACTCGCCTTGGCGCCACTGCTGGAGGGGGGTGGAGAGCTGGCTGCTTGTCATGTAGTAGCTGTCCCGTCTTTGTATCTCAGGGAAGGTTCTCTGGGGCTGGCTCTGTGTCCCCATGTACCAGGTGACAGTCACGGACGTTGGGTGGTACCCAGTTATCAAGCATGCCAGGACCACAGGGCTGTTATCCTTTGGGTGTCTGCACCCTGATATGATGGGGAACACATCCGGAGCCTTGGTGGGTG");

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

  cout<<"Done!!!"<<endl;

  return 0;
}
