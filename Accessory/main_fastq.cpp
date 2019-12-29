#include <iostream>
#include <string.h>

#include "FASTQ.hpp"
#include "FastqHandler.hpp"
#include "../SequenceHandlerBarcode.hpp"
#include "FastaHandler.hpp"
using namespace std;

string ifname("/home/feng/Feng/hg/LocalAlign/Accessory/SCRS-01-1_S1_L002_R1_001_first100.fastq.gz");
int main(int argc, char* argv[])
{
  if(argc ==2)
	  ifname.assign(argv[1]);

  cout<<"****Testing Poisson distribution******"<<endl;
//string c("hello world");
  char c[1000];
  cout<<"the string is: "<<c<<endl;
  cout<<"the len is :"<<strlen(c)<<endl;
  cout<<"the char is:"<<(int)(c[0])<<endl;
  
  //start checking
  cout<<"Start checking the fastq class......."<<endl;
  SequenceString ss("seq1 ", "ATCG");
  cout<<ss.toString()<<endl;
  
  Fastq fq("seq 1",ss, "@+AF");
  cout<<fq.toString(true)<<endl;
  
  ofstream ofs("text.txt");
  fq.Serialize(ofs);
  ofs.close();
  
  ifstream ifs("text.txt");
  Fastq fq_d;
  cout<<"before deserialization:"<<endl;
  cout<<"fq_d:"<<fq_d.toString()<<endl;
  fq_d.Deserialize(ifs);
  ifs.close();
  
  cout<<"after deserialization:"<<endl;
  cout<<"fq_d:"<<fq_d.toString(true)<<endl;
  string _fname("feng.tar.fastq.gz1");
  bool gzflag=_fname.substr(_fname.size()-3, 3)==".gz";
  cout<<"file name input:"<<_fname<<endl;
  if(gzflag)
	  cout<<"\tYes, found it correctly"<<endl;
  else
	  cout<<"\tNo, did not find it."<<endl;
  
  //Start testing the readfastq function
  cout<<"---------start testing read fastq (either gz'ed or regular fastq) ......."<<endl;
  vector<Fastq> vec;
  _fname.assign(ifname);
  
  unsigned int num=ReadFastq(_fname, vec, false);
  cout<<"total number of sequences read: "<<num<<endl;
  
  vector<SequenceString> v_seq;
  //now we need to print out for debugging.
  for(unsigned i=0;i<vec.size();i++)
  {
	  cout<<"i:"<<i<<"--"<<vec.at(i).toString()<<endl;
	  v_seq.push_back(vec.at(i).GetSequenceString());
  }
  cout<<"DONE"<<endl;
  cout<<"sequence string 0:"<<v_seq.at(0).toString()<<endl;
  //
  cout<<"=====Start testing get barcodes ..........."<<endl;
  
  
  //first get the indexes from the header of the sequences.
  vector<SequenceString> v_index1;
  vector<SequenceString> v_index2;
  
  unsigned lenOfBarcode=8;
  bool dualIndex=true;
  
  cout<<"Reading indexes from Sequence names......."<<endl;
  unsigned numOfSeq=ReadIndexFromSequenceName(v_seq, v_index1, v_index2, lenOfBarcode, dualIndex);

  cout<<"\tnumber of sequences processed:"<<numOfSeq<<endl;
  
  /*for(unsigned i=0;i<v_index1.size();i++)
  {
	  cout<<v_index1.at(i).toString()<<endl;
	  cout<<v_index2.at(i).toString()<<endl;
  }*/
  //now we have the indexes, pick the unique ones as barcodes from them
  vector<SequenceString> v_BarSeq1;
  vector<SequenceString> v_BarSeq2;
  vector<unsigned> v_count;
  cout<<"Reading barcodes from indexes......."<<endl;
  unsigned numOfBar=GetBarcodes(v_index1, v_index2, dualIndex, v_BarSeq1, v_BarSeq2, v_count);
  cout<<"\tnumber of bar:"<<numOfBar<<endl;
  
  WriteFasta("index1.fasta", v_index1);
  WriteFasta("index2.fasta", v_index2);
  WriteFasta("bar1.fasta", v_BarSeq1);
  WriteFasta("bar2.fasta", v_BarSeq2);
  WriteTextFile("count.txt", v_count);
  
  return 0;
}
