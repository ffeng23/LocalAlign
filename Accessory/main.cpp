#include <iostream>
#include <string.h>
#include "Poisson.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "FileHandler.hpp"
#include "FASTQ.hpp"
#include "FastaHandler.hpp"
#include "FastqHandler.hpp"

using namespace std;
//string ifname("Undetermined_S0_L001_R1_001.fastq.gz");
string ifname {"Undetermined_S0_L001_R1_001_first5000.fastq.gz"};
int main(int argc, char* argv[])
{
	if(argc>=2)
	{
		ifname.assign(argv[1]);
	}
  cout<<"****Testing Poisson distribution******"<<endl;

  Poisson pois;

  //now calling the function
  double p=pois.Density(10, 500);
  cout<<"calling pois.Density(10,11):"<<p<<endl;
  cout<<"***DONE*****"<<endl;

  char c[1000];
  cout<<"testing the uninitialized character array:"<<endl;
  cout<<"the string is: "<<c<<endl;
  cout<<"the len is :"<<strlen(c)<<endl;
  cout<<"the char is:"<<(int)(c[0])<<endl;

  cout<<"Testing memcpy........."<<endl;
  unsigned dst[]{1,2,3,4,5,7,8,9,10};
  cout<<dst[1]<<endl;
  cout<<"Copy/shift 2-3 to 3-4......"<<endl;
  memcpy(&(dst[1]), &(dst[2]), sizeof(unsigned)*2);
	
  cout<<"[0]:"<<dst[0]<<";[1]:"<<dst[1]<<";[2]:"<<dst[2]<<";[3]:"<<dst[3]<<";[4]:"<<dst[4]<<endl;
  
  GZ_UTILITY gu{0};
  gu.strm={0};
  gu.z_stream_init=false;
  cout<<"gu."<<gu.z_stream_init<<endl;
  init_gu(gu);
  /*GZ_UTILITY2 gu2{0};
  int x=10;
  cout<<"gu2===="<<gu2.z_stream_init<<endl;
  
  gu2.gzip_in[0]='c';
  if(gu2.gzip_in==nullptr)
	  cout<<"null poointer.........."<<endl;
  cout<<"x address:"<<&x<<"; c:"<<&(gu2.gzip_in[0])<<endl;
  cout<<"gu array address:"<<(char*)(gu2.gzip_in)<<endl;
  printf("address by printf:%p\n", gu2.gzip_in);
//    printf("SIZE_MAX = %zu\n", SIZE_MAX);
  */
  
  //now testing reading the files, 
  cout<<"Testing reading fasta/fasta.gz/fastq/fastq.gz files......."<<endl;
  vector<SequenceString> vec;
  unsigned num=readFile2SeqStrVector(ifname, vec);
  cout<<"finished....."<<endl;
  cout<<"number of sequences read in:"<<num<<endl;
  for(unsigned i=0;i<num;i++)
  {
	  if(i>20)
		  break;
	  cout<<"i:"<<i<<vec.at(i).toString()<<endl;
  }
  
  cout<<"Done.........."<<endl;
  
  //now testing the gz file, 
  //read the first a few byte
  FILE* f=fopen(ifname.c_str(),"rb" );
  unsigned char c_byte='c';
  fread(&c_byte, sizeof(char), (size_t)1, f);
  
  cout<<"the first byte is "<<(int)c_byte<<endl;
  fread(&c_byte,  sizeof(char), (size_t)1, f);
  cout<<"the second byte is "<<(int)c_byte<<endl;
  fclose(f);
  //
  f=fopen("WL02Nano_barcode1_noPhix.fasta", "rb");
  fread(&c_byte, sizeof(char), (size_t)1, f);
  
  cout<<"the first byte of fasta file :"<<c_byte<<endl;
  fclose(f);
  
  //
  f=fopen("SCRS-01-1_S1_L002_R1_001_first100.fastq", "rb");
  fread(&c_byte, sizeof(char), (size_t)1, f);
  
  cout<<"the first byte of fastq file : "<<c_byte<<endl;
  fclose(f);
  
  //now testing the read filetype
  cout<<"------------start testing check file name-------------"<<endl;
  string fn("Undetermined_S0_L001_R1_001_first5000_2618First.fastq");
  cout<<"is_text():"<<fn<<"--"<<is_text(fn)<<endl;
  cout<<"is_fastq():"<<fn<<"--"<<is_fastq(fn)<<endl;
  cout<<"is_fasta():"<<fn<<"--"<<is_fasta(fn)<<endl;
  FileType ft=getFileType_deep(fn);
  cout<<"check filetype deep():"<<fn<<"----"<<ft<<endl;
  
  ft=check_gzFileType(fn);
  
  cout<<"check_filetype():"<<fn<<"--"<<(ft==UNKNOWN)<<endl;  
  
  fn=move("Undetermined_S0_L001_R1_001.fastq_bar0.fasta");
  cout<<"is_text():"<<fn<<"--"<<is_text(fn)<<endl;
  cout<<"is_fastq():"<<fn<<"--"<<is_fastq(fn)<<endl;
  cout<<"is_fasta():"<<fn<<"--"<<is_fasta(fn)<<endl;
  ft=getFileType_deep(fn);
  cout<<"check filetype deep():"<<fn<<"----"<<ft<<endl;
  
  fn=move("Undetermined_S0_L001_R1_001_first5000_2618First.fastq.gz");
  cout<<"is_text():"<<fn<<"--"<<is_text(fn)<<endl;
  cout<<"is_fastq():"<<fn<<"--"<<is_fastq(fn)<<endl;
  cout<<"is_fasta():"<<fn<<"--"<<is_fasta(fn)<<endl;
  ft=getFileType_deep(fn);
  cout<<"check filetype deep():"<<fn<<"----"<<ft<<endl;
  
  ft=check_gzFileType(fn);
  
  cout<<"check_filetype():"<<fn<<"--"<<(ft==GZ)<<endl;
    cout<<"check_filetype():"<<fn<<"--"<<(ft==GZ_FASTA)<<endl;
	cout<<"check_filetype():"<<fn<<"--"<<(ft==GZ_FASTQ)<<endl;
	
  //fn=move("Undetermined_S0_L001_R1_001.fastq_bar0.fasta.gz");
    fn=move("Undetermined_S0_L001_R1_001.fastq_bar1.fasta");
  cout<<"is_text():"<<fn<<"--2"<<is_text(fn)<<endl;
  cout<<"is_fastq():"<<fn<<"--2"<<is_fastq(fn)<<endl;
  cout<<"is_fasta():"<<fn<<"--2"<<is_fasta(fn)<<endl;
  ft=check_gzFileType(fn);
  
  cout<<"check_filetype():"<<fn<<"--2"<<(ft==GZ)<<endl;
    cout<<"check_filetype():"<<fn<<"--2"<<(ft==GZ_FASTA)<<endl;
	cout<<"check_filetype():"<<fn<<"--2"<<(ft==GZ_FASTQ)<<endl;	
  
  //call the outer function.
  ft=getFileType_deep(fn);
  cout<<"check filetype deep():"<<fn<<"----"<<ft<<endl;
  
  ft=getFileType(fn,false);
  cout<<"check filetype():"<<fn<<"********"<<ft<<endl;
  ft=getFileType(fn,true);
  cout<<"check filetype by name():"<<fn<<"********"<<ft<<endl;
  
  //start testing the concatenation.
  cout<<"-----------------start testing the concatenation-----------"<<endl;
  //open file 
  //GZ_UTILITY gu2;
	FILE* fb {nullptr};
	//fb=gzOpen_B(fn, gu, "rb");
	fb=fopen(fn.c_str(), "r");
	if(!fb)
	{
		cout<<"can not open file:"<<fn<<endl;
	}
  //now let's read one record.
  bool ok=true;
  SequenceString ss;
  int lc=0;
  string hangover;
  while (ok)
  {
	  ok=get1FastaSeq(fb,ss , hangover, gu, false);
	  if(ok||(ss.GetName().length()!=0&&ss.GetSequence().length()!=0)){
		lc++;
	  
		//cout<<"record "<<lc<<":"<<ss.toString()<<endl;//if(!ok)
	  }
	  /*else
	  {
		  cout<<"breaking......."<<endl;
		  break;
	  }*/
	 
  }
  cout<<"total number of reads:"<<lc<<"........."<<endl;
  //gzClose_B(fb,gu);
  fclose(fb);
  
  //---fasta multi-line fasta
  cout<<is_fasta("top_1_clone1073_seqs.picked.fa.IgG3_11nt_pre.fa")<<endl;
  cout<<getFileType("top_1_clone1073_seqs.picked.fa.IgG3_11nt_pre.fa",true)<<endl;
  
  //-----------testing the read1FastqSeq--------
  cout<<"++++++++++++++++testing read single fastq sequence+++++++++"<<endl;
  fn.assign("SCRS-01-1_S1_L002_R1_001_first98a.fastq");
//  fb=gzOpen_B(fn, gu, "rb");
  	fb=fopen(fn.c_str(), "r");
	if(!fb)
	{
		cout<<"can not open file:"<<fn<<endl;
	}
  //now let's read one record.
  ok=true;
  cout<<"left over:"<<ss.toString()<<endl;
  //SequenceString ss;
  lc=0;
  string qs;
  while (ok)
  {
	  ok=get1FastqSeq(fb,ss , qs, gu, false);
	  if(ok||(ss.GetName().length()!=0&&ss.GetSequence().length()!=0)){
		lc++;
	  
		cout<<"record "<<lc<<":"<<ss.toString();//if(!ok)
			cout<<"+\n"<<qs<<endl;
	  }
	  /*else
	  {
		  cout<<"breaking......."<<endl;
		  break;
	  }*/
	 
  }
  cout<<"total number of reads:"<<lc<<"........."<<endl;
//  gzClose_B(fb,gu);
  fclose(fb);
  
  return 0;
}
