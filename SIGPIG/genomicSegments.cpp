#include <vector>
#include <exception>
#include <algorithm>

#include "genomicSegments.hpp"
#include "../FastaHandler.hpp"
#include "GenomicJ.hpp"
#include "GenomicV.hpp"
#include "GenomicD.hpp"
#include "../string_ext.hpp"

bool GenomicVCompare_bySequenceName(const GenomicV& ss1,const GenomicV& ss2)
{
  int ret=ss1.Get_Name().compare(ss2.Get_Name());
  

  if(ret<=0)
    {
      return true;
    }
  else
    return false;
}


unsigned ReadGenomicV(const string& _fastaFileName, GenomicV** _gseg)
{
  //first read the files into vector
  vector<SequenceString> seq;
  unsigned totalNumber=ReadFasta(_fastaFileName, seq, true); 
  //get pointer to array of genomic segments
  *_gseg=new GenomicV[totalNumber];

  GenomicV* p_arr=*_gseg;
  vector<string> gene_vec;
  //cout<<"here it is"<<endl;
  //now start reading into the segments
  for(unsigned int i=0;i<totalNumber;i++)
    {
      string seqName;
      try
	{
	  seqName=ParseSequenceName(seq.at(i).GetName());
	}
      catch (const char* msg)
	{
	  throw msg;
	}
      catch (...)
	{
	  throw "unknow error, quit";
	}
      seq[i].SetName(seqName);//set the sequence string name to short one

      //now, we parse gene name
      string gene;
      try
	{
	  gene=ParseGeneName(seqName);
	}
      catch (const char* msg)
	{
	  throw msg;
	}
      catch (...)
	{
	  throw "unknow error, quit";
	}
      //seq[i].SetName(seqName);//set the sequence string name to short one
      
      //now need to figure out the gene index;

      //now go on further to parse the gene name and allele number
      p_arr[i].Set_Seq(seq.at(i));
      p_arr[i].Set_Gene(gene);
      gene_vec.push_back(gene);
    }
  
  //now we want to figure out gene_index
  //get unique gene vector first
  vector< string> unique_gene_vec;
  vector<unsigned> n_allele_vec;
  sort(gene_vec.begin(),gene_vec.end(), stringCompare_ext);
  
  //here, we need to sort the sequence string array, in order to 
  //later no compare to get the number of alleles and allele number 
  //for each gene
  sort(p_arr,p_arr+totalNumber,GenomicVCompare_bySequenceName);

  //now with sorted gene name vector
  //we need to get the unique genes
  unique_gene_vec.push_back(gene_vec.at(0));

  unsigned running_gene_index=0;
  unsigned running_allele=1; //this is the number of this current allele, it is starting from 1, not zero

  unsigned running_total_n_allele=1;
  
  p_arr[0].Set_GeneIndex(running_gene_index);
  p_arr[0].Set_Allele(running_allele);
  for(unsigned i=1;i<totalNumber;i++)
    {
      //now we need first get the unique genes
      //for each identifical one, we go next
      running_allele++;
      running_total_n_allele++;
      if(gene_vec.at(i).compare(gene_vec.at(i-1))!=0)
	{
	  //found a unique gene
	  running_gene_index++;
	  running_allele=1;
	  unique_gene_vec.push_back(gene_vec.at(i));
	  n_allele_vec.push_back(running_total_n_allele);
	  running_total_n_allele=1;
	}
      p_arr[i].Set_GeneIndex(running_gene_index);
      p_arr[0].Set_Allele(running_allele);
    }

  return totalNumber; 
}


//here we assume the name is arranged with a fix format
//for example >Jxxxx|IGHJ1*01|Homo Sapiens|F|......
string ParseSequenceName(const string& _seqName)
{
  //we split and return the name
  string buffer[500];
  int num=split_ext(_seqName,buffer, '|');
  if(num<2)
    {//something wrong 
      throw "bad format of the sequence name";
    }
  
  return buffer[1];//return second field 
}

//here we assume the name is arranged with a fix format
//for example > IGHD1*01
string ParseSequenceNameD(const string& _seqName)
{
  //we strip off the '> '
  string ret=_seqName;
  return ret.replace(0,2,"");
  //return buffer[1];//return second field 
}


string ParseGeneName(const string& _seqName)
{
  //we split and return the gene name
  string buffer[100];
  int num=split_ext(_seqName,buffer, '*');
  if(num<2)
    {//something wrong 
      throw "bad format of the gene sequence name";
    }
  return buffer[0];//return second field
}

bool GenomicJCompare_bySequenceName(const GenomicJ& ss1,const GenomicJ& ss2)
{
  int ret=ss1.Get_Name().compare(ss2.Get_Name());
  

  if(ret<=0)
    {
      return true;
    }
  else
    return false;
}


bool GenomicDCompare_bySequenceName(const GenomicD& ss1,const GenomicD& ss2)
{
  int ret=ss1.Get_Name().compare(ss2.Get_Name());
  

  if(ret<=0)
    {
      return true;
    }
  else
    return false;
}


//see the definition above
//the one to read or J
unsigned ReadGenomicJ(const string& _fastaFileName, GenomicJ** _gseg)
{
  //first read the files into vector
  vector<SequenceString> seq;
  unsigned totalNumber=ReadFasta(_fastaFileName, seq, true); 
  //get pointer to array of genomic segments
  *_gseg=new GenomicJ[totalNumber];

  GenomicJ* p_arr=*_gseg;
  vector<string> gene_vec;
  //cout<<"here it is"<<endl;
  //now start reading into the segments
  for(unsigned int i=0;i<totalNumber;i++)
    {
      string seqName;
      try
	{
	  seqName=ParseSequenceName(seq.at(i).GetName());
	}
      catch (const char* msg)
	{
	  throw msg;
	}
      catch (...)
	{
	  throw "unknow error, quit";
	}
      seq[i].SetName(seqName);//set the sequence string name to short one

      //now, we parse gene name
      string gene;
      try
	{
	  gene=ParseGeneName(seqName);
	}
      catch (const char* msg)
	{
	  throw msg;
	}
      catch (...)
	{
	  throw "unknow error, quit";
	}
      //seq[i].SetName(seqName);//set the sequence string name to short one
      
      //now need to figure out the gene index;

      //now go on further to parse the gene name and allele number
      p_arr[i].Set_Seq(seq.at(i));
      p_arr[i].Set_Gene(gene);
      gene_vec.push_back(gene);
    }
  
  //now we want to figure out gene_index
  //get unique gene vector first
  vector< string> unique_gene_vec;
  vector<unsigned> n_allele_vec;
  sort(gene_vec.begin(),gene_vec.end(), stringCompare_ext);
  
  //here, we need to sort the sequence string array, in order to 
  //later no compare to get the number of alleles and allele number 
  //for each gene
  sort(p_arr,p_arr+totalNumber,GenomicJCompare_bySequenceName);

  //now with sorted gene name vector
  //we need to get the unique genes
  unique_gene_vec.push_back(gene_vec.at(0));

  unsigned running_gene_index=0;
  unsigned running_allele=1; //this is the number of this current allele, it is starting from 1, not zero

  unsigned running_total_n_allele=1;
  
  p_arr[0].Set_GeneIndex(running_gene_index);
  p_arr[0].Set_Allele(running_allele);
  for(unsigned i=1;i<totalNumber;i++)
    {
      //now we need first get the unique genes
      //for each identifical one, we go next
      running_allele++;
      running_total_n_allele++;
      if(gene_vec.at(i).compare(gene_vec.at(i-1))!=0)
	{
	  //found a unique gene
	  running_gene_index++;
	  running_allele=1;
	  unique_gene_vec.push_back(gene_vec.at(i));
	  n_allele_vec.push_back(running_total_n_allele);
	  running_total_n_allele=1;
	}
      p_arr[i].Set_GeneIndex(running_gene_index);
      p_arr[0].Set_Allele(running_allele);
    }

  //now we have unique gene vector, we can go ahead to figure what is the total number of alleles for each gene


  return totalNumber;
}


unsigned ReadGenomicD(const string& _fastaFileName, GenomicD** _gseg)
{
 
  //first read the files into vector
  vector<SequenceString> seq;
  unsigned totalNumber=ReadFasta(_fastaFileName, seq, true); 
  //get pointer to array of genomic segments
  *_gseg=new GenomicD[totalNumber];

  GenomicD* p_arr=*_gseg;
  vector<string> gene_vec;
  //cout<<"here it is"<<endl;
  //now start reading into the segments
  for(unsigned int i=0;i<totalNumber;i++)
    {
      string seqName;
      try
	{
	  seqName=ParseSequenceNameD(seq.at(i).GetName());
	}
      catch (const char* msg)
	{
	  cout<<"seqaname:::"<<msg<<endl;
	  throw msg;
	}
      catch (...)
	{
	  throw "unknow error, quit";
	}
      seq[i].SetName(seqName);//set the sequence string name to short one

      //now, we parse gene name
      string gene;
      try
	{
	  gene=ParseGeneName(seqName);
	}
      catch (const char* msg)
	{
	  throw msg;
	}
      catch (...)
	{
	  throw "unknow error, quit";
	}
      //seq[i].SetName(seqName);//set the sequence string name to short one
      
      //now need to figure out the gene index;

      //now go on further to parse the gene name and allele number
      p_arr[i].Set_Seq(seq.at(i));
      p_arr[i].Set_Gene(gene);
      gene_vec.push_back(gene);
    }
  
  //now we want to figure out gene_index
  //get unique gene vector first
  vector< string> unique_gene_vec;
  vector<unsigned> n_allele_vec;
  sort(gene_vec.begin(),gene_vec.end(), stringCompare_ext);
  
  //now we need to have the array sorted in order to later on check
  //total number of allele and allele number
  sort(p_arr, p_arr+totalNumber, GenomicDCompare_bySequenceName);

  //now with sorted gene name vector
  //we need to get the unique genes
  unique_gene_vec.push_back(gene_vec.at(0));

  unsigned running_gene_index=0;
  unsigned running_allele=1; //this is the number of this current allele, it is starting from 1, not zero

  unsigned running_total_n_allele=1;
  
  p_arr[0].Set_GeneIndex(running_gene_index);
  p_arr[0].Set_Allele(running_allele);
  for(unsigned i=1;i<totalNumber;i++)
    {
      //now we need first get the unique genes
      //for each identifical one, we go next
      running_allele++;
      running_total_n_allele++;
      if(gene_vec.at(i).compare(gene_vec.at(i-1))!=0)
	{
	  //found a unique gene
	  running_gene_index++;
	  running_allele=1;
	  unique_gene_vec.push_back(gene_vec.at(i));
	  n_allele_vec.push_back(running_total_n_allele);
	  running_total_n_allele=1;
	}
      p_arr[i].Set_GeneIndex(running_gene_index);
      p_arr[0].Set_Allele(running_allele);
    }

  //now we have unique gene vector, we can go ahead to figure what is the total number of alleles for each gene


  return totalNumber; 
}


void FlipSequenceString(SequenceString& _ss)
{
  string temp=_ss.GetSequence();
  temp=flipStr(temp);
    _ss.SetSequence(temp);
}

vector<string> DetermineOutputFileNames(const string& _outFileNameBase, const unsigned& _NPerFile, const unsigned& _totalNSeq)
{
  vector<string> vecFileNames;

  unsigned NFiles;
  if(_totalNSeq/_NPerFile*_NPerFile==_totalNSeq)
    {
      NFiles=_totalNSeq/_NPerFile;
    }
  else
    {
      NFiles=_totalNSeq/_NPerFile+1;
    }
  
  for(unsigned i=0;i<NFiles;i++)
    {
      ostringstream convert;
      convert<<i;
      string temp(_outFileNameBase);
      temp.append("_");
      temp.append(convert.str());
      temp.append(".align");
      vecFileNames.push_back(temp);
    }
  return vecFileNames;
}



