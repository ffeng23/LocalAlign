#include "DataIO.hpp"

#include <iostream>
#include <fstream>
static string  version("0.1"); //for serialization

//static method to write the version information
static void WriteVersionInfo(const string& _version, ofstream& _ofs)
{ 
  if(!_ofs.is_open())
    {
      cout<<"******ERROR: can not open file, quit..."<<endl;
      cerr<<"*****error: can not open file, quit..."<<endl;
      //return false;
      exit(-1);
    }
//start serialize the version number
  unsigned len=_version.length()+1;
  const char* p_char=(char*)&len;
  _ofs.write(p_char, sizeof(unsigned));
  p_char=version.c_str();
  _ofs.write(p_char, version.length()+1);
} 
static void DoSerilizationStep(const Alignment_Object& _v_align, 
			       const Alignment_D& _d_align, 
			       const Alignment_Object& _j_align,
			       const SequenceString& _seq, const unsigned& _step,
			       ofstream& _ofs)
{
  if(!_ofs.is_open())
    {
      cout<<"******ERROR: can not open file, quit..."<<endl;
      cerr<<"******ERROR: can not open file, quit..."<<endl;
      exit(-1);
      //return false;
    }
  const char* p_char;
  //serialize sequence
  _seq.Serialize(_ofs);
  //serialize V
  _v_align.Serialize(_ofs);

  //serialize D
  _d_align.Serialize(_ofs);
  //serialize J
  _j_align.Serialize(_ofs);
  //serialize number of written so far
  p_char=(char*)&_step;
  _ofs.write(p_char, sizeof(unsigned));
}
//=========>serialization of alignment datas
bool DoSerialization
(bool const* const* _vdj_align_ok_arrays, Alignment_Object const* const* _v_align_arrays,
 Alignment_D const* const* _d_align_arrays, Alignment_Object const* const* _j_align_arrays,
 vector<SequenceString>::const_iterator* _it_arrays , 
 const unsigned* _numOfAlignedSequenceForThread, 
 const unsigned* _numOfInputSequencesForThread, const unsigned& _numOfThread, 
 const string& _fileName
 )
//NOTE about the structure of the serialization
//1)first we store the version "ver 0.01"
//2) the total number of aligment sets, an unsigned, 
//3) each set contains seq, v_align, 
//           d align and v align.
//4) an unsigned used to indicating the number of set so far.
//       this number could be used to indicating the integraty of the serialization file
//       when doing deserialization.
//
{
  //open up a file first
  ofstream ofs(_fileName.c_str(), std::ios::binary);
  if(!ofs.is_open())
    {
      cout<<"******ERROR: can not open file, quit..."<<endl;
      return false;
    }
  const char* p_char;//pointer used to direct tht writing to the file stream
  //start serialize the version number
  /*unsigned len=version.length()+1;
  p_char=(char*)&len;
  ofs.write(p_char, sizeof(unsigned));
  p_char=version.c_str();
  ofs.write(p_char, version.length()+1);
  */
  WriteVersionInfo(version, ofs);
  //start serialize the total number of alignment sets
  //start getting the totalNumber of align ok
  unsigned alignedTotal1=0, alignedTotal2=0;
  for(unsigned i=0;i<_numOfThread;i++)
    {
      alignedTotal1+=_numOfAlignedSequenceForThread[i];
      for(unsigned j=0;j<_numOfInputSequencesForThread[i];j++)
	{
	  if(_vdj_align_ok_arrays[i][j])
	    alignedTotal2++;
	}
    }
  //now for safety purpose, we need to check for equality
  if(alignedTotal1!=alignedTotal2)
    {
      cout<<"when processing the aligned data, the number of aligned is inconsistent, quit"<<endl;
      return false;
    }
  p_char=(char*)&alignedTotal2;
  ofs.write(p_char, sizeof(unsigned));

  //looping to write the alignment set
  unsigned runningCount=0;
  for(unsigned i=0;i<_numOfThread;i++)
    {
      //advance(_it_array[
      for(unsigned j=0;j<_numOfInputSequencesForThread[i];j++)
	{
	  if(_vdj_align_ok_arrays[i][j])
	    {
	      runningCount++;
	      /*
	      //serialize sequence
	      (*_it_arrays[i]).Serialize(ofs);
	      //serialize V
	      _v_align_arrays[i][j].Serialize(ofs);

	      //serialize D
	      _d_align_arrays[i][j].Serialize(ofs);
	      //serialize J
	      _j_align_arrays[i][j].Serialize(ofs);
	      //serialize number of written so far
	      p_char=(char*)&runningCount;
	      ofs.write(p_char, sizeof(unsigned));
	      */
	      DoSerilizationStep(_v_align_arrays[i][j], _d_align_arrays[i][j],
				 _j_align_arrays[i][j], *_it_arrays[i],
				 runningCount, ofs);
	    }//a ok align, do writing, end of if block
	  _it_arrays[i]++;
	}//end of each alignment of each thread, for loop
      
    }//end of thread for loop
  
  
  //close the file
  ofs.close();
  return true;

}



//this method is used to only select first "_numOfAlignments" of aligns from
//the original file.
bool DoSerialization
(const Alignment_Object * _v_align_arrays,
 const Alignment_D * _d_align_arrays, const Alignment_Object * _j_align_arrays,
 const SequenceString* _seq , const unsigned& _totalNumOfAligns,
 const string& _fileName,
 const unsigned& _numOfAlignments
 )
{
  //open up a file first
  ofstream ofs(_fileName.c_str(), std::ios::binary);
  if(!ofs.is_open())
    {
      cout<<"******ERROR: can not open file, quit..."<<endl;
      return false;
    }
  const char* p_char;//pointer used to direct tht writing to the file stream
  
  WriteVersionInfo(version, ofs);

  //start doing the serilization
  //check to make sure we are deal with the correct numbers
  unsigned num=_numOfAlignments;
  if(_numOfAlignments>_totalNumOfAligns)
    {
      num=_totalNumOfAligns;
    }

  p_char=(char*)&num;
  ofs.write(p_char, sizeof(unsigned));  

  
  //for loop
  for(unsigned i=0;i<num;i++)
    {
      DoSerilizationStep(_v_align_arrays[i], _d_align_arrays[i], _j_align_arrays[i],
			 _seq[i], i+1, ofs);
    }
  
  ofs.close();
  return true;
}

//=========>deserialization of alignment datas
//reading the alignment from the disk, 
//the caller DON'T initialize the output arrays, but only declare it.
//also the caller need to delete/clean up the memory.
//
//Note: for the outside caller, please please make sure that 
//the passed in input variable (*_v_align, *_d_align and *_j_align)
//    is null; that means those variable don't point to any memory. 
//    or most importantly there is no memory leak. The caller has
//    deallocated the memory those variable pointed to previously
//******************************************************
bool DoDeserialization
(const string& _fileName, unsigned& _numOfAligned,
 Alignment_Object** _v_align, Alignment_D** _d_align, Alignment_Object** _j_align,
 SequenceString** _seq
 )
{
  //open up a file first
  ifstream ifs(_fileName.c_str(), std::ios::binary);
  if(!ifs.is_open())
    {
      cout<<"******ERROR: can not open file, quit..."<<endl;
      return false;
    }
  //need to follow the structure of the serialization
  //first need to check for the version information, if the version information
  //do't match have to issue warning........
  unsigned len;
  char* p_char=(char*)&len;
  ifs.read(p_char, sizeof(unsigned));
  p_char=new char[len];
  ifs.read(p_char, len);
  string v_read(p_char);
  delete[] p_char;
  //check the version information
  if(v_read.compare(version)!=0)
    {
      cout<<"version is not right, stop reading........"<<endl;
      ifs.close();
      return false;
    }
  //read the total number of alignments
  p_char=(char*)&_numOfAligned;
  ifs.read(p_char, sizeof(unsigned));

  if((signed)_numOfAligned<0)
    {
      //someting wrong
      cout<<"impossible number of alignment object store, quit ("<<_numOfAligned<<")"<<endl;
      return false;
    }
  if(*_seq!=NULL)
    {
      delete [] *_seq;
    }
  if(*_v_align!=NULL)
    {
      delete[] *_v_align;
    }
  if(*_d_align!=NULL)
    {
      delete[] *_d_align;
    }
  if(*_j_align!=NULL)
    {
      delete[] *_j_align;
    }
  
  //initialize memory
  *_seq=new SequenceString[_numOfAligned];
  *_v_align=new Alignment_Object[_numOfAligned];
  *_d_align=new Alignment_D[_numOfAligned];
  *_j_align=new Alignment_Object[_numOfAligned];
  unsigned count;
  //now start reading
  for(unsigned i=0;i<_numOfAligned;i++)
    {
      //read sequence
      (*_seq)[i].Deserialize(ifs);
      //read v
      (*_v_align)[i].Deserialize(ifs);
      //read d
      (*_d_align)[i].Deserialize(ifs);
      //read j
      (*_j_align)[i].Deserialize(ifs);
      //read count
      p_char=(char*)&count;
      ifs.read(p_char, sizeof(unsigned));
      if(count!=i+1)
	{
	  cout<<"the file is not kept intact ,something wrong, please check its integrity"<<endl;
	  return false;
	}
    }//end of for loop

    //clean up
  ifs.close();
  return true;
}

//read in the genomic segments from the files into variables.
bool DoGenomicTemplateReading
(const string& _genomicVFileName,
 const string& _genomicDFileName,
 const string& _genomicJFileName,
 /*output*/
 GenomicV** _genV, GenomicD** _genD, GenomicJ** _genJ,
 unsigned& _totalNumV, unsigned& _totalNumD, unsigned& _totalNumJ
 )
{
  _totalNumJ=  ReadGenomicJ(_genomicJFileName, _genJ);
  
  cout<<_totalNumJ<<" J genomic segments are read in."<<endl; 
  /*cout<<"\tshowing the J seq 1:"<<genJ[1].Get_Sequence()<<endl;
  cout<<totalNumJ<<" J genomic segments are read in."<<endl; 
  for(unsigned i=0;i<totalNumJ;i++)
    {
      cout<<i<<":"<<genJ[i].Get_Seq().toString()<<endl;
      cout<<"\t==>geneIndex:"<<genJ[i].Get_GeneIndex()<<endl;
      cout<<"\t==>n_allele:"<<genJ[i].Get_n_alleles()<<endl;
      cout<<"\t==>allele:"<<genJ[i].Get_Allele()<<endl;
    }
  */
  _totalNumD=  ReadGenomicD(_genomicDFileName, _genD);
  /*
  cout<<_totalNumD<<" D genomic segments are read in."<<endl; 
  for(unsigned i=0;i<_totalNumD;i++)
    {
      cout<<i<<":"<<(*_genD)[i].Get_Seq().toString()<<endl;
      cout<<"\t==>geneIndex:"<<(*_genD)[i].Get_GeneIndex()<<endl;
      cout<<"\t==>n_allele:"<<(*_genD)[i].Get_n_alleles()<<endl;
      cout<<"\t==>allele:"<<(*_genD[i]).Get_Allele()<<endl;
    }
  */
  _totalNumV=  ReadGenomicV(_genomicVFileName ,_genV);
  
  cout<<_totalNumV<<" V genomic segments are read in."<<endl; 
  /*for(unsigned i=0;i<totalNumV;i++)
    {
      cout<<i<<":"<<genV[i].Get_Seq().toString()<<endl;
      cout<<"\t==>geneIndex:"<<genV[i].Get_GeneIndex()<<endl;
      cout<<"\t==>n_allele:"<<genV[i].Get_n_alleles()<<endl;
      cout<<"\t==>allele:"<<genV[i].Get_Allele()<<endl;
      }*/

  /*  cout<<100<<":"<<genV[100].Get_Seq().toString()<<endl;
  cout<<216<<":"<<genV[216].Get_Seq().toString()<<endl;
  cout<<217<<":"<<genV[217].Get_Seq().toString()<<endl;
  */
  return true;
}



