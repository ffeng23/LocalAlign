//#include "SequenceHandler.hpp"
#include "Accessory/FastaHandler.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>

#include "SequenceHandlerCheckUPM_GCGvsATG.hpp"
#include "OverlapAlignment.hpp"
#include "AlignmentString.hpp"
#include "LocalAlignment.hpp"
#include "SequenceHandlerCommon.hpp"


void CheckingUPMGCGvsATC( vector<SequenceString>& _vecReverse, vector<SequenceString>& _vecSeq, 
		     ScoreMatrix* _sm, const double& _gapOpen, const double& _gapExtension,
			 const double& _mismatchRateThreshold, const unsigned _minimumOverlapLength, 
		      const string& _mapReverse_fname, const string& _mapNone_fname )
{
  unsigned int numOfSeqsUnit=20000;
  unsigned int timeOfWriting=1;
  
  //now we need to prepare the output vector, now we will prepare the outfile for each entry in _vecForward
  //they are the alleles for different isotypes (IgG/M), different subtypes (IgG1/2/3/4) and alleles IgG1*01/*02/*03, etc
  //
  
  vector<SequenceString> vec_mapReverse;
  vector<SequenceString> vec_mapReverse_nonATG;
  vector<SequenceString> vec_mapNone;
  unsigned int stat_ATG=0;
  unsigned int stat_GCG=0;
  unsigned int stat_reverseMapped=0;
  
  //now, we need to prepare the output files
  AlignmentString* tempAS_arr;
  unsigned int numOfLocalAlignments=2;
  string strPattern;
  string strSubject;
  bool foundATGFlag=false;
  
  double mismatch_rate=0.0;
  cout<<"\n";
  for(unsigned int i=0;i<_vecSeq.size();i++)
    {
      foundATGFlag=false;
      //printing progress
      if(i/numOfSeqsUnit *numOfSeqsUnit==i)
	{
	  cout<<"..."<<i<<"/"<<_vecSeq.size();
	  flush(cout);
	  //need to write to file.
	  
	}
      /*double bestForwardScore= -10000000;
	AlignmentString	bestForwardAlign;*/
      double bestReverseScore= -1000000;
      AlignmentString bestReverseAlign;
      
      //unsigned int 	bestForwardIndex=0;
      unsigned int	bestReverseIndex=0;
	
      //bool 	foundForwardFlag=false;
      bool	foundReverseFlag=false;
      //bool      foundCrossOverFlag=false;
      //bool foundBreakOutsideConstantFlag=false;

      //cout<<"map reverse set"<<endl;
      //reverse side should be mapped on the end of the reads, need to reverse complement the sequence too
      for(unsigned int k=0;k< _vecReverse.size();k++)
	{
	  //cout<<"start doing k:"<<k<<endl;
	  SequenceString reverseComplementReverse=ReverseComplement(_vecReverse.at(k));
	  //cout<<"after revcomp"<<endl;
	  //OverlapAlignment ola ( &(_vecSeq.at(i)), &reverseComplementReverse, _sm, _gapOpen, _gapExtension,1);
	  LocalAlignment ola ( &(_vecSeq.at(i)), &reverseComplementReverse, _sm, _gapOpen, _gapExtension,1,numOfLocalAlignments);
	  //cout<<"after alignment:"<<k<<endl;
	  //need to get the mismatch rate
	  tempAS_arr=ola.GetAlignmentArr();
	  //need to get the mismatch rate
	  strPattern=tempAS_arr[0].GetPattern(true);
	  strSubject=tempAS_arr[0].GetSubject(true);
	  mismatch_rate=1-CompareStrings(strPattern, strSubject)/((double)strPattern.length());
	  	
	  //#check the best score
	  if(tempAS_arr[0].GetScore()>bestReverseScore&&mismatch_rate>_mismatchRateThreshold&&strPattern.length()>_minimumOverlapLength)
	    {
	      //##we need to see the offset too, only important to the pattern, here the 
	      //	##the pattern is the long one, we align the primer sequence against
	      //	###the primer should on the both ends, not too far,
	      //if( 
	      // (_vecSeq.at(i).GetLength()-tempAS.GetPatternIndexEnd())<_offsetReverse )
	      //{//we good
		  bestReverseScore=tempAS_arr[0].GetScore();
		  bestReverseAlign=tempAS_arr[0];//with name
		  bestReverseIndex=k;
		  foundReverseFlag=true;
		  //}
	    }
	}

      //cout<<"Done with mapping:i="<<i<<endl;
      
      //doing stats, check for the GCG vs ATG
      
      //cout<<"bestReverseAlign.GetSubjectIndexStart():"<<bestReverseAlign.GetSubjectIndexStart()<<endl;
      if(foundReverseFlag)
	{
	  stat_reverseMapped++;
	  if(bestReverseAlign.GetSubjectIndexStart()<=3) //it has to start in the beginning
	    {
	      //check for GCG or ATG
	      strPattern=bestReverseAlign.GetPattern();
	      if(strPattern.find("ATG")<=3)
		{
		  foundATGFlag=true;
		  stat_ATG++;
		}
	      else
		{
		  if(strPattern.find("GCG")<=3)
		    stat_GCG++;
		}
	    }
	}
      
      
      //***************************output*****************
      //prepare the aligned for debugging purpose
      //#the forward 
      //string leadingSpaceOriginal("");
      //string leadingSpaceForward("");
      string leadingSpaceReverse("");
      string replaceOne;

      SequenceString tempLstR;
      SequenceString tempLstSeq;

      //#the revverse
      //cout<<"reverse set output read"<<endl;
	if(foundReverseFlag)
	  {
	    //#now we need to add the aligned sequence to replace the original one
	    //#we assume the original one is longer than the aligned one, it has to be
	    //	 #this is the intial part
	    SequenceString rcReverseSeq=ReverseComplement(_vecReverse.at(bestReverseIndex));
	    replaceOne=rcReverseSeq.GetSequence().substr(0, bestReverseAlign.GetSubjectIndexStart());
	    //cout<<"\t\t*****subject start:"<<bestReverseAlign.GetSubjectIndexStart()<<";subject end:"<<bestReverseAlign.GetSubjectIndexEnd()<<endl;
	    //cout<<"\t\t***replaceOne:"<<replaceOne<<endl;
	    //#aligned part
	    replaceOne.append( bestReverseAlign.GetSubject(true));
	    //cout<<"\t\t***replaceOne22:"<<replaceOne<<endl;
	    //#last part unaligned
	    replaceOne.append( rcReverseSeq.GetSequence().substr(bestReverseAlign.GetSubjectIndexEnd()+1, rcReverseSeq.GetLength()));
	    //cout<<"\t\t***replaceOne333:"<<replaceOne<<endl;
	    tempLstR.SetName(_vecReverse.at(bestReverseIndex).GetName());
	    //tempLstR.SetSequence(replaceOne);
	    //#now we need to figure out how the leading space to put in front of reverse one
	    if(bestReverseAlign.GetPatternIndexStart() >= bestReverseAlign.GetSubjectIndexStart())
	      {
		unsigned int templen=bestReverseAlign.GetPatternIndexStart()- bestReverseAlign.GetSubjectIndexStart();
		for(unsigned int p=0;p<templen;p++)
		  {
		    leadingSpaceReverse.push_back('-');
		  }
		tempLstR.SetSequence(leadingSpaceReverse+replaceOne);
	      }
	    else
	      {//#here, in this case, the adaptor+primer is longer than the seqs just by alignment, then we need to simply remove some leading part of the adaptor primer
		tempLstR.SetSequence(replaceOne.substr( bestReverseAlign.GetSubjectIndexStart()-bestReverseAlign.GetPatternIndexStart(), replaceOne.length()));
	      }
	  }
	else
	  {
	    tempLstR.SetSequence("");
	    tempLstR.SetName("NoMatch");
	  }
	//cout<<"end of reverse on, tempLstR.seq:"<<tempLstR.toString()<<endl;
	//#now we need to take care of the read sequence alignment string
	//#on the reverse part first
	//cout<<"seq string output read...."<<endl;
	//unsigned int spaceCarryOverFTR=0;//####this is the leading space for reverse adaptor primer, because the insertion in the forward alignment
	tempLstSeq.SetSequence(_vecSeq.at(i).GetSequence());
	tempLstSeq.SetName(_vecSeq.at(i).GetName());
	if(foundReverseFlag)
	{
	  //cout<<"1.."<<endl;
	  replaceOne=_vecSeq.at(i).GetSequence().substr(0, bestReverseAlign.GetPatternIndexStart());
	  //#aligned part
	  //cout<<"2.."<<endl;
	  replaceOne.append( bestReverseAlign.GetPattern(true));
	  //#last part unaligned

	  //cout<<"3..."<<endl;
	  replaceOne.append( _vecSeq.at(i).GetSequence().substr(bestReverseAlign.GetPatternIndexEnd()+1));//nchar(ff[[i]])), sep="");
	  tempLstSeq.SetSequence(replaceOne);
	  //tempLstSeq.SetName(_vecSeq.at(i).GetName());
	}
	
	
	//******************write to vectors***************
	//#now put the sequences to the correct vectors
	//cout<<"\tready to output strings........"<<endl;
	vector<SequenceString>* p_vec_map;
	if(foundReverseFlag)
	  {
	    if(foundATGFlag)
	      p_vec_map=&vec_mapReverse;
	    else
	      p_vec_map=&vec_mapReverse_nonATG;
	  }
	else
	  {
	    
		p_vec_map=&vec_mapNone;
		
	  }
	
	//now we need to store the data into the correct vectors/files
	p_vec_map->push_back(tempLstSeq); 
	//p_vec_map->push_back(tempLstF);
	p_vec_map->push_back(tempLstR);
	
	//cout<<"start writing output........"<<endl;
	if(i>=timeOfWriting*numOfSeqsUnit) //#write once every 20000 sequences
	  {
	    timeOfWriting++;
	    string t_fileName;
	  //cout<<"i round:"<<i<<endl;
	  
	  //reverse
	      if(vec_mapReverse.size()>0)
		{
		  //cout<<"\t------writing files at i-----:"<<s<<endl;
		  //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
		  		  
		  t_fileName=_mapReverse_fname;
		  WriteFasta(t_fileName, vec_mapReverse,100, ofstream::app);
		  //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
		  vec_mapReverse.clear();
		}
	    
	  //none
	      if(vec_mapNone.size()>0)
		{
		  //cout<<"\t------writing files at i-----:"<<s<<endl;
		  //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
		  		  
		  t_fileName=_mapNone_fname;
		  WriteFasta(t_fileName, vec_mapNone,100, ofstream::app);
		  //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
		  vec_mapNone.clear();
		}
	      //nonATG
	      if(vec_mapReverse_nonATG.size()>0)
		{
		  //cout<<"\t------writing files at i-----:"<<s<<endl;
		  //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
		  		  
		  t_fileName=_mapReverse_fname+"_nonATG";
		  WriteFasta(t_fileName, vec_mapReverse_nonATG,100, ofstream::app);
		  //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
		  vec_mapReverse_nonATG.clear();
		}
	  
	  //cout<<"done map both"<<endl;
	}//end of each 1000 seqs read write
      
    }//end of for loop of sequence data vec

  //one last writing
  cout<<"start writing last round of output........"<<endl;
  //reverse	  
  string t_fileName;
  if(vec_mapReverse.size()>0)
    {
      //cout<<"\t------writing files at i-----:"<<s<<endl;
      //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
      
      t_fileName=_mapReverse_fname;
      WriteFasta(t_fileName, vec_mapReverse,100, ofstream::app);
      //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
      vec_mapReverse.clear();
    }
  
  //none
  if(vec_mapNone.size()>0)
    {
      //cout<<"\t------writing files at i-----:"<<s<<endl;
      //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
		  		  
      t_fileName=_mapNone_fname;
      WriteFasta(t_fileName, vec_mapNone,100, ofstream::app);
      //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
      vec_mapNone.clear();
    }

  if(vec_mapReverse_nonATG.size()>0)
    {
      //cout<<"\t------writing files at i-----:"<<s<<endl;
      //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
		  		  
      t_fileName=_mapReverse_fname+"_nonATG";
      WriteFasta(t_fileName, vec_mapReverse_nonATG,100, ofstream::app);
      //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
      vec_mapReverse_nonATG.clear();
    }
  
  cout<<"done with one last write"<<endl;
  
  cout<<"---------Summary---------"<<endl;
  cout<<"\ttotal sequences: "<<_vecSeq.size()<<"\n"
      <<"\tmap reverse sequences:"<<stat_reverseMapped <<"\n"
      <<"\t\"ATG\" sequences:"<<stat_ATG<<"\n"
      <<"\t\"GCG\" sequences:"<<stat_GCG<<"\n";
    
  cout<<endl;
  //cleanup
  
}

