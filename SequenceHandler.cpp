#include "SequenceHandler.hpp"
#include "FastaHandler.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "OverlapAlignment.hpp"
#include "AlignmentString.hpp"


void ProcessAdaptorSequences(const string& _adaptorFile, const string& _barcodeFile, const string& _forwardFile, 
			     const string& _reverseFile, vector<SequenceString>& _vecForward, 
			     vector<SequenceString>& _vecReverse )
{
  //first, we need to read the sequeces
  vector<SequenceString> adaptor;
  vector<SequenceString> barcode;
  vector<SequenceString> forward;
  vector<SequenceString> reverse;
  
  //
  ReadFasta(_adaptorFile, adaptor);
  ReadFasta(_barcodeFile, barcode);
  ReadFasta(_forwardFile, forward);
  ReadFasta(_reverseFile, reverse);

  SequenceString adaptorA(adaptor.at(0).GetName(), adaptor.at(0).GetSequence());
  SequenceString adaptorB(adaptor.at(1).GetName(), adaptor.at(1).GetSequence());

  //go combine them together
  //forward set -- with adaptorA
  string tempStr;
  string tempName;
  for(unsigned int i=0;i<barcode.size();i++)
    {
      for(unsigned int j=0;j<forward.size();j++)
	{
	  tempStr=adaptorA.GetSequence();
	  tempStr+=barcode.at(i).GetSequence();
	  tempStr+=forward.at(j).GetSequence();
	  tempName=adaptorA.GetName();
	  tempName+=barcode.at(i).GetName();
	  tempName+=forward.at(j).GetName();
	  //
	  SequenceString tempSS(tempName, tempStr);
	  _vecForward.push_back(tempSS);
	}
    }

  //reverse set
  //with adaptor B
  for(unsigned int i=0;i<barcode.size();i++)
    {
      for(unsigned int j=0;j<reverse.size();j++)
	{
	  tempStr=adaptorB.GetSequence();
	  tempStr+=barcode.at(i).GetSequence();
	  tempStr+=reverse.at(j).GetSequence();
	  tempName=adaptorB.GetName();
	  tempName+=barcode.at(i).GetName();
	  tempName+=reverse.at(j).GetName();
	  //
	  SequenceString tempSS(tempName, tempStr);
	  _vecReverse.push_back(tempSS);
	}
    }
  //we are doing here
}


SequenceString ReverseComplement(SequenceString& seq)
{
  string tempStr=seq.GetSequence();
  SequenceString temp(seq.GetName(), "");
  string tempStrReturn("");
  //cout<<"temStr (seq get sequence):"<<tempStr<<endl;
  //cout<<"length:"<<tempStr.length()<<endl;
  for(unsigned int i=tempStr.length();i>0;i--)
    {
      //cout<<"\tloop i="<<i<<endl;
      switch(tempStr.at(i-1))
	{
	case 'A':
	case 'a':
	  tempStrReturn.push_back('T');
	  break;
	case 'T':
	case 't':
	  tempStrReturn.push_back('A');
	  break;
	case 'U':
	case 'u':
	  tempStrReturn.push_back('A');
	  break;
	case 'G':
	case 'g':
	  tempStrReturn.push_back('C');
	  break;
	case 'C':
	case 'c':
	  tempStrReturn.push_back('G');
	  break;
	case 'Y':
	case 'y':
	  tempStrReturn.push_back('R');
	  break;
	case 'R':
	case 'r':
	  tempStrReturn.push_back('Y');
	  break;
	case 'S':
	case 's':
	  tempStrReturn.push_back('S');
	  break;
	case 'W':
	case 'w':
	  tempStrReturn.push_back('W');
	  break;
	case 'K':
	case 'k':
	  tempStrReturn.push_back('M');
	  break;
	case 'M':
	case 'm':
	  tempStrReturn.push_back('K');
	  break;
	case 'B':
	case 'b':
	  tempStrReturn.push_back('V');
	  break;
	case 'D':
	case 'd':
	  tempStrReturn.push_back('H');
	  break;
	case 'H':
	case 'h':
	  tempStrReturn.push_back('D');
	  break;
	case 'V':
	case 'v':
	  tempStrReturn.push_back('B');
	  break;
	case 'N':
	case 'n':
	default:
	  tempStrReturn.push_back('N');
	  break;
	}
      
    }

  temp.SetSequence(tempStrReturn);

  return temp;
}


void MappingAdaptors(vector<SequenceString>& _vecForward, vector<SequenceString>& _vecReverse, vector<SequenceString>& _vecSeq, 
		     ScoreMatrix* _sm, const double& _gapOpen, const double& _gapExtension, const unsigned int& _trim,
		     const double& _mismatchRateThreshold, const unsigned _minimumOverlapLength, const unsigned int& _offsetForward, const unsigned int& _offsetReverse, 
		     const string& _mapBoth_fname, const string& _mapForward_fname,
		     const string& _mapReverse_fname, const string& _mapNone_fname)
{
  int numOfSeqsUnit=5000;
  ofstream ofs_both(_mapBoth_fname.c_str());
  ofstream ofs_forward(_mapForward_fname.c_str());
  ofstream ofs_reverse(_mapReverse_fname.c_str());
  ofstream ofs_none(_mapNone_fname.c_str());

  if(!ofs_both.is_open()||!ofs_forward.is_open()||!ofs_reverse.is_open()||!ofs_none.is_open())
    {
      cout<<">>>>>>ERROR:the output file \"*map*.fasta\" can not be opened, quit....\n";
      exit(-1);
    }
  
  vector<SequenceString> vec_mapBoth;
  vector<SequenceString> vec_mapForward;
  vector<SequenceString> vec_mapReverse;
  vector<SequenceString> vec_mapNone;

  vector<int> vec_both_len;
      vector<int> vec_forward_len;
      vector<int> vec_reverse_len;
      vector<int> vec_none_len;

  AlignmentString tempAS;
      string strPattern;
      string strSubject;
      
      double mismatch_rate=0.0;
      
  for(unsigned int i=0;i<_vecSeq.size();i++)
    {
      //printing progress
      if(i/numOfSeqsUnit *numOfSeqsUnit==i)
	{
	  cout<<"..."<<i<<"/"<<_vecSeq.size();
	  flush(cout);
	  //need to write to file.
	  
	}
      double bestForwardScore= -10000000;
      AlignmentString	bestForwardAlign;
      double bestReverseScore= -1000000;
      AlignmentString bestReverseAlign;
	
      unsigned int 	bestForwardIndex=0;
      unsigned int	bestReverseIndex=0;
	
      bool 	foundForwardFlag=false;
      bool	foundReverseFlag=false;
      //cout<<"Map forward set"<<endl;
      //forward has to be mapped to the beginning!!!
      for(unsigned int j=0;j<_vecForward.size();j++)
	{
	  //cout<<"forward set:"<<j<<endl;
	  //cout<<_vecForward.at(j).toString()<<endl;
	  OverlapAlignment ola (&(_vecSeq.at(i)), &(_vecForward.at(j)), _sm, _gapOpen, _gapExtension,1);
	  tempAS=ola.GetAlignment();
	  //cout<<"\talignment score:"<<tempAS.GetScore()<<endl;
	  //need to get the mismatch rate
	  strPattern=tempAS.GetPattern(true);
	  strSubject=tempAS.GetSubject(true);
	  //cout<<"\tstrPattern:"<<strPattern<<endl;
	  //cout<<"\tstrSubject:"<<strSubject<<endl;

	  mismatch_rate=1-CompareStrings(strPattern, strSubject)/((double)strPattern.length());
	  //cout<<"\tcompare:"<<CompareStrings(strPattern, strSubject)<<";length():"<<strPattern.length()<<endl;
	  //cout<<"\tmismatch_rate:"<<mismatch_rate<<endl;
		
	  //check the best score and mismatch rate
	  if(tempAS.GetScore()>bestForwardScore&&mismatch_rate>_mismatchRateThreshold&&strPattern.length()>_minimumOverlapLength)
		{
		  //we need to see the offset too, only important to the pattern, here the 
		  //the pattern is the long one, we align the primer sequence against
		  //the primer should be in the beginning, not too far,
		  if(tempAS.GetPatternIndexStart()<_offsetForward )
		    {//we good
		      //cout<<"\t***get one bigger"<<endl;
		      bestForwardScore=tempAS.GetScore();
		      bestForwardAlign=tempAS;//with name
		      bestForwardIndex=j;
		      foundForwardFlag=true;
		    }
		}
	}
      //cout<<"map reverse set"<<endl;
      //reverse side should be mapped on the end of the reads, need to reverse complement the sequence too
      for(unsigned int k=0;k< _vecReverse.size();k++)
	{
	  //cout<<"start doing k:"<<k<<endl;
	  SequenceString reverseComplementReverse=ReverseComplement(_vecReverse.at(k));
	  //cout<<"after revcomp"<<endl;
	  OverlapAlignment ola ( &(_vecSeq.at(i)), &reverseComplementReverse, _sm, _gapOpen, _gapExtension,1);
	  //cout<<"after alignment:"<<k<<endl;
	  //need to get the mismatch rate
	  tempAS=ola.GetAlignment();
	  //need to get the mismatch rate
	  strPattern=tempAS.GetPattern(true);
	  strSubject=tempAS.GetSubject(true);
	  mismatch_rate=1-CompareStrings(strPattern, strSubject)/((double)strPattern.length());
	  	
	  //#check the best score
	  if(tempAS.GetScore()>bestReverseScore&&mismatch_rate>_mismatchRateThreshold&&strPattern.length()>_minimumOverlapLength)
	    {
	      //##we need to see the offset too, only important to the pattern, here the 
	      //	##the pattern is the long one, we align the primer sequence against
	      //	###the primer should on the both ends, not too far,
	      if( 
		 (_vecSeq.at(i).GetLength()-tempAS.GetPatternIndexEnd())<_offsetReverse )
		{//we good
		  bestReverseScore=tempAS.GetScore();
		  bestReverseAlign=tempAS;//with name
		  bestReverseIndex=k;
		  foundReverseFlag=true;
		}
	    }
	}
	
      //cout<<"Done with mapping"<<endl;
      //###done with mapping, now we need to make the output ready
      
      //#the forward 
      string leadingSpaceOriginal("");
      string leadingSpaceForward("");
      string leadingSpaceReverse("");
      string replaceOne;
      SequenceString tempLstF;
      SequenceString tempLstR;
      SequenceString tempLstSeq;
      
      //cout<<"forwardSet Output read"<<endl;
      if(foundForwardFlag)
	{
	  //###we need to figure out the leading spaces in front of the original sequences
	  unsigned int startOriginal=bestForwardAlign.GetPatternIndexStart();
	  unsigned int startSubject=bestForwardAlign.GetSubjectIndexStart();
	  if(startSubject>startOriginal)
	    {
	      //leadingSpaceOriginalsg=as.character();
	      for(unsigned int i=0;i<startSubject-startOriginal;i++)
		{
		  leadingSpaceOriginal.push_back('-');
		}
	    }
	  else
	    {
	      for(unsigned int i=0;i<startOriginal-startSubject;i++)
		{
		  leadingSpaceForward.push_back('-');
		}
	      
	      //leadingSpaceForward<-paste(rep("-",startOriginal-startSubject), collapse = '')
	    }
	  
	  //now we need to add the aligned sequence to replace the original one
	  //#we assume the original one is longer than the aligned one, it has to be
	  // #this is the intial part
	  replaceOne=_vecForward.at(bestForwardIndex).GetSequence().substr(0, bestForwardAlign.GetSubjectIndexStart());
	  //#aligned part
	  replaceOne.append( bestForwardAlign.GetSubject(true));
	  //#last part unaligned
	  replaceOne.append(_vecForward.at(bestForwardIndex).GetSequence().substr(bestForwardAlign.GetSubjectIndexEnd()+1
										   // , _vecForward.at(bestForwardIndex).GetSequence()
										   ) 
			    );
	  
	  
	  //tempLstF.SetSequence(list(desc=ForwardNameSet[[bestForwardIndex]], seq=replaceOne);
	  
	  tempLstF.SetSequence(leadingSpaceForward+replaceOne);
	  tempLstF.SetName(_vecForward.at(bestForwardIndex).GetName());
	  
	}
      else
	{
	  //#no need to add leading space
	  tempLstF.SetName("NoMatch");
	  tempLstF.SetSequence("");
	}
      //############here to do!!!!!!!!!
      
      //tempLstF$seq<-paste(leadingSpaceForward, tempLstF$s
      
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
	unsigned int spaceCarryOverFTR=0;//####this is the leading space for reverse adaptor primer, because the insertion in the forward alignment
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
	if(foundForwardFlag)
	  {
	    //cout<<"1.."<<endl;
	    replaceOne=replaceOne.substr(0, bestForwardAlign.GetPatternIndexStart());
	    //#aligned part
	    //cout<<"2.."<<endl;
	    replaceOne.append( bestForwardAlign.GetPattern(true));
	    //#last part unaligned
	    //cout<<"3.."<<endl;
	    //cout<<"\ttempLstSeq.GetSequence().length():"<<tempLstSeq.GetSequence().length()<<endl;
	    //cout<<"\tbestForwardAlign.GetPatternIndexEnd():"<<bestForwardAlign.GetPatternIndexEnd()<<endl;
	    replaceOne.append( tempLstSeq.GetSequence().substr(bestForwardAlign.GetPatternIndexEnd()+1));// nchar(ff[[i]])));
	    //cout<<"4.."<<endl;
	    tempLstSeq.SetSequence(replaceOne);
	    spaceCarryOverFTR=bestForwardAlign.GetPattern(true).length()- bestForwardAlign.GetPattern(false).length();
	    //cout<<"***with gap:"<<bestForwardAlign.GetPattern(true)<<";without gap:"<<bestForwardAlign.GetPattern(false)<<endl;
	    //cout<<"length is "<<bestForwardAlign.GetPattern(true).length()<<":"<<bestForwardAlign.GetPattern(false).length()<<endl;
	    //cout<<"&&&&&&&&&&&&&carry over FTPR is :"<<spaceCarryOverFTR<<endl;
	  }
	string tempStrSpaces;
	
	for(unsigned int m = 0 ; m < spaceCarryOverFTR; m++)
	  {
	    tempStrSpaces.append("-");
	  }
	tempLstR.SetSequence(leadingSpaceOriginal+tempStrSpaces+tempLstR.GetSequence());//<-paste(tempStr, as.character(tempLstR$seq), sep="");
	//tempLstR$seq<-paste(leadingSpaceOriginal, tempStr, sep="");
	tempLstSeq.SetSequence(leadingSpaceOriginal+tempLstSeq.GetSequence());
	
	//#now put the sequences to the correct vectors
	//cout<<"ready to output strings........"<<endl;
	if(foundForwardFlag)
	{
	  if(foundReverseFlag)
	    {
	      /*tempStrSet<-DNAStringSet(ff[[i]]);
		names(tempStrSet)[1]=names(ff)[i];
		mappedReadsBoth<-append(mappedReadsBoth, tempStrSet);
		tempStrSet<-DNAStringSet(tempLstF$seq);
		names(tempStrSet)[1]=tempLstF$desc;
		mappedReadsBoth<-append(mappedReadsBoth, tempStrSet);
		tempStrSet<-DNAStringSet( tempLstR$seq);
		names(tempStrSet)[1]=tempLstR$desc;*/
	       vec_mapBoth.push_back(tempLstSeq); 
	       vec_mapBoth.push_back(tempLstF);
	       vec_mapBoth.push_back(tempLstR);
	       
	       vec_both_len.push_back(_vecSeq.at(i).GetSequence().length());
	      //mappedBothLen_arr<-c(mappedBothLen_arr,tempLen);
	    }
	  else
	    {
	      /*tempStrSet<-DNAStringSet(ff[[i]]);
			names(tempStrSet)[1]=names(ff)[i];
			mappedReadsForward<-append(mappedReadsForward, tempStrSet);
			
			tempStrSet<-DNAStringSet(tempLstF$seq);
			names(tempStrSet)[1]=tempLstF$desc;
			mappedReadsForward<-append(mappedReadsForward, tempStrSet);
			
			tempStrSet<-DNAStringSet( tempLstR$seq);
			names(tempStrSet)[1]=tempLstR$desc;
			mappedReadsForward<-append(mappedReadsForward, tempStrSet);
	      */
	      vec_mapForward.push_back(tempLstSeq); 
	      vec_mapForward.push_back(tempLstF);
	      vec_mapForward.push_back(tempLstR);
	      vec_forward_len.push_back(_vecSeq.at(i).GetSequence().length());
	      //#mappedReadsForward<-c(mappedReadsForward, list(ff[[i]], tempLstF, tempLstR));
		//	mappedForwardLen_arr<-c(mappedForwardLen_arr,tempLen);
	    }
	}
	else
	  {
	    if(foundReverseFlag)
	      {
		/*tempStrSet<-DNAStringSet(ff[[i]]);
		names(tempStrSet)[1]=names(ff)[i];
		mappedReadsReverse<-append(mappedReadsReverse, tempStrSet);
		
		tempStrSet<-DNAStringSet(tempLstF$seq);
		names(tempStrSet)[1]=tempLstF$desc;
		mappedReadsReverse<-append(mappedReadsReverse, tempStrSet);
		
			tempStrSet<-DNAStringSet( tempLstR$seq);
			names(tempStrSet)[1]=tempLstR$desc;
			mappedReadsReverse<-append(mappedReadsReverse, tempStrSet);
		*/
		//	#mappedReadsReverse<-c(mappedReadsReverse, list(ff[[i]], tempLstF, tempLstR));
		//	mappedReverseLen_arr<-c(mappedReverseLen_arr,tempLen);
		vec_mapReverse.push_back(tempLstSeq); 
		vec_mapReverse.push_back(tempLstF);
		vec_mapReverse.push_back(tempLstR);
		
		vec_reverse_len.push_back(_vecSeq.at(i).GetSequence().length());
	      }
	    else
	      {
		/*tempStrSet<-DNAStringSet(ff[[i]]);
		names(tempStrSet)[1]=names(ff)[i];
		mappedReadsNone<-append(mappedReadsNone, tempStrSet);
			
			tempStrSet<-DNAStringSet(tempLstF$seq);
			names(tempStrSet)[1]=tempLstF$desc;
			mappedReadsNone<-append(mappedReadsNone, tempStrSet);
			
			tempStrSet<-DNAStringSet( tempLstR$seq);
			names(tempStrSet)[1]=tempLstR$desc;
			mappedReadsNone<-append(mappedReadsNone, tempStrSet);
		*/
		vec_mapNone.push_back(tempLstSeq); 
		vec_mapNone.push_back(tempLstF);
		vec_mapNone.push_back(tempLstR);
		//#mappedReadsNone<-c(mappedReadsNone, list(ff[[i]], tempLstF, tempLstR));
		//	mappedNoneLen_arr<-c(mappedNoneLen_arr,tempLen);
		vec_none_len.push_back(_vecSeq.at(i).GetSequence().length());
	      }
	}
	
	if(i%numOfSeqsUnit==0) //#write once every 1000 sequences
	{
	  if(vec_mapBoth.size()>0)
	    {
	      WriteFasta(_mapBoth_fname, vec_mapBoth,100, ofstream::app);
	      //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
	      vec_mapBoth.clear();
	      WriteTextFile(_mapBoth_fname+"_state.txt", vec_both_len, ' ', 1,ofstream::app);
	      vec_both_len.clear();
	    }
	  if(vec_mapForward.size()>0)
	    {
	      WriteFasta( _mapForward_fname, vec_mapForward,100, ofstream::app);
	      //#fileCounter_mpForward<-fileCounter_mpForward+1;
	      vec_mapForward.clear();
	      WriteTextFile(_mapForward_fname+"_state.txt", vec_forward_len, ' ', 1, ofstream::app);
	      vec_forward_len.clear();
	    }
	    if( vec_mapReverse.size()>0)
	      {
		WriteFasta(_mapReverse_fname, vec_mapReverse, 100, ofstream::app);
		//#fileCounter_mpReverse<-fileCounter_mpReverse+1;
		vec_mapReverse.clear();
	      WriteTextFile(_mapReverse_fname+"_state.txt", vec_reverse_len, ' ', 1, ofstream::app);
	      vec_reverse_len.clear();
	      }
	    if(vec_mapNone.size()>0)
	      {
		WriteFasta(  _mapNone_fname,vec_mapNone, 100, ofstream::app);
		//#fileCounter_mpNone<-fileCounter_mpNone+1;
		vec_mapNone.clear();
		WriteTextFile(_mapNone_fname+"_state.txt", vec_none_len, ' ', 1, ofstream::app);
		vec_none_len.clear();
	      }

	    //#also write out the stats
	    /*write(mappedBothLen_arr, "mappedBoth_stats.txt", append=TRUE);
	    mappedBothLen_arr<-c();
	    write(mappedForwardLen_arr, "mappedForward_stats.txt", append=TRUE);
	    mappedForwardLen_arr<-c();
	    write(mappedReverseLen_arr, "mappedReverse_stats.txt",append=TRUE);
	    mappedReveseLen_arr<-c();
	    write(mappedNoneLen_arr, "mappedNone_stats.txt",append=TRUE);
	    mappedNoneLen_arr<-c();*/

	    }//end of each 1000 seqs read write
      

    }//end of for loop of sequence data vec
  
  if(vec_mapBoth.size()>0)
    {
      WriteFasta(_mapBoth_fname, vec_mapBoth,100, ofstream::app);
      //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
      vec_mapBoth.clear();
      WriteTextFile(_mapBoth_fname+"_state.txt", vec_both_len, ' ', 1,ofstream::app);
      vec_both_len.clear();
	    }
  if(vec_mapForward.size()>0)
    {
      WriteFasta( _mapForward_fname, vec_mapForward,100, ofstream::app);
	      //#fileCounter_mpForward<-fileCounter_mpForward+1;
      vec_mapForward.clear();
      WriteTextFile(_mapForward_fname+"_state.txt", vec_forward_len, ' ', 1, ofstream::app);
      vec_forward_len.clear();
    }
  if( vec_mapReverse.size()>0)
    {
      WriteFasta(_mapReverse_fname, vec_mapReverse, 100, ofstream::app);
      //#fileCounter_mpReverse<-fileCounter_mpReverse+1;
      vec_mapReverse.clear();
      WriteTextFile(_mapReverse_fname+"_state.txt", vec_reverse_len, ' ', 1, ofstream::app);
      vec_reverse_len.clear();
    }
  if(vec_mapNone.size()>0)
    {
      WriteFasta(  _mapNone_fname,vec_mapNone, 100, ofstream::app);
      //#fileCounter_mpNone<-fileCounter_mpNone+1;
      vec_mapNone.clear();
      WriteTextFile(_mapNone_fname+"_state.txt", vec_none_len, ' ', 1, ofstream::app);
      vec_none_len.clear();
    }

  cout<<endl;
  ofs_both.close();
  ofs_forward.close();
  ofs_reverse.close();
  ofs_none.close();
}


unsigned int CompareStrings(const string& str1, const string& str2)
{
  unsigned int ret=0;
  unsigned int len1, len2, len;
  len1=str1.length();len2=str2.length();
  len=len1;

  
  if(len>len2)
    {
      len=len2;
      ret=len1-len2;
    }
  else
    {
      ret=len2-len1;
    }
  for(unsigned int i=0;i<len;i++)
    {
      
      if(str1.at(i)!=str2.at(i))
	{
	  ret++;
	}
    }
  
  return ret;
}
