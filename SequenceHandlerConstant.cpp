#include "SequenceHandler.hpp"
#include "FastaHandler.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include "SequenceHandlerConstant.hpp"

#include "OverlapAlignment.hpp"
#include "AlignmentString.hpp"
#include "LocalAlignment.hpp"

//****global variable***********
//variable vector will hold by isotypes the sequence stringa aligned for output
vector<vector<SequenceString> > g_vec_mapBoth;
vector<vector<SequenceString> > g_vec_mapForward;
vector<vector<SequenceString> > g_vec_mapReverse;
vector<vector<SequenceString> > g_vec_mapNone;

vector<vector<SequenceString> > g_vec_mapBoth_trim;
vector<vector<SequenceString> > g_vec_mapForward_trim;
vector<vector<SequenceString> > g_vec_mapReverse_trim;
vector<vector<SequenceString> > g_vec_mapNone_trim;


//variable vector will hold by isotypes the sequence stirng aligned for output 
vector<vector<unsigned int> > g_vec_len_mapBoth;
vector<vector<unsigned int> > g_vec_len_mapForward;
vector<vector<unsigned int> > g_vec_len_mapReverse;
vector<vector<unsigned int> > g_vec_len_mapNone;

//forward primer set
vector<SequenceString> g_vec_primer_isotype; //forward one
vector<SequenceString> g_vec_primer_upm; //reverse one

//by isotype output flag
bool g_by_isotype_flag=false;
bool g_trim_flag=false;

//**********************************************

/*void SetUpByIsotypeOutputFlag(const bool& _f)
{
  g_by_isotype_flag=_f;
}

void SetUpTrimFlag(const bool& _f)
{
  g_trim_flag=_f;
}

//in this function, we also need to set up the primer information
void ProcessAdaptorSequences(const string& _adaptorFile, const string& _barcodeFile, const string& _forwardFile, 
			     const string& _reverseFile, vector<SequenceString>& _vecForward, 
			     vector<SequenceString>& _vecReverse )
{
  //first, we need to read the sequeces
  vector<SequenceString> adaptor;
  vector<SequenceString> barcode;
  //vector<SequenceString> forward;
  //vector<SequenceString> reverse;
  
  //
  ReadFasta(_adaptorFile, adaptor);
  ReadFasta(_barcodeFile, barcode);
  ReadFasta(_forwardFile, g_vec_primer_isotype);
  ReadFasta(_reverseFile, g_vec_primer_upm);

  SequenceString adaptorA(adaptor.at(0).GetName(), adaptor.at(0).GetSequence());
  SequenceString adaptorB(adaptor.at(1).GetName(), adaptor.at(1).GetSequence());

  //go combine them together
  //forward set -- with adaptorA
  string tempStr;
  string tempName;
  for(unsigned int i=0;i<barcode.size();i++)
    {
      for(unsigned int j=0;j<g_vec_primer_isotype.size();j++)
	{
	  tempStr=adaptorA.GetSequence();
	  tempStr+=barcode.at(i).GetSequence();
	  tempStr+=g_vec_primer_isotype.at(j).GetSequence();
	  tempName=adaptorA.GetName();
	  tempName+=barcode.at(i).GetName();
	  tempName+=g_vec_primer_isotype.at(j).GetName();
	  //
	  SequenceString tempSS(tempName, tempStr);
	  _vecForward.push_back(tempSS);
	}
    }

  //reverse set
  //with adaptor B
  for(unsigned int i=0;i<barcode.size();i++)
    {
      for(unsigned int j=0;j<g_vec_primer_upm.size();j++)
	{
	  tempStr=adaptorB.GetSequence();
	  tempStr+=barcode.at(i).GetSequence();
	  tempStr+=g_vec_primer_upm.at(j).GetSequence();
	  tempName=adaptorB.GetName();
	  tempName+=barcode.at(i).GetName();
	  tempName+=g_vec_primer_upm.at(j).GetName();
	  //
	  SequenceString tempSS(tempName, tempStr);
	  _vecReverse.push_back(tempSS);
	}
    }
  //we are seting up the forward primer information, IgG/M/D, we will 
  //need this information to write output by different isotype
  //initialize the vectors
  if(g_by_isotype_flag)
    {
      for(unsigned int i=0;i<g_vec_primer_isotype.size();i++)
	{
	  vector<SequenceString> t_both;
	  g_vec_mapBoth.push_back(t_both);
	  vector<SequenceString> t_forward;
	  g_vec_mapForward.push_back(t_forward);
	  vector<SequenceString> t_reverse;
	  g_vec_mapReverse.push_back(t_reverse);
	  vector<SequenceString> t_none;
	  g_vec_mapNone.push_back(t_none);

	  vector<SequenceString> t_both_t;
	  g_vec_mapBoth_trim.push_back(t_both_t);
	  vector<SequenceString> t_forward_t;
	  g_vec_mapForward_trim.push_back(t_forward_t);
	  vector<SequenceString> t_reverse_t;
	  g_vec_mapReverse_trim.push_back(t_reverse_t);
	  vector<SequenceString> t_none_t;
	  g_vec_mapNone_trim.push_back(t_none_t);
      

	  vector<unsigned int> t_len_both;
	  g_vec_len_mapBoth.push_back(t_len_both);
	  vector<unsigned int> t_len_forward;
	  g_vec_len_mapForward.push_back(t_len_forward);
	  vector<unsigned int> t_len_reverse;
	  g_vec_len_mapReverse.push_back(t_len_reverse);
	  vector<unsigned int> t_len_none;
	  g_vec_len_mapNone.push_back(t_len_none);
	}
    }
  else
    {
      vector<SequenceString> t_both;
      g_vec_mapBoth.push_back(t_both);
      vector<SequenceString> t_forward;
      g_vec_mapForward.push_back(t_forward);
      vector<SequenceString> t_reverse;
      g_vec_mapReverse.push_back(t_reverse);
      vector<SequenceString> t_none;
      g_vec_mapNone.push_back(t_none);
      
      vector<SequenceString> t_both_t;
      g_vec_mapBoth_trim.push_back(t_both_t);
      vector<SequenceString> t_forward_t;
      g_vec_mapForward_trim.push_back(t_forward_t);
      vector<SequenceString> t_reverse_t;
      g_vec_mapReverse_trim.push_back(t_reverse_t);
      vector<SequenceString> t_none_t;
      g_vec_mapNone_trim.push_back(t_none_t);
      
      
      vector<unsigned int> t_len_both;
      g_vec_len_mapBoth.push_back(t_len_both);
      vector<unsigned int> t_len_forward;
      g_vec_len_mapForward.push_back(t_len_forward);
      vector<unsigned int> t_len_reverse;
      g_vec_len_mapReverse.push_back(t_len_reverse);
      vector<unsigned int> t_len_none;
      g_vec_len_mapNone.push_back(t_len_none);
    }

  //for revere/upm arrays
  

}
*/

//the _adaptorname is the name for the aligned sequence we will be 
// _adaptorName is of format "adaptor A:Barcode:IgM/G/D"
static unsigned int LookUpVectorIndex(const string& _adaptorName, vector<SequenceString>& _vecPrimer)
{
  //cout<<"&&&&&&looking up"<< _adaptorName <<endl;
  //use the name as the key
  for(unsigned int i=0;i<_vecPrimer.size();i++)
    {
      //cout<<"\t%%%%%%"<<_vecPrimer.at(i).GetName()<<endl;
      unsigned found;
      found=_adaptorName.find(_vecPrimer.at(i).GetName()); //will return -1( string::npos) if can not find it.
      //cout<<"\t%%%%%%%%%%%"<<(signed)found<<":npos: "<<string::npos<<endl;
      //if(found!=std::string::npos)
      if((signed)found!=-1)
	{
	  return i;
	}
    }
  cout<<"***ERROR: can not find which catogory to store data (IgM/G/D)"<<endl;
  return -1;
}

void MappingConstants(vector<SequenceString>& _vecForward, vector<SequenceString>& _vecReverse, vector<SequenceString>& _vecSeq, 
		     ScoreMatrix* _sm, const double& _gapOpen, const double& _gapExtension,
			 const double& _mismatchRateThreshold, const unsigned _minimumOverlapLength, /*const unsigned int& _offsetForward, const unsigned int& _offsetReverse,*/ 
		     const string& _mapBoth_fname, const string& _mapForward_fname,
		     const string& _mapReverse_fname, const string& _mapNone_fname)
{
  unsigned int numOfSeqsUnit=20000;
  unsigned int timeOfWriting=1;
  
  //now we need to prepare the output vector, now we will prepare the outfile for each entry in _vecForward
  //they are the alleles for different isotypes (IgG/M), different subtypes (IgG1/2/3/4) and alleles IgG1*01/*02/*03, etc
  //
  
  vector<SequenceString>* pt_vec_mapBoth=new vector<SequenceString>[_vecFoward.size()];
  vector<SequenceString>* pt_vec_mapForward=new vector<SequenceString>[_vecFoward.size()];
  vector<SequenceString> vec_mapReverse;
  vector<SequenceString> vec_mapNone;
  
  vector<SequenceString> vec_crossOver;//this is one holding some extra sequences that after alignment, has constant and upm mapping cross over each other.
    // most of time, we don't expect cross over between them, but might have small chances there are. we want to write down them and take a look to see 
    //what is going on.
    
  /*vector<SequenceString> vec_mapBothTrim;
  vector<SequenceString> vec_mapForwardTrim;
  vector<SequenceString> vec_mapReverseTrim;
  vector<SequenceString> vec_mapNoneTrim;
  */

  
  vector<unsigned int>* pt_vec_mapBoth_pos_end=new vector<unsigned int>[_vecForward.size()];
  vector<unsigned int>* pt_vec_mapForward_pos_end=new vector<unsigned int>[_vecForward.size()];
  
  //now, we need to prepare the output files
  
  AlignmentString tempAS;
  AlignmentString* tempAS_arr;
  unsigned int numOfLocalAlignments=10;
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
      bool      foundCrossOverFlag=false;
      //cout<<"Map forward set"<<endl;
      //forward has to be mapped to the beginning!!!
      for(unsigned int j=0;j<_vecForward.size();j++)
	{
	  //cout<<"forward set:"<<j<<endl;
	  //cout<<_vecForward.at(j).toString()<<endl;
	  LocalAlignment lal (&(_vecSeq.at(i)), &(_vecForward.at(j)), _sm, _gapOpen, _gapExtension,1,numOfLocalAlignments);
	  tempAS_arr=ola.GetAlignment();
	  //cout<<"\talignment score:"<<tempAS.GetScore()<<endl;
	  //need to get the mismatch rate
	  strPattern=tempAS_arr[0].GetPattern(true);
	  strSubject=tempAS_arr[0].GetSubject(true);
	  //cout<<"\tstrPattern:"<<strPattern<<endl;
	  //cout<<"\tstrSubject:"<<strSubject<<endl;

	  mismatch_rate=1-CompareStrings(strPattern, strSubject)/((double)strPattern.length());
	  //cout<<"\tcompare:"<<CompareStrings(strPattern, strSubject)<<";length():"<<strPattern.length()<<endl;
	  //cout<<"\tmismatch_rate:"<<mismatch_rate<<endl;
		
	  //check the best score and mismatch rate
	  if(tempAS_arr[0].GetScore()>bestForwardScore&&mismatch_rate>_mismatchRateThreshold&&strPattern.length()>_minimumOverlapLength)
	    {
	      //we need to see the offset too, only important to the pattern, here the 
	      //the pattern is the long one, we align the primer sequence against
	      //the primer should be in the beginning, not too far,
	      //if(tempAS.GetPatternIndexStart()<_offsetForward )
	      //{//we good
		  //cout<<"\t***get one bigger"<<endl;
		  bestForwardScore=tempAS_arr[0].GetScore();
		  bestForwardAlign=tempAS_arr[0];//with name
		  bestForwardIndex=j;
		  foundForwardFlag=true;
		  //}
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
	      //if( 
	      // (_vecSeq.at(i).GetLength()-tempAS.GetPatternIndexEnd())<_offsetReverse )
	      //{//we good
		  bestReverseScore=tempAS.GetScore();
		  bestReverseAlign=tempAS;//with name
		  bestReverseIndex=k;
		  foundReverseFlag=true;
		  //}
	    }
	}

      //cout<<"Done with mapping:i="<<i<<endl;
      //###done with mapping, now we need to check out primer dimer
      

      //for constant region, we only need to get the constant alignment and will write where the constant breaks
      //if there are breaks within the constant region. we are anticipating more to be full length???
      
      //first, we need to check the alignment, ?????local alignment or overlap alignment????? 
      //  **********************************       ---sequence read
      //    ==^^^^^^^^^^^^===                      ---foward primer
      //             =====!!!!!!!!!!!!!==          ----reverse primer
      //
      //= : mismatch
      //* : sequence read
      //^ : forward primer match
      //! : reverse primer match
      
      //so the algorithm is we check the alignment, write down where the constant alignment ends
      //we also want to distinguish between ends within or outside the constant region
      
      //now we want to write down data by each different allele, we can always manually combine them into different categories.

      //********WE NEED THIS LATER**************************************
      /*if(!foundReverseFlag||!foundForwardFlag)
	{
	  continue;
	}
      */
      //now, start the algorithm for the constant regions
      if(foundForwardFlag&&foundReverseFlag&&(bestForwardAlign.GetPatternIndexEnd()<bestReverseAlign.GetPatternIndexStart())) //constant end
	{//here, we might want to write out for debugging, we are not expecting the cross over between this 
	  foundCrossOverFlag=true;
	  //continue;
	}
      //cout<<"\t***found a dimer"<<endl;
      //now we are finding the primer dimer.
      //write down the stats
      //first look up what this sequence is from IgM or IgG or upm
      //isotype first
      
      //unsigned int foundIndex=LookUpVectorIndex(_vecForward.at(bestForwardIndex).GetName(), g_vec_primer_isotype);
	  unsigned int foundIndex=bestForwardIndex;//here we don't need to look up, because we did not precess adaptor/primer, we simply get sort them by each allel
	//here, to keep track of positions, because we are running local alignment, so we can simply keep track the end position of alignment.
	//we also have to check whether 
	 if(foundForwardFlag)
	 {
	   if(foundReverseFlag)
	   {
	     if(!foundCrossOverFlag)
	     {
		pt_vec_mapBoth_pos_end[foundIndex].push_back(bestForwardAlign.GetSubjectIndexEnd());
	     }
	     else
	     {
		//empty
	       //here, for now we only go ahead write out the sequences for debugging, will do something more???
	     }
	   }
	   else
	   {
	     //found forward only, we will put down the stats
	     pt_vec_mapForward_pos_end[foundIndex].push_back(bestForwardAlign.GetSubjectIndexEnd());
	   }
	 }
	 else  //if not found the constant region anyway, we will not write any stats, but simply go aheed for output
	 {
	    //empty
	 }
	
      //we don't save upm primer dimer positions for this
      //foundIndex=LookUpVectorIndex(_vecReverse.at(bestReverseIndex).GetName(), g_vec_primer_upm);
      //pt_vec_upm_primerDimer_pos_cross[foundIndex].push_back(bestReverseAlign.GetSubjectIndexStart()+(bestForwardAlign.GetPatternIndexEnd()-bestReverseAlign.GetPatternIndexStart()));

      //now we are done with stats recording.
      
      //***************************output*****************
      //prepare the aligned for debugging purpose
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
	//cout<<"\tready to output strings........"<<endl;
	vector<SequenceString>* p_vec_map;
	//vector<unsigned int>* p_vec_len_map;
	//vector<SequenceString>* p_vec_map_trim=NULL;
	
	if(foundForwardFlag)
	{
	  if(foundReverseFlag)
	    {
	      //check to see whether we need to store the data by isotype
	      if(g_by_isotype_flag)
		{
		  unsigned int found=LookUpVectorIndex(tempLstF.GetName(), g_vec_primer_isotype);
		  //cout<<"\t**********writing to vec now"<<found<<endl;
		  if(found != string::npos)
		    {
		      p_vec_map=&g_vec_mapBoth.at(found);  //g_vec_mapBoth is initilized in the processingAdaptor function above
		      
		    }
		}
	      
	     
	    }
	  
	}
	
	//now we need to store the data into the correct vectors/files
	p_vec_map->push_back(tempLstSeq); 
	p_vec_map->push_back(tempLstF);
	p_vec_map->push_back(tempLstR);
	
		
	//cout<<"start writing output........"<<endl;
	if(i>=timeOfWriting*numOfSeqsUnit) //#write once every 20000 sequences
	{
	  //cout<<"i round:"<<i<<endl;
	  //mapBoth
	  for(unsigned int s=0;s<g_vec_mapBoth.size();s++)
	    {
	      if(g_vec_mapBoth.at(s).size()>0)
		{
		  //cout<<"\t------writing files at i-----:"<<s<<endl;
		  //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
		  string t_fileName=_mapBoth_fname;
		  if(g_by_isotype_flag)
		    t_fileName=_mapBoth_fname+g_vec_primer_isotype.at(s).GetName();
		  WriteFasta(t_fileName, g_vec_mapBoth.at(s),100, ofstream::app);
		  //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
		  g_vec_mapBoth.at(s).clear();
		  //cout<<"vecBoth cleared:"<<vec_mapBoth.size()<<endl;
		  WriteTextFile(t_fileName+"_primerdimerStat.txt", pt_vec_isotype_primerDimer_pos_cross[s], ' ', 1,ofstream::app);
		  pt_vec_isotype_primerDimer_pos_cross[s].clear();//g_vec_len_mapBoth.at(s).clear();
		  /*if(g_trim_flag)
		    {
		      WriteFasta(t_fileName+"_trim.fas", g_vec_mapBoth_trim.at(s), 100,ofstream::app);
		      g_vec_mapBoth_trim.at(s).clear();
		      }*/
		}
	      
	    }
	  for(unsigned int t=0;t<g_vec_primer_upm.size();t++)
	    {
	      if(pt_vec_upm_primerDimer_pos_cross[t].size()>0)
		{
		  //cout<<"\t.........writing upm"<<endl;
		  string t_fileName=_mapBoth_fname;
		  if(g_by_isotype_flag)
		    t_fileName=_mapBoth_fname+g_vec_primer_upm.at(t).GetName();
		  WriteTextFile(t_fileName+"_primerdimerStat.txt", pt_vec_upm_primerDimer_pos_cross[t], ' ', 1,ofstream::app);
		  pt_vec_upm_primerDimer_pos_cross[t].clear();
		}
	    }
	  //cout<<"done map both"<<endl;
	}//end of each 1000 seqs read write
      
    }//end of for loop of sequence data vec

  //one last writing
  //cout<<"start writing output........"<<endl;
	
  for(unsigned int s=0;s<g_vec_mapBoth.size();s++)
    {
      if(g_vec_mapBoth.at(s).size()>0)
	{
	  //cout<<"\t------writing files lastly:"<<s<<endl;
	  //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
	  string t_fileName=_mapBoth_fname;
	  if(g_by_isotype_flag)
	    t_fileName=_mapBoth_fname+g_vec_primer_isotype.at(s).GetName();
	  WriteFasta(t_fileName, g_vec_mapBoth.at(s),100, ofstream::app);
	  //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
	  //g_vec_mapBoth.at(s).clear();
	  //cout<<"vecBoth cleared:"<<vec_mapBoth.size()<<endl;
	  WriteTextFile(t_fileName+"_primerdimerStat.txt", pt_vec_isotype_primerDimer_pos_cross[s], ' ', 1,ofstream::app);
	  //pt_vec_isotype_primerDimer_pos_cross[s].clear();//g_vec_len_mapBoth.at(s).clear();
	  
	}
      
    }
  for(unsigned int t=0;t<g_vec_primer_upm.size();t++)
    {
      //cout<<"\t%%%%%%%%tryingn......."<<endl;
      if(pt_vec_upm_primerDimer_pos_cross[t].size()>0)
	{
	  //cout<<"....writing upm lastly"<<endl;
	  string t_fileName=_mapBoth_fname;
	  if(g_by_isotype_flag)
		    t_fileName=_mapBoth_fname+g_vec_primer_upm.at(t).GetName();
	  WriteTextFile(t_fileName+"_primerdimerStat.txt", pt_vec_upm_primerDimer_pos_cross[t], ' ', 1,ofstream::app);
	  //pt_vec_upm_primerDimer_pos_cross[t].clear();
	}
    }
  cout<<"done map both"<<endl;
  
  cout<<endl;

}

