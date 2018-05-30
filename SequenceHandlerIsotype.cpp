#include "SequenceHandler.hpp"
#include "FastaHandler.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>

#include "SequenceHandlerIsotype.hpp"
#include "OverlapAlignment.hpp"
#include "AlignmentString.hpp"
#include "LocalAlignment.hpp"



//**********************************************

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
void MappingIsotypes(vector<SequenceString>& _vecSeq, /*this is the sequence data that we want to find isotypes*/
		     vector<SequenceString>& _vecIsotype, /*this is the isotype sequences*/
		     const mapType& type, /*indicating whether it is 5'prime or 3' prime*/
		     ScoreMatrix* _sm, const double& _gapOpen, const double& _gapExtension,
		     const double& _matchRateThreshold, const unsigned _minimumOverlapLength,
		     const unsigned int& _offset, //const unsigned int& _offsetReverse, 
		     const string& _map_fname,
		     const string& _unmap_fname//,
		     )
{
  unsigned int numOfSeqsUnit=20000;
  unsigned int timeOfWriting=1;
  int* numOfWritesDone_mapped =new int [_vecIsotype.size()+1];
  for(unsigned int i=0;i<_vecIsotype.size()+1;i++)
    {
      numOfWritesDone_mapped[i]=0;
    }

  int numOfWritesDone_unmapped=0;
  ios_base::openmode mode=ofstream::trunc; //by default, we need to clear out the file.

  //now we need to prepare the output vector, now we will prepare the outfile for each entry in _vecForward
  //they are the alleles for different isotypes (IgG/M), different subtypes (IgG1/2/3/4) and alleles IgG1*01/*02/*03, etc
  //
  vector<SequenceString>* pt_vec_mapped=new vector<SequenceString>[_vecIsotype.size()+1]; //one more for holding unknown isotype
  //vector<SequenceString>* pt_vec_mapForward=new vector<SequenceString>[_vecForward.size()+1];//one more for holding unkown isotype
  //vector<SequenceString> vec_mapReverse;
  vector<SequenceString> vec_mapNone;
  
  //holding the position where within the constant region the break is.
  vector<unsigned int>* pt_vec_map_pos_end=new vector<unsigned int>[_vecIsotype.size()];
  //vector<unsigned int>* pt_vec_mapForward_pos_end=new vector<unsigned int>[_vecIsotype.size()];

  //vector<unsigned int[2]> vec_mapCrossOver_pos_end;
  //vector<unsigned int>vec_mapBreakOutsideConstant_pos_end;

  //now, we need to prepare the output files
  AlignmentString tempAS;
  AlignmentString* tempAS_arr;
  unsigned int numOfLocalAlignments=5; //we most likely will do overlap alignment. so this might not be necessary
  string strPattern;
  string strSubject;
  
  double match_rate=0.0;
  cout<<"\n";
  for(unsigned int i=0;i<_vecSeq.size();i++)
    {
      //printing progress
      if(i/numOfSeqsUnit *numOfSeqsUnit==i)
	{
	  cout<<"..."<<i<<"/"<<_vecSeq.size();
	  flush(cout);
	  
	  //need to write to file.
	}
      
      //here we separate the two different type of mapping, in case we need to do outputing differently
      double best5PrimeScore= -10000000;
      AlignmentString	best5PrimeAlign;
      double best3PrimeScore= -1000000;
      AlignmentString best3PrimeAlign;

      //The following two are used to remember for the best alignment for the Isotype seq index
      unsigned int 	best5PrimeIndex=0; //what is this for? to remeber index of local alignment? do we need it for overlap alignment
      unsigned int	best3PrimeIndex=0;
	
      bool 	found5PrimeFlag=false;
      bool	found3PrimeFlag=false;
      //cout<<"Map isotype set"<<endl;
      for(unsigned int j=0;j<_vecIsotype.size();j++)
	{
	  //now we need to separate out the two different mapping case, 5' and 3'
	  if(type==FivePrime)
	    {
	      //cout<<"forward set:"<<j<<endl;
	      //cout<<_vecForward.at(j).toString()<<endl;
	      LocalAlignment lal (&(_vecSeq.at(i)), &(_vecIsotype.at(j)), _sm, _gapOpen, _gapExtension,1,numOfLocalAlignments);
	      tempAS_arr=lal.GetAlignmentArr();
	      //cout<<"\talignment score:"<<tempAS.GetScore()<<endl;
	      //need to get the match rate
	      strPattern=tempAS_arr[0].GetPattern(true);
	      strSubject=tempAS_arr[0].GetSubject(true);
	      //cout<<"\tstrPattern:"<<strPattern<<endl;
	      //cout<<"\tstrSubject:"<<strSubject<<endl;

	      match_rate=1-CompareStrings(strPattern, strSubject)/((double)strPattern.length());
	      //cout<<"\tcompare:"<<CompareStrings(strPattern, strSubject)<<";length():"<<strPattern.length()<<endl;
	      //cout<<"\tmatch_rate:"<<match_rate<<endl;
		
	      //check the best score and match rate
	      if(tempAS_arr[0].GetScore()>best5PrimeScore&&match_rate>_matchRateThreshold&&strPattern.length()>_minimumOverlapLength)
		{
		  //we need to see the offset too, only important to the pattern, here the 
		  //the pattern is the long one, we align the primer sequence against
		  //the primer should be in the beginning, not too far,
		  if(tempAS_arr[0].GetPatternIndexStart()<_offset )
		    {//we good
		      //cout<<"\t***get one bigger"<<endl;
		      best5PrimeScore=tempAS_arr[0].GetScore();
		      best5PrimeAlign=tempAS_arr[0];//with name
		      best5PrimeIndex=j;
		      found5PrimeFlag=true;
		    }
		}
	    }
	  else //in this case this is 3' mapping
	    {
	      //cout<<"map reverse set"<<endl;
	      //3' should be mapped on the end of the reads, need to reverse complement the sequence too
	      //for(unsigned int k=0;k< _vecReverse.size();k++)
	      //{
		  //cout<<"start doing k:"<<k<<endl;
		  SequenceString reverseComplementReverse=ReverseComplement(_vecIsotype.at(j));
		  //cout<<"after revcomp"<<endl;
		  //OverlapAlignment ola ( &(_vecSeq.at(i)), &reverseComplementReverse, _sm, _gapOpen, _gapExtension,1);
		  LocalAlignment ola ( &(_vecSeq.at(i)), &reverseComplementReverse, _sm, _gapOpen, _gapExtension,1,numOfLocalAlignments);
		  //cout<<"after alignment:"<<k<<endl;
		  //need to get the match rate
		  tempAS=ola.GetAlignment();
		  //need to get the match rate
		  strPattern=tempAS.GetPattern(true);
		  strSubject=tempAS.GetSubject(true);
		  match_rate=1-CompareStrings(strPattern, strSubject)/((double)strPattern.length());
	  	
		  //#check the best score
		  if(tempAS.GetScore()>best3PrimeScore&&match_rate>_matchRateThreshold&&strPattern.length()>_minimumOverlapLength)
		    {
		      //##we need to see the offset too, only important to the pattern, here the 
		      //	##the pattern is the long one, we align the primer sequence against
		      //	###the primer should on the both ends, not too far,
		      if( 
			 (_vecSeq.at(i).GetLength()-tempAS.GetPatternIndexEnd())<_offset )
		      {//we good
			best3PrimeScore=tempAS.GetScore();
			best3PrimeAlign=tempAS;//with name
			best3PrimeIndex=j;
			found3PrimeFlag=true;
		      }
		    }
	    }//for 3 prime mapping
	}//loop through the isotype sequences 
      
      //cout<<"Done with mapping:i="<<i<<endl;
      //###done with mapping, now we need to check out constant region
      

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
      
      //write down the stats, for where the alignment breaks (the end position).
      //write to the variable of pt_vec_map_pos_end 
      //first look up what this sequence is from IgM or IgG 
      //isotype
      //Here we only list two cases where we got alignment, it could
	 //be the case that we don't have alignemtn, in which case, we simply skip the stats
      //unsigned int foundIndex=LookUpVectorIndex(_vecForward.at(bestForwardIndex).GetName(), g_vec_primer_isotype);
      unsigned int foundIndex=best5PrimeIndex;//here we don't need to look up, because we did not precess adaptor/primer, we simply get sort them by each allel
      if(found3PrimeFlag)
	{
	  foundIndex=best3PrimeIndex;
	}

      //here, to keep track of positions, because we are running local alignment, so we can simply keep track the end position of alignment.
      //we also have to check whether the break is in the constant
	 if(found5PrimeFlag)
	 {	   
	   pt_vec_map_pos_end[foundIndex].push_back(best5PrimeAlign.GetSubjectIndexEnd());
	 }
	 
	 if(found3PrimeFlag)
	   { //this part, we need to test, to see whether this is correct
	     //the rationale is that we revcomp the alignment, then the break part for this revComp is that we get the beginning of the alignment and then used the total length (-1) to get the correct break position (end pos)
	   pt_vec_map_pos_end[foundIndex].push_back(_vecIsotype.at(foundIndex).GetLength()-1 - best3PrimeAlign.GetSubjectIndexStart());
	 }
	 
      //now we are done with stats recording.
      
      //***************************output*****************
      //prepare the aligned for debugging purpose
      //#the forward 
      string leadingSpaceOriginal("");
      string leadingSpace5Prime("");
      string leadingSpace3Prime("");
      string replaceOne;
      SequenceString tempLst5Prime;
      SequenceString tempLst3Prime;
      SequenceString tempLstSeq;
            
      //cout<<"forwardSet Output read"<<endl;
      if(found5PrimeFlag)
	{
	  //###we need to figure out the leading spaces in front of the original sequences
	  unsigned int startOriginal=best5PrimeAlign.GetPatternIndexStart();
	  unsigned int startSubject=best5PrimeAlign.GetSubjectIndexStart();
	  if(startSubject>startOriginal)
	    {
	      //leadingSpaceOriginalsg=as.character();
	      for(unsigned int p=0;p<startSubject-startOriginal;p++)
		{
		  leadingSpaceOriginal.push_back('-');
		}
	    }
	  else
	    {
	      for(unsigned int p=0;p<startOriginal-startSubject;p++)
		{
		  leadingSpace5Prime.push_back('-');
		}
	    }
	  
	  //now we need to add the aligned sequence to replace the original one
	  //#we assume the original one is longer than the aligned one, it has to be
	  // #this is the intial part
	  replaceOne=_vecIsotype.at(best5PrimeIndex).GetSequence().substr(0, best5PrimeAlign.GetSubjectIndexStart());
	  //#aligned part
	  replaceOne.append( best5PrimeAlign.GetSubject(true));
	  //#last part unaligned
	  replaceOne.append(_vecIsotype.at(best5PrimeIndex).GetSequence().substr(best5PrimeAlign.GetSubjectIndexEnd()+1)
			    );
	  	  
	  tempLst5Prime.SetSequence(leadingSpace5Prime+replaceOne);
	  tempLst5Prime.SetName(_vecIsotype.at(best5PrimeIndex).GetName());
	}
      else
	{
	  //#no need to add leading space
	  tempLst5Prime.SetName("NoMatch");
	  tempLst5Prime.SetSequence("");
	}
      //############here to do!!!!!!!!!
                  
      //#the revverse
      //cout<<"reverse set output read"<<endl;
	if(found3PrimeFlag)
	  {
	    //#now we need to add the aligned sequence to replace the original one
	    //#we assume the original one is longer than the aligned one, it has to be
	    //	 #this is the intial part
	    SequenceString rc3PrimeSeq=ReverseComplement(_vecIsotype.at(best3PrimeIndex));
	    replaceOne=rc3PrimeSeq.GetSequence().substr(0, best3PrimeAlign.GetSubjectIndexStart());
	    //cout<<"\t\t*****subject start:"<<bestReverseAlign.GetSubjectIndexStart()<<";subject end:"<<bestReverseAlign.GetSubjectIndexEnd()<<endl;
	    //cout<<"\t\t***replaceOne:"<<replaceOne<<endl;
	    //#aligned part
	    replaceOne.append(best3PrimeAlign.GetSubject(true));
	    //cout<<"\t\t***replaceOne22:"<<replaceOne<<endl;
	    //#last part unaligned
	    replaceOne.append( rc3PrimeSeq.GetSequence().substr(best3PrimeAlign.GetSubjectIndexEnd()+1, rc3PrimeSeq.GetLength()));
	    //cout<<"\t\t***replaceOne333:"<<replaceOne<<endl;
	    tempLst3Prime.SetName(_vecIsotype.at(best3PrimeIndex).GetName());
	    //tempLstR.SetSequence(replaceOne);
	    //#now we need to figure out how the leading space to put in front of reverse one
	    if(best3PrimeAlign.GetPatternIndexStart() >= best3PrimeAlign.GetSubjectIndexStart())
	      {
		unsigned int templen=best3PrimeAlign.GetPatternIndexStart()- best3PrimeAlign.GetSubjectIndexStart();
		for(unsigned int p=0;p<templen;p++)
		  {
		    leadingSpace3Prime.push_back('-');
		  }
		tempLst3Prime.SetSequence(leadingSpace3Prime+replaceOne);
	      }
	    else
	      {
		unsigned int templen=best3PrimeAlign.GetSubjectIndexStart()-best3PrimeAlign.GetPatternIndexStart() ;
		for(unsigned int p=0;p<templen;p++)
		{
		  leadingSpaceOriginal.push_back('-');
		}

		//#here, in this case, the adaptor+primer is longer than the seqs just by alignment, then we need to simply remove some leading part of the adaptor primer
		//tempLst3Prime.SetSequence(replaceOne.substr( best3PrimeAlign.GetSubjectIndexStart()-best3PrimeAlign.GetPatternIndexStart(),
		//					     replaceOne.length()));
	      }
	  }
	else
	  {
	    tempLst3Prime.SetSequence("");
	    tempLst3Prime.SetName("NoMatch");
	  }
	//cout<<"end of reverse on, tempLstR.seq:"<<tempLstR.toString()<<endl;

	//#now we need to take care of the read sequence alignment string
	//#on the reverse part first
	//cout<<"seq string output read...."<<endl;
	//unsigned int spaceCarryOverFTR=0;//####this is the leading space for reverse adaptor primer, because the insertion in the forward alignment
	tempLstSeq.SetSequence(_vecSeq.at(i).GetSequence());
	tempLstSeq.SetName(_vecSeq.at(i).GetName());
	replaceOne=_vecSeq.at(i).GetSequence();
	if(found3PrimeFlag)
	{
	  //cout<<"1.."<<endl;
	  replaceOne=_vecSeq.at(i).GetSequence().substr(0, best3PrimeAlign.GetPatternIndexStart());
	  //#aligned part
	  //cout<<"2.."<<endl;
	  replaceOne.append( best3PrimeAlign.GetPattern(true));
	  //#last part unaligned

	  //cout<<"3..."<<endl;
	  replaceOne.append( _vecSeq.at(i).GetSequence().substr(best3PrimeAlign.GetPatternIndexEnd()+1));//nchar(ff[[i]])), sep="");
	  tempLstSeq.SetSequence(replaceOne);
	  //tempLstSeq.SetName(_vecSeq.at(i).GetName());
	}
	if(found5PrimeFlag)
	  {
	    //cout<<"1.."<<endl;
	    replaceOne=replaceOne.substr(0, best5PrimeAlign.GetPatternIndexStart());
	    //#aligned part
	    //cout<<"2.."<<endl;
	    replaceOne.append( best5PrimeAlign.GetPattern(true));
	    //#last part unaligned
	    //cout<<"3.."<<endl;
	    //cout<<"\ttempLstSeq.GetSequence().length():"<<tempLstSeq.GetSequence().length()<<endl;
	    //cout<<"\tbestForwardAlign.GetPatternIndexEnd():"<<bestForwardAlign.GetPatternIndexEnd()<<endl;
	    replaceOne.append( tempLstSeq.GetSequence().substr(best5PrimeAlign.GetPatternIndexEnd()+1));// nchar(ff[[i]])));
	    //cout<<"4.."<<endl;
	    tempLstSeq.SetSequence(replaceOne);
	    //spaceCarryOverFTR=bestForwardAlign.GetPattern(true).length()- bestForwardAlign.GetPattern(false).length();

	    //cout<<"***with gap:"<<bestForwardAlign.GetPattern(true)<<";without gap:"<<bestForwardAlign.GetPattern(false)<<endl;
	    //cout<<"length is "<<bestForwardAlign.GetPattern(true).length()<<":"<<bestForwardAlign.GetPattern(false).length()<<endl;
	    //cout<<"&&&&&&&&&&&&&carry over FTPR is :"<<spaceCarryOverFTR<<endl;
	  }
	//in case there is no matching on either side, we do nothing, since tempLstSeq has been assigned in the beginning
	
	//string tempStrSpaces;
	
	//for(unsigned int m = 0 ; m < spaceCarryOverFTR; m++)
	//  {
	//    tempStrSpaces.append("-");
	//  }
	//tempLstR.SetSequence(leadingSpaceOriginal+tempLstR.GetSequence());//<-paste(tempStr, as.character(tempLstR$seq), sep="");
	//tempLstR$seq<-paste(leadingSpaceOriginal, tempStr, sep="");
	tempLstSeq.SetSequence(leadingSpaceOriginal+tempLstSeq.GetSequence());
	
	//******************write to vectors***************
	//#now put the sequences to the correct vectors
	//cout<<"\tready to output strings........"<<endl;
	vector<SequenceString> * p_vec_map;
       /*if(foundCrossOverFlag)
	  {
	    p_vec_map=&vec_mapCrossOver;
	    stringstream ss;
	    ss<<tempLstF.GetName()<<":"<<bestForwardAlign.GetSubjectIndexStart()<<":"
	      <<bestForwardAlign.GetSubjectIndexEnd()
	      <<":"<<bestReverseAlign.GetSubjectIndexStart()
	      <<":"<<bestReverseAlign.GetSubjectIndexEnd();
	    tempLstF.SetName(ss.str());
	  }
	else  //not crossover flag
	  {
	    if(foundBreakOutsideConstantFlag)
	      {
		p_vec_map=&vec_mapBreakOutsideConstant;
		stringstream ss;
		ss<<tempLstF.GetName()<<":"<<bestForwardAlign.GetSubjectIndexStart()<<":"
		  <<bestForwardAlign.GetSubjectIndexEnd();

		tempLstF.SetName(ss.str());
	      }
	    else  //'if' foundBreakOutsideConstantFlag
	    {
		if(foundForwardFlag)
		  {
		    unsigned int found=LookUpVectorIndex(tempLstF.GetName(), _vecForward);
		    if(foundReverseFlag)
		      {

			//check to see whether we need to store the data by isotype
			//if(g_by_isotype_flag)
			//	{
			
			//cout<<"\t**********writing to vec now"<<found<<endl;

			//here we always store by isotype
			if(((signed)found) != -1)//string::npos)
			  {
			    p_vec_map=&pt_vec_mapBoth[found];  //g_vec_mapBoth is initilized in the processingAdaptor function above
			  }
			else
			  {
			    cout<<"***WARNING:Can not find the forward primer name, something is wrong"<<endl;
			    p_vec_map=&pt_vec_mapBoth[_vecForward.size()];
			  }
			stringstream ss;
			ss<<tempLstF.GetName()<<":"<<bestForwardAlign.GetSubjectIndexStart()<<":"
			  <<bestForwardAlign.GetSubjectIndexEnd()
			  <<":"<<bestReverseAlign.GetSubjectIndexStart()
			  <<":"<<bestReverseAlign.GetSubjectIndexEnd();
			tempLstF.SetName(ss.str());						
		      }
		    else
		      {
			if(((signed)found) != -1)//string::npos)
			  {
			    p_vec_map=&pt_vec_mapForward[found];  //g_vec_mapBoth is initilized in the processingAdaptor function above
			  }
			else
			  {
			    cout<<"***ERROR:Can not find the forward primer name, something is wrong"<<endl;
			    p_vec_map=&pt_vec_mapForward[_vecForward.size()];
			  }
			stringstream ss;
			ss<<tempLstF.GetName()<<":"<<bestForwardAlign.GetSubjectIndexStart()<<":"
			  <<bestForwardAlign.GetSubjectIndexEnd();
			
			tempLstF.SetName(ss.str());
		      }
		  }
		else
		  {
		    if((foundReverseFlag))
		      {
			p_vec_map=&vec_mapReverse;
		      }
		      else
			p_vec_map=&vec_mapReverse;
		  }
	      }
	  }
	  
	*/
        p_vec_map=&vec_mapNone;
	if(found5PrimeFlag)
	  {
	    unsigned int found=LookUpVectorIndex(tempLst5Prime.GetName(), _vecIsotype);
	    
	    
	    //check to see whether we need to store the data by isotype
	    //if(g_by_isotype_flag)
	    //	{
	    
	    //cout<<"\t**********writing to vec now"<<found<<endl;
	    
	    //here we always store by isotype
	    if(((signed)found) != -1)//string::npos)
	      {
		p_vec_map=&pt_vec_mapped[found];  //g_vec_mapBoth is initilized in the processingAdaptor function above
	      }
	    else
	      {
		cout<<"***WARNING:Can not find the forward primer name, something is wrong"<<endl;
		p_vec_map=&pt_vec_mapped[_vecIsotype.size()];
	      }
	    stringstream ss;
	    ss<<tempLst5Prime.GetName()<<":"<<best5PrimeAlign.GetSubjectIndexStart()<<":"
	      <<best5PrimeAlign.GetSubjectIndexEnd();
	    
	    tempLst5Prime.SetName(ss.str());						
	  }
	
	if(found3PrimeFlag)
	  {
	    unsigned int found=LookUpVectorIndex(tempLst3Prime.GetName(), _vecIsotype);
	    if(((signed)found) != -1)//string::npos)
	      {
		p_vec_map=&pt_vec_mapped[found];  //g_vec_mapBoth is initilized in the processingAdaptor function above
	      }
	    else
	      {
		cout<<"***ERROR:Can not find the forward primer name, something is wrong"<<endl;
		p_vec_map=&pt_vec_mapped[_vecIsotype.size()];
	      }
	    stringstream ss;
	    ss<<tempLst3Prime.GetName()<<":"<<best3PrimeAlign.GetSubjectIndexStart()<<":"
	      <<best3PrimeAlign.GetSubjectIndexEnd();
	    
	    tempLst3Prime.SetName(ss.str());
	  }
	
	//now we need to store the data into the correct vectors/files
	p_vec_map->push_back(tempLstSeq);
	if(found5PrimeFlag)
	  p_vec_map->push_back(tempLst5Prime);
	if(found3PrimeFlag)
	  p_vec_map->push_back(tempLst3Prime);
		
	//cout<<"start writing output........"<<endl;
	//need to figure out how to (trunc or app) write the output at each round.
	
	if(i>=timeOfWriting*numOfSeqsUnit) //#write once every 20000 sequences
	  {
	    
	    //if(timeOfWriting>1)
	    //  {
	    //	mode=ofstream::app;
	    //  }
	    timeOfWriting++;
	    string t_fileName;
	    //cout<<"i round:"<<i<<endl;
	    //mapped
	    for(unsigned int s=0;s<=_vecIsotype.size();s++)
	      {
		if(pt_vec_mapped[s].size()>0)
		  {
		    if(numOfWritesDone_mapped[s]==0)
		      {
			mode=ofstream::trunc;
		      }
		    else
		      {
			mode=ofstream::app;
		      }
		    //cout<<"\t------writing files at i-----:"<<s<<endl;
		    //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
		    if(s<_vecIsotype.size())
		      {
			t_fileName=_map_fname+ _vecIsotype.at(s).GetName()+".fasta";
		      }
		    else
		      {
			t_fileName=_map_fname+ "notFoundIsotype.fasta";
		      }
		    WriteFasta(t_fileName, pt_vec_mapped[s],100, mode);
		    //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
		    pt_vec_mapped[s].clear();
		    //cout<<"vecBoth cleared:"<<vec_mapBoth.size()<<endl;
		    WriteTextFile(t_fileName+"_ConstantStat.txt", pt_vec_map_pos_end[s], ' ', 1,mode);
		    numOfWritesDone_mapped[s]++;
		    pt_vec_map_pos_end[s].clear();//g_vec_len_mapBoth.at(s).clear();
		    /*if(g_trim_flag)
		      {
		      WriteFasta(t_fileName+"_trim.fas", g_vec_mapBoth_trim.at(s), 100,ofstream::app);
		      g_vec_mapBoth_trim.at(s).clear();
		      }*/
		  }
	      }
	    
	    //none
	    if(vec_mapNone.size()>0)
	      {
		//cout<<"\t------writing files at i-----:"<<s<<endl;
		//cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
		if(numOfWritesDone_unmapped==0)
		  {
		    mode=ofstream::trunc;
		  }
		else
		  {
		    mode=ofstream::app;
		  }
		t_fileName=_unmap_fname;
		WriteFasta(t_fileName, vec_mapNone,100, mode);
		numOfWritesDone_unmapped++;
		//#fileCounter_mpBoth<-fileCounter_mpBoth+1;
		vec_mapNone.clear();
	      }
	    
	    //cout<<"done map both"<<endl;
	  }//end of each 1000 seqs read write
      
    }//end of for loop of sequence data vec

  //one last writing
  cout<<"start writing last round of output........"<<endl;
    
  string t_fileName;
  //mapBoth
  for(unsigned int s=0;s<_vecIsotype.size();s++)
    {
      if(pt_vec_mapped[s].size()>0)
	{
	  //cout<<"\t------writing files at i-----:"<<s<<endl;
	  //cout<<"numOfWritesDone:"<<numOfWritesDone_mapped[s]<<endl;
	  if(numOfWritesDone_mapped[s]==0)
	    {
	      mode=ofstream::trunc;
	    }
	  else
	    {
	      mode=ofstream::app;
	    }
	  
	  if(s<_vecIsotype.size())
	    {
	      t_fileName=_map_fname+ _vecIsotype.at(s).GetName()+".fasta";
	    }
	  else
	    {
	      t_fileName=_map_fname+ "notFoundIsotype.fasta";
	    }
	  //cout<<"writing fasta"<<endl;
	  // t_fileName=_mapBoth_fname+ _vecForward.at(s).GetName();
	  WriteFasta(t_fileName, pt_vec_mapped[s],100, mode);
	  //cout<<"writing info"<<endl;
	  //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
	  pt_vec_mapped[s].clear();
	  WriteTextFile(t_fileName+"_ConstantStat.txt", pt_vec_map_pos_end[s], ' ', 1,mode);
	  pt_vec_map_pos_end[s].clear();//g_vec_len_mapBoth.at(s).clear();
	  //cout<<"vecBoth cleared:"<<pt_vec_mapped[s].size()<<endl;
	  //WriteTextFile(t_fileName+"_primerdimerStat.txt", pt_vec_isotype_primerDimer_pos_cross[s], ' ', 1,ofstream::app);
	  //pt_vec_isotype_primerDimer_pos_cross[s].clear();//g_vec_len_mapBoth.at(s).clear();
	  /*if(g_trim_flag)
	    {
	    WriteFasta(t_fileName+"_trim.fas", g_vec_mapBoth_trim.at(s), 100,ofstream::app);
	    g_vec_mapBoth_trim.at(s).clear();
	    }*/
	}
    }
  
  //none
  if(vec_mapNone.size()>0)
    {
      //cout<<"\t------writing files at unmapped last-----:"<<endl;
      //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
      if(numOfWritesDone_unmapped==0)
	{
	  mode=ofstream::trunc;
	}
      else
	{
	  mode=ofstream::app;
	}
      t_fileName=_unmap_fname;
      WriteFasta(t_fileName, vec_mapNone,100, mode);
      //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
      vec_mapNone.clear();
    }
  cout<<"done with one last write"<<endl;
  
  cout<<endl;
  //cleanup

  delete [] pt_vec_mapped;
  //delete [] pt_vec_mapForward;

  delete [] pt_vec_map_pos_end;
  //delete [] pt_vec_mapForward_pos_end;

  delete [] numOfWritesDone_mapped;
}//end of function of map isotype

