#include "SequenceHandler.hpp"
#include "FastaHandler.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>

#include "SequenceHandlerBarcode.hpp"
//#include "OverlapAlignment.hpp"
//#include "AlignmentString.hpp"
//#include "LocalAlignment.hpp"

/*Helper function to get index sequences from the sequence name
 *we will pass by sequence the vectors and the read 
 *the indexs from the name of each sequences.
 *we assume the sequences is like 
 *     xxx:xxx..xx:xxxx+xxxxx
 * the last field is indexes. could be dual or single index
 * in either case, we will assume one the R1 sequence data
 *file is good enough for index information and will not
 * need R2 sequence file. For a bad one, we will output
 *When processing the indexes, they could be longer, but
 *if they are shorter, we will pad with Ns in the end(??
 *not sure whether this a good way, but wish we won't
 *see this happens a lot).
 *
 *_vecIndex1 and _vecIndex2 holding the indexes from the sequenes
 *return total number of sequences processed
 *
 */

unsigned int ReadIndexFromSequenceName(vector<SequenceString>& _vecSeq,  /*this is the sequence data r1*/
			      vector<SequenceString>& _vecIndex1, /*this is the sequence index data r1*/
			      vector<SequenceString>& _vecIndex2, /*this is the sequence index data r2*/
			      const unsigned& lenOfBarcode,
			      const bool& dualIndex /*indicating whether to do dual indexes*/
			      )
{
  
  //clear up the spaces and make them ready 
  _vecIndex1.clear();
  _vecIndex2.clear();
  //cout<<"****inside function........_vecSeq.size():"<<_vecSeq.size()<<endl;
  _vecIndex1.reserve(_vecSeq.size());
  if(dualIndex)
    {
      _vecIndex2.reserve(_vecSeq.size());
    }
  //cout<<"****inside function........_vecIndex1.size():"<<_vecIndex1.size()<<endl;

  //Start processing sequences and peel of the indexes
  SequenceString seq;
  string seqName;
  string unknown(lenOfBarcode, 'N');
  string index;
  string index2;
  //cout<<"******before for loop......."<<endl;
  for(unsigned i=0;i<_vecSeq.size();i++)
    {
      //cout<<"\t********loop:"<<i<<endl;
      seq=_vecSeq.at(i);
      //get name
      seqName=seq.GetName();
      //parse it
      unsigned int pt= seqName.find_last_of(':');
      
      if(pt!=string::npos)
	{
	  //can find one 
	  index=seqName.substr(pt+1);
	  if(dualIndex)
	    {
	      //parse again
	      pt=index.find_last_of('+');
	      if(pt==string::npos)
		{
		  //can not find it assuming missing index2
		  index2=unknown;
		}
	      else
		{
		  index2=index.substr(pt+1);
		  index=index.substr(0,pt);
		}
	    }
	}
      else
	{
	  index=unknown;
	  index2=unknown;
	}
      if(index.length()<lenOfBarcode)
	{
	  index.append( lenOfBarcode-index.length(),'N');
	}
      if(dualIndex&&index2.length()<lenOfBarcode)
	{
	  index2.append( lenOfBarcode-index2.length(),'N');
	}
      
      
      //now writting out the output
      _vecIndex1.push_back(SequenceString(seqName, index));
      if(dualIndex)
	{
	  _vecIndex2.push_back(SequenceString(seqName,index2));
	}
    }
  //cout<<"*******end "<<endl;
  return _vecSeq.size();
}


/*in this function, we process the index read or sequence for indexes(the actual read
 *from sequencing). Then we will do simply matching/instead of alignment for
 *demuxing or simply identifying. 
 *the input is from the caller (see the input information in the caller 
 *(NGSMapping_Demux.cpp)
 *The output depends on whether we want to do demuxing. If not,
 *we simply output the stats about how many sequences are 
 *for each index and also the indexes for each barcode
 *. note the latter is different from the demuxing. For
 *demuxing, we mean to have real sequences and put them
 *into different file according to the barcode. The sequences
 *might not contain the index. In either case, we will
 * output "demux'ed" indexes.
 */

void MappingBarcodes(vector<SequenceString>& _vecSeq1, /*this is the sequence data r1*/
		    vector<SequenceString>& _vecSeq2, /*this is the sequence data r2*/
		     vector<SequenceString>& _vecIndex1, /*this is the sequence index data r1*/
		     vector<SequenceString>& _vecIndex2, /*this is the sequence index data r2*/
		     vector<SequenceString>& _vecBarSeq1, /*this is the barcode sequences R1*/
		     vector<SequenceString>& _vecBarSeq2, /*this is the barcode sequences R2*/
		     const unsigned& misMatchNum, /*total number of mismatches allowed*/
		     const bool& pairEnd, /*indicating whether to do pairEnd*/
		     const bool& dualIndex, /*indicating whether to do dual indexes*/
		     const bool& indexFromSeq, /*indicating whether to do index from header of sequenences*/
		     const mapType& type, /*indicating whether it is 5'prime or 3' prime*/
		     const bool& demux,/*indicating whether to do dumux as out, otherwise only output stats*/
		     const MatchMatrix* mm,
			const string& _indexR1_fname,// output file name
		    const string& _indexR2_fname,// output file name
			const string& _seqR1_fname,// output file name, demux
			const string& _seqR2_fname// output file name, demux
		      )
{
  //first get everything ready for input
  //checking the conditions
  unsigned lenOfBarcode=0;
  if(_vecBarSeq1.size()<=0)
    {
      //error
      cout<<"***ERROR: barcodes (Read 1) contain zero sequence"<<endl;
      exit(-1);
    }
  if(dualIndex&&_vecBarSeq2.size()<=0)
    {
      cout<<"***ERROR: barcodes (Read 2) contain zero sequence"<<endl;
      exit(-1);
    }
  lenOfBarcode=_vecBarSeq1.at(0).GetSequence().length();
  if(indexFromSeq)
    {
	  cout<<"Parsing the index from seq file......"<<endl;
      unsigned temp=ReadIndexFromSequenceName(_vecSeq1, _vecIndex1,_vecIndex2,lenOfBarcode,dualIndex);
	  cout<<"\t"<<temp <<" indexes are read ........"<<endl;
	  cout<<"Done!!!"<<endl;
    }
  //
  if(pairEnd&&demux&&_vecSeq2.size()<=0)
    {
      //error, no seqs for 
      cout<<"***ERROR: barcodes (Read 2) contain zero sequence, but is required to use demux and pair end reads. "<<endl;
      exit(-1);
    }
  if(dualIndex&&(_vecIndex1.size()!=_vecIndex2.size()))
    {
      //error
      cout<<"***ERROR: index (Read 1) and read2  contain different number of sequences"<<endl;
      exit(-1);
    }
  if(_vecIndex1.size()<=0)
    {
      //error
      cout<<"***ERROR: index (Read 1) contain zero sequence"<<endl;
      exit(-1);
    }
  if(demux&&pairEnd&&(_vecSeq2.size()!=_vecSeq1.size()))
    {
      //error
      cout<<"***ERROR: sequence R1 and R2 contain different number of sequences, but demux and pairEnd required"<<endl;
      exit(-1);
    }
  if(demux&&(_vecIndex1.size()!=_vecSeq1.size()))
    {
      //error
      cout<<"***ERROR: index (Read 1) and sequence (Read 1) contain different number of  sequences"<<endl;
      exit(-1);
    }

	cout<<"Start doing mapping.........."<<endl;
  //now everything is OK. let's do the job.
  unsigned maxNumSeqToWrite=20000;
  double bestScore=lenOfBarcode+10; //small is better
  unsigned bestBar=0;
  unsigned count=0;
  SequenceString index1;
  SequenceString index2;
  
  //cout<<"before looping, the size of vec bar seq1:"<<_vecBarSeq1.size()<<endl;
  vector<SequenceString>* mapIndex1=new vector<SequenceString>[_vecBarSeq1.size()+1]();
  vector<SequenceString>* mapIndex2=new vector<SequenceString>[_vecBarSeq1.size()+1]();
  
  vector<SequenceString>* mapSeq1=new vector<SequenceString>[_vecBarSeq1.size()+1]();
  vector<SequenceString>* mapSeq2=new vector<SequenceString>[_vecBarSeq1.size()+1]();
  //cout<<"*****testing initialization:"<<endl;
	//cout<<"\t mapIndex1 at 1 size:"<<mapIndex1[1].size()<<endl;
  //plus one for the size is because we need store the unmapped sequences
  int* numOfWritesDone =new int [_vecBarSeq1.size()+1]();
  
  vector<vector<double> > stats(_vecBarSeq1.size()+1);//=new unsigned int [_vecBarSeq1.size()+1];
  for(unsigned int i=0;i<_vecBarSeq1.size()+1;i++)
    {
      numOfWritesDone[i]=0;
	  //vector<double> temp;
	  //temp.push_back(0);
	  stats.at(i).push_back(0);
    }
	
	//cout<<"size of stats:"<<stats.size()<<";size of stats[2]:"<<stats[2].size()<<endl;

  ios_base::openmode mode=ofstream::trunc;
  
  //int numOfWritesDone_unmapped=0;
  for(unsigned i=0;i<_vecIndex1.size();i++)//only loop through read 1 index
    {
		if(i%10000==0)
		{
			cout<<"loop "<< i <<"/"<<_vecIndex1.size()<<"......."<<endl;
		}
		index1=_vecIndex1.at(i);
		//cout<<"\tindex1:"<<index1.toString()<<endl;
		if(dualIndex)
		{
		  index2=_vecIndex2.at(i);
		  //cout<<"\tindex2:"<<index2.toString()<<endl;
		}
		bestScore=lenOfBarcode+10;
		//find the best match barcode
		for(unsigned j=0;j<_vecBarSeq1.size();j++)
		{
			//cout<<"\t*****Borcode "<<j<<":"<< _vecBarSeq1.at(j).toString()<<endl;
		  double temp_score=MatchBarcodes(index1,_vecBarSeq1.at(j), mm);
		  //cout<<"\t\tmatch score:"<<temp_score<<endl;
		  if(dualIndex)
			{
				if(type==FivePrime)
					temp_score+=MatchBarcodes(index2,_vecBarSeq2.at(j), mm);
				else
				{
					temp_score+=MatchBarcodes(index2,ReverseComplement(_vecBarSeq2.at(j)), mm);
				}
				//cout<<"\t*****Borcode "<<j<<":"<< _vecBarSeq2.at(j).toString()<<endl;
				//cout<<"\t\ttempscore:"<<temp_score<<endl;
			}
		  if(temp_score <bestScore)
			{//found a better one.
					bestScore=temp_score;
					bestBar=j;
			}
			//cout<<"Best score so far:"<<bestScore<<endl;
		}
		//now we are done with searching.
		/*if(i%1==0)
		{
			cout<<"\tdone search with score of "<<bestScore <<" index of "<<bestBar<<"......."<<endl;
		}*/
		
		if(bestScore>misMatchNum)//the best one still not good, we don't find a good one
		{
			//cout<<"\t\t XXXXXXnot mapped"<<endl;
			bestBar=_vecBarSeq1.size();
		}
		//else 
		//{
		//	//skip and do the next round
		//	continue;
		//}
		//cout<<"\t\tbest Bar:"<<bestBar<<endl;
		mapIndex1[bestBar].push_back(index1);
		
		if(dualIndex)
		{
			//cout<<"doing dual Index"<<endl;
			mapIndex2[bestBar].push_back(index2);
		}
		
		if(demux)
		{
			mapSeq1[bestBar].push_back(_vecSeq1.at(i));
			if(pairEnd)
			{
				mapSeq2[bestBar].push_back(_vecSeq2.at(i));
			}
		}
		count++;//remember how many we done so far.
		
		//cout<<"---------start checking for writing to disk"<<endl;
		string t_fileNameR1;
		string t_fileNameR2;
		stringstream oss;
		//we need to see whether we need to write out open
		if(count>=maxNumSeqToWrite||i==_vecIndex1.size()-1) //do writing either every 10000 seq or at the last one.
		{//do writing
			//determine the mode of file open
		  for(unsigned int s=0;s<=_vecBarSeq1.size();s++)
	      {
			cout<<"Writing output ........."<<endl;  
			if(mapIndex1[s].size()>0)
			  {
				unsigned runningSum=stats.at(s).at(0); 
				runningSum +=mapIndex1[s].size();
				stats.at(s).at(0)=(double)runningSum;				
				if(numOfWritesDone[s]==0)
				  {
					mode=ofstream::trunc;
					//writeHeader=true;
				  }
				else
				  {
					mode=ofstream::app;
					//writeHeader=false;
				  }
				cout<<"\t------writing files at i-----:"<<s<<endl;
				//cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
				string suffix;
				oss<<s;
				suffix="_bar"+oss.str();
				if(s<_vecBarSeq1.size())
				  {
					t_fileNameR1=_indexR1_fname+ suffix;
					if(dualIndex)
					{
						cout<<"\t\t----dual index write Read2 indexes (barcode mapped)"<<endl;
						//t_fileNameR1=t_fileNameR1+"_"+_vecBarSeq2.at(s).GetSequence();
						t_fileNameR2=_indexR2_fname+suffix;
						//t_fileNameR2=t_fileNameR2+"_"+_vecBarSeq2.at(s).GetSequence();
					}
					//vector<double> temp(pt_vec_map_overlap[s].begin(),pt_vec_map_overlap[s].end());
					//vec_stats[0]= temp;
					//temp.clear();

					//WriteTextTableFile(t_fileName+"_ConstantStat.txt", vec_stats, ' ', writeHeader,mode, header);
					
					//pt_vec_map_overlap[s].clear();				
				  }
				else
				  {
					t_fileNameR1=_indexR1_fname+ "notFound";
					if(dualIndex)
					{
						cout<<"\t\t----dual index write Read2 indexes (not found)"<<endl;
						t_fileNameR2=_indexR2_fname+ "notFound";
					}
				  }
				cout<<"\t\t------ filename:"<< t_fileNameR1<<endl;
				WriteFasta(t_fileNameR1+".fasta", mapIndex1[s],100, mode);
				mapIndex1[s].clear();
				if(dualIndex)
				{
					cout<<"\t\t------ filename(REade2):"<< t_fileNameR2<<endl;

					WriteFasta(t_fileNameR2+".fasta", mapIndex2[s],100, mode);
					mapIndex2[s].clear();
				}
				
				if(demux)  
				{
					if(s<_vecBarSeq1.size())
					{
						t_fileNameR1=_seqR1_fname+ suffix;
						if(pairEnd)
						{
							t_fileNameR2=_seqR2_fname+ suffix;
						}
						/*if(dualIndex)
						{
							t_fileNameR1=t_fileNameR1+"_"+ _vecBarSeq2.at(s).GetSequence();
							if(pairEnd)
							{
								t_fileNameR2=t_fileNameR2+"_"+ _vecBarSeq2.at(s).GetSequence();
							}
						}*/
					}
					else
					{
						t_fileNameR1=_seqR1_fname+ "_notFound";
						if(pairEnd)
						{
							t_fileNameR2=_seqR2_fname+ "_notFound";
						}
					}
					//pt_vec_demux[s].clear();
					//if(
					WriteFasta(t_fileNameR1+".fasta", mapSeq1[s],100, mode);
					mapSeq1[s].clear();
					if(pairEnd)
					{
						WriteFasta(t_fileNameR2+".fasta", mapSeq2[s],100, mode);
						mapSeq2[s].clear();
					}
				}
				
				numOfWritesDone[s]++;
				
			}//if each barcode needs writing  
			oss.str("");//clear out the stringstream for next run of barcode number
		  }//for barcode types loop
		  cout<<"Done writting ......."<<endl;
		  count=0;//reset it for the new round of write
		}
		//cout<<"**************Done for checking index ....."<<endl;
    }//end of for loop of index
	
	cout<<"finish checking for all....ready to do stats........"<<endl;
	//write the stats for the last round
	vector<string> header;
	header.push_back("numOfSeqs");
	WriteTextTableFile(_indexR1_fname+"_Stat1.txt", stats, '\t', true,ios_base::trunc, header);
	
  delete [] numOfWritesDone;
  delete [] mapIndex1;
  delete [] mapIndex2;
  delete [] mapSeq1;
  delete [] mapSeq2;
  return ;
}

unsigned int GetBarcodes(vector<SequenceString>& _vecIndex1, /*this is the sequence index data r1*/
		     vector<SequenceString>& _vecIndex2, /*this is the sequence index data r2*/
		     const bool& dualIndex, /*indicating whether to do dual indexes, Index1 and Index2 are available*/
			 vector<SequenceString>& _vecBarSeq1, /*this is the barcode sequences R1*/
		     vector<SequenceString>& _vecBarSeq2, /*this is the barcode sequences R2*/
			 vector<unsigned>& _vecCount /*show the number of sequences holding the barcodes*/
			 )
{
	/*
	//need to turn the vector Barcodes into list for each insertion.
	//leave the index vector as vector 
	list<SequenceString> listIndex1;
	list<SequenceString> listIndex2;
	for(vector<SequenceString>::it=_vecIndex1.cbegin();it!=cend();it++)
	{
		listIndex1.push_back(*it);
	}
	
	if(dualIndex)
	{
		for(vector<SequenceString>::it=_vecIndex2.cbegin();it!=cend();it++)
		{
			listIndex2.push_back(*it);
		}
	}
	
	list<SequenceString> listBarSeq1;
	list<SequenceString> listBarSeq2;
	list<SequenceString> listCount;*/
	cout<<"processing indexes:"<<flush;
	//go through the indexes and insert into the output Barcode vector
	for(unsigned i=0;i<_vecIndex1.size();i++)
	{
		if(i%10000==0)
		{
			cout<<i<<"/"<<_vecIndex1.size()<<"...."<<flush;
		}
		//check whether it has been doing insertion sort kind of operation
		//working lists of barcodes.
		SequenceString r1(_vecIndex1.at(i));
		SequenceString r2;
		if(dualIndex)
			r2=_vecIndex2.at(i);
		//cout<<"i:"<<i<<endl;
		size_t pos=FindInsertionPosition(r1,r2,dualIndex,
					_vecBarSeq1, _vecBarSeq2, _vecCount);
		//cout<<"\tafter position"<<endl;
		if(pos!=string::npos)
			InsertRecordAt(r1,r2,pos, dualIndex, _vecBarSeq1,_vecBarSeq2,_vecCount);
		//else //otherwise we are, because this one is alreay in the barcode vec, and we have update the count in the FindInsertionPosition function
		//cout<<"\tdone for one loop"<<endl;
	}
	
	/*
	//turning lists back into vectors
	for(list<SequenceString>::it=listBarSeq1.cbegin();it!=cend();it++)
	{
		_vecBarSeq1.push_back(*it);
	}
	
	if(dualIndex)
	{
		for(list<SequenceString>::it=listBarSeq2.cbegin();it!=cend();it++)
		{
			_vecBarSeq2.push_back(*it);
		}
	}
	*/
	return _vecBarSeq1.size();				
}

//finding where to insert the "new" barcode
//returning an index to it. if this is not new, update the count and 
//return string::npos 
size_t FindInsertionPosition(const SequenceString& r1,const SequenceString& r2, 
		     const bool& dualIndex, /*indicating whether to do dual indexes, Index1 and Index2 are available*/
			 vector<SequenceString>& _vecBarSeq1, /*this is the barcode sequences R1*/
		     vector<SequenceString>& _vecBarSeq2, /*this is the barcode sequences R2*/
			 vector<unsigned>& _vecCount 
			 )
{
	//now we first need to compare the sequences and find the location to insert
	size_t start=0;
	size_t end=_vecBarSeq1.size()-1;
	
	size_t middle;//=(start+end)/2;
	string ss;
	size_t in_pos=0;
	
	SequenceString r=r1;
	if(dualIndex)
	{
		r.SetSequence(r1.GetSequence()+r2.GetSequence());
	}
	for(;end>=start&&_vecBarSeq1.size()!=0;)//check for cases where we only have 0 elements
	{						//end==start means we have one element, start ==0 and end ==0. 1 element
							//_listBarSeq1.size()==0, no elements
							//end<start occurs only when we have zero element in here.
							//mean we only take care non-zero case in the loop.
							//we need take care 1 or more element case in here 
		middle=(start+end)/2;//in case one element, middle == start
		//do comparison
		ss.assign(_vecBarSeq1.at(middle).GetSequence());
		if(dualIndex)
		{
			ss.append(_vecBarSeq2.at(middle).GetSequence());
		}
		int c=ss.compare(r.GetSequence());
		if(c==0) //we are done
		{
			//update the
			_vecCount.at(middle)+=1;
			return string::npos;
		}
		else
		{
			if(c<0)
			{
				if(end==start+1)//only possible case where we have more than one elements in the vector 
				{
					//in this case, the middle is the start, this means the inserting string is bigger 
					//(i.e. the ss string is smaller than r). this is most likely the case.
					//now we need to compare this one to the "end" string
					ss.assign(_vecBarSeq1.at(end).GetSequence());
					if(dualIndex)
					{
						ss.append(_vecBarSeq2.at(end).GetSequence());
					}
					//check for "end" string 
					int c_2=ss.compare(r.GetSequence());
					if(c_2==0)
					{
						_vecCount.at(end)+=1;//insert after "end" string 
						return string::npos;
					}
					else //for un-equal cases 
					{
						if(c_2<0)  //the end position string also smaller
						{
							in_pos=end+1;
							break;
						}
						else //c_2>0
						{
							in_pos=end;
							break;
						}
					}						
						
				}
				if(end==start)//only possible case is where we have only 1 element in the vector
				{
					in_pos=middle+1; //also equivalently end+1 or start+1
					break;
				}
				start=middle;
			}
			else //c>0
			{
				if(end==start+1)//only possible case where we have more than one elements in the vector 
				{
					//in this case start==middle. so start is bigger, which means the all bigger 
					//than the inserting string, assuming we are dealing with a sorted vector .
					in_pos=start;
					break;
				}
				if(end==start)//only possible case is where we have only 1 element in the vector
				{
					in_pos=start;
					break;
				}
				//if we are here, we might have many elements.
				end=middle;
			}
		}
	}
	return in_pos;
}

//doing the real insertion at the specific position "index"
//
void InsertRecordAt(const SequenceString& r1, const SequenceString& r2, const size_t& index, 
		     const bool& dualIndex, /*indicating whether to do dual indexes, Index1 and Index2 are available*/
			 vector<SequenceString>& _vecBarSeq1, /*this is the barcode sequences R1*/
		     vector<SequenceString>& _vecBarSeq2, /*this is the barcode sequences R2*/
			 vector<unsigned>& _vecCount)
{
	//in this function, we will for sure insert a new record.
	//do it in the begin
	_vecBarSeq1.push_back(r1);
	
	if(dualIndex)
	{
		_vecBarSeq2.push_back(r2);
	}
	_vecCount.push_back(1);
	
	//cout<<"index:"<<index<<";size:"<<_vecBarSeq1.size()<<endl;
	//now move it. from the back
	for(unsigned i=_vecBarSeq1.size()-1;i>index;i--)
	{
		//move i-1 to i
		_vecBarSeq1.at(i)=_vecBarSeq1.at(i-1);
		if(dualIndex)
			_vecBarSeq2.at(i)=_vecBarSeq2.at(i-1);
		_vecCount.at(i)=_vecCount.at(i-1);
	}
	//insert record to the index
	_vecBarSeq1.at(index)=r1;
	if(dualIndex)
	{
		_vecBarSeq2.at(index)=r2;
	}
	_vecCount.at(index)=1;
	//done
}