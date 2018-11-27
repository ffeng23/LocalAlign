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

  _vecIndex1.reserve(_vecSeq.size());
  if(dualIndex)
    {
      _vecIndex2.reserve(_vecSeq.size());
    }

  //Start processing sequences and peel of the indexes
  SequenceString seq;
  string seqName;
  string unknown(lenOfBarcode, 'N');
  string index;
  string index2;
  for(unsigned i=0;i<_vecSeq.size();i++)
    {
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
      //now writting out the output
      _vecIndex1.at(i)=SequenceString(seqName, index);
      if(dualIndex)
	{
	  _vecIndex2.at(i)=SequenceString(seqName,index2);
	}
    }
  return _vecSeq.size();
}


/*in this function, we process the index read or sequence for indexes(the actual read
 *from sequencing. Then we will do simply matching/instead of alignment for
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
		     const string& _out_fname//, const string& _map_fname
		      )
{
  return ;
}
