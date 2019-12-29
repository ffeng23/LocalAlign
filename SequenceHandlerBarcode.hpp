#ifndef SEQUENCEHANDLER_BARCODE_HPP
#define SEQUENCEHANDLER_BARCODE_HPP

#include <vector>
#include <string>
#include "Accessory/SequenceString.hpp"
#include "Accessory/FileHandler.hpp"

#include "score.hpp"
#include <sstream>

using namespace std;

enum mapType{ FivePrime, ThreePrime};

//IN this function, we take care to do the barcode manipulation
//in here we will separate the sequences (_vecSeq1 and _vecSeq2) into different file based on the barcods 
//these vectors could be empty based on input type
// vector: Seq1 and Seq2 are the sequences Read 1 or 2.
//			index1 and index2 are the read index of data
//          barcode1 and barcode2 are the expected barcode user feed in.
//      only barcode1 is required. all others could be input.
//   IndexFromSeq: indicating whether to get indexes from the header of the sequences.
// 
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
/*For mapType, we really mean to reverse complement Barcode2 to match with index2 sequence
 *Not affecting index1 or barcode 1.
 */  
void MappingBarcodes(vector<SequenceString>& _vecSeq1, /*this is the sequence data r1*/
					 vector<SequenceString>& _vecSeq2, /*this is the sequence data r2*/
		     vector<SequenceString>& _vecIndex1, /*this is the sequence index data r1*/
		     vector<SequenceString>& _vecIndex2, /*this is the sequence index data r2*/
		     vector<SequenceString>& _vecBarSeq1, /*this is the barcode sequences R1*/
		     vector<SequenceString>& _vecBarSeq2, /*this is the barcode sequences R2*/
		     const unsigned& misMatchNum, /*total number of mismatches allowed*/
		     const bool& pairEnd, /*indicating whether to do pairEnd, Seq1 and Seq2 both will be needed*/
		     const bool& dualIndex, /*indicating whether to do dual indexes, Index1 and Index2 are available*/
		     const bool& indexFromSeq, /*indicating whether to do index from header of sequenences*/
		     const mapType& type, /*indicating whether it is 5'prime or 3' prime*/
		     const bool& demux,/*indicating whether to do dumux as out, otherwise only output stats*/
		     const MatchMatrix* mm, /*match matrix for scoring the matching of indexes to barcode*/
			const string& _indexR1_fname,// output file name
		    const string& _indexR2_fname,// output file name
			const string& _seqR1_fname,// output file name, demux
			const string& _seqR2_fname,// output file name, demux
			vector<string>& _vecSeq1_Q, //input vector for quality, could be empty.if no emptry, we will write fastq
			vector<string>& _vecSeq2_Q, //input vector for quality, could be empty.if no emptry, we will write fastq
			vector<string>& _vecIndex1_Q, //input vector for quality, could be empty.if no emptry, we will write fastq
			vector<string>& _vecIndex2_Q	//input vector for quality, could be empty.if no emptry, we will write fastq
		      ); 

//a helper function to peel the indexes from the sequences. should not be called directly by the outside caller.
//
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
				       );

//a helper function to read the barcodes from the file. The barcode files must be of the fasta format.
//
/*Helper function to get barcode sequences from the fasta file.
 *this function is different from reading fasta file into a vector by its checking the 
 *barcode length. If the length is smaller than the specified length, an error will be
 *issued. If the length is bigger, then we will chop it to read only the first specified
 *portion. 
 *
 *Input:_fname , the fasta file name of barcodes. 
		_vectBar the holder of the barcodes 
 *		 lenOfBarcode
 *
 *Output: the number of barcodes read in. 
 */
unsigned int ReadBarcodeFromFile(const string& _fname,
					vector<SequenceString>& _vecBar,  /*this is the sequence data r1*/
			      const unsigned& lenOfBarcode
				);


/*function to collect barcodes from the indexes, assuming we don't know the expected barcodes 
	we will also generate stats/counts of each barcode
	dualIndex, is to indicate whether this is dualindex. if not, the _vecIndex2 and _vecBarSeq2 are empty
	_vecBarSeq1, _vecIndex1 and _vecCount are outputs 
 */					   
unsigned int GetBarcodes(vector<SequenceString>& _vecIndex1, /*this is the sequence index data r1*/
		     vector<SequenceString>& _vecIndex2, /*this is the sequence index data r2*/
		     const bool& dualIndex, /*indicating whether to do dual indexes, Index1 and Index2 are available*/
			 vector<SequenceString>& _vecBarSeq1, /*this is the barcode sequences R1*/
		     vector<SequenceString>& _vecBarSeq2, /*this is the barcode sequences R2*/
			 vector<unsigned>& _vecCount /*show the number of sequences holding the barcodes*/
			 );
//finding where to insert the "new" barcode
//returning an index to it. if this is not new, update the count and 
//return string::npos 
size_t FindInsertionPosition(const SequenceString& r1,const SequenceString& r2, 
		     const bool& dualIndex, /*indicating whether to do dual indexes, Index1 and Index2 are available*/
			 vector<SequenceString>& _vecBarSeq1, /*this is the barcode sequences R1*/
		     vector<SequenceString>& _vecBarSeq2, /*this is the barcode sequences R2*/
			 vector<unsigned>& _vecCount 
			 );

//doing the real insertion in the vectors, need to shift everything to keep the vector sorted.
void InsertRecordAt(const SequenceString& r1, const SequenceString& r2, const size_t& index, 
		     const bool& dualIndex, /*indicating whether to do dual indexes, Index1 and Index2 are available*/
			 vector<SequenceString>& _vecBarSeq1, /*this is the barcode sequences R1*/
		     vector<SequenceString>& _vecBarSeq2, /*this is the barcode sequences R2*/
			 vector<unsigned>& _vecCount);
			 
//for the group of 3 functions below, GetBarcodes2, FindInsertionPosition2, InsertRecordAt2,
//we are doing basically the identically manipulations by insertion sort kind of algorithm
//the things we changed are using extra index array to take care of insertion aiming to 
//increase the performance. In the original function, we using vector, the insertion is 
//slow, therefore we try to copy over manually. But still this is slow. So now we try to
//now add an index array. We don't insert or work on vecotr of SequenceString, but instead
//we do it on array of indexes (unsigned). We pre-allocate the array and using memcpy to 
//"insert/shift" elements. Hope we will make it faster.

unsigned int GetBarcodes2(vector<SequenceString>& _vecIndex1, /*this is the sequence index data r1*/
		     vector<SequenceString>& _vecIndex2, /*this is the sequence index data r2*/
		     const bool& dualIndex, /*indicating whether to do dual indexes, Index1 and Index2 are available*/
			 vector<SequenceString>& _vecBarSeq1, /*this is the barcode sequences R1*/
		     vector<SequenceString>& _vecBarSeq2, /*this is the barcode sequences R2*/
			 vector<unsigned>& _vecCount /*show the number of sequences holding the barcodes*/
			 );

#endif


