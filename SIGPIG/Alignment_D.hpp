#ifndef ALIGNMENT_D_HPP
#define ALIGNMENT_D_HPP
#include <fstream>
#include <cstring>
#include "AlignmentSettings.hpp"

//defining the alignment object for D gene segment

class Alignment_D
{
public:
  Alignment_D();

  //copy constructor
  Alignment_D(const Alignment_D& _aod);

  //assignment operator
  Alignment_D& operator = (const Alignment_D& _aod);

  //virtual destructor, not necessary???
  virtual ~Alignment_D();
  void ResetData();

  string toString() const;
  
  /*initialize the array elements before doing the alignment
   *at this point, we still don't know the numOfAligned in each
   *D allele. so inside this we only intialize the first level
   *n_D_alleles level;
   */
  bool initialize(const unsigned& _n_D_alleles);

  void Serialize(ofstream& _ofs) const;
  void Deserialize(ifstream& _ifs);

  //======>start to define members
  unsigned n_D_alleles;//total number of D alleles
  unsigned D_max_errors;
  unsigned* numOfAligned;//of length n_D_alleles, indicating the length of align
                         //for each D alleles
    
  unsigned** align_length; //2D vector, n_D_alleles x numOfAligned 
  // 1D:numOfGenD alleles; 2D:numOfAligned;
  
  double** score; //same dimension as above, holding the alignment score
  unsigned** n_errors; //same dimension as above
  unsigned*** error_positions;//3D,
  //1D,fixed, numOf alleles; 2D:numOfAligned;
  //3D,having a size of n_errors[i][j];
  
  unsigned** align_position_left;//alleles by numOfAligned
  //only record in here the starting index of Seq. the 
  //starting align position of target is in deletions array

  unsigned** align_position_right;//same as above

  unsigned** deletions_left; //alleles by numOfAligned
  unsigned** deletions_right; //alleles by numOfAligned

  unsigned*** p_region_max_length_left;//allele x numOfAligned x  max_number of deletions
  unsigned*** p_region_max_length_right; //above
  
  unsigned*** excess_error_positions_left;//3D,
                    //alleles x numOfAligned x max_negative_excess_errors

  unsigned*** excess_error_positions_right;//same as above
  
  unsigned* allele_order;//this is sorted allele array, contains all
  //D segments/alleles, but in an order according to the aligned score.
};


//input:
// _seq, the input seq that is being aligned against. constant reference 
//_genDs, the V gen segments arrays.
//_numOfDSegs, the number of D gene segments in the array
//_V_end, the position where the previous V alignemnt stops
//_J_start, the position where the previous J alignment starts, 
//_scoreMatrix, the matrix do the alignment. pointer to the scorematrix
//_D_maximum_deletion, maximum deletion
//_negative_excess_deletion_max, maximum excess deletion, usually 3
//_max_align, the maximum number of alignments for the D align, 200 usually.
//
//output:
// _D, the Alignment_obj pointer V holding the alignment details. The caller need to allocate the memory
//bool, indicating whether the alignment is successful.
bool match_D(const SequenceString& _seq,
	      const GenomicD* _genDs, const unsigned& _numOfDSegs,
	      const unsigned& _V_end, const unsigned& _J_start,
	     const unsigned& _flank_length,
	      const ScoreMatrix* _ScoreMatrix, 
	     const unsigned& _D_maximum_deletion, 
	     const unsigned& _negative_excess_deletion_max, 
	     const unsigned& _max_align,
	     /*output*/ Alignment_D& _D
 );


/*find the error positions for two sequence aligned.
 * the output array has to be allocated by the caller outside
 *output:
 *    return number of errors
 *    error_positions: are the relative error position index,
 *           relative to the beginning of the start of the array
 */
unsigned findErrors(const string& _seq1, const string& _seq2,
	   const unsigned& pos_start1, const unsigned& pos_end1,
		    const unsigned& pos_start2, const unsigned& pos_end2,
		    const unsigned& max_n_errors, const int& _flank_offset,
		    /*output*/unsigned* error_position
		    );

void DeterminePalindromAndExcessError_D
( const SequenceString& _seq, const GenomicD* _genDs,
  /*const unsigned* _ok_order,*/ unsigned** _deletions_left,
  unsigned** _deletions_right,
  const unsigned& _negative_excess_deletions_max, 
  const unsigned& _D_maximum_deletion,
  const unsigned* const * _align_length, const unsigned& _numOfDSegs,
  const unsigned* _numOfAligned, 
  const unsigned* const* _align_position_left, const unsigned* const* _align_position_right,
  /*output*/ unsigned*** _p_region_max_length_left, 
  unsigned*** _p_region_max_length_right,
  unsigned*** _excess_error_positions_left,
  unsigned*** _excess_error_positions_right
  );
#endif
