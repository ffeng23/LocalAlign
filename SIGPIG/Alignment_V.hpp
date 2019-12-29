#ifndef ALIGNMENT_V_HPP
#define ALIGNMENT_V_HPP

#include "Alignment.hpp"
#include "MatrixFunctions.hpp"
void DeterminePalindromAndExcessError_V
( const SequenceString& _seq, const GenomicV* _genVs, const unsigned* _ok_order,
  const unsigned* _min_deletions, 
  const unsigned& _negative_excess_deletions_max, const unsigned& _V_maximum_deletion,
  const unsigned* _align_length, const unsigned& _numOfAligned, unsigned** _align_positions,
  /*output*/ unsigned** _p_region_max_length, unsigned** _excess_error_position  
  );


//defining the functions doing alignment

//the function to do the V alignments
//input:
// _seq, the input seq that is being aligned against. constant reference 
//_genVs, the V gen segments arrays.
//_numOfVSegs, the number of v gene segments in the array
//_V_minimum_alignment_length, length of v mininum alignement
//_V_maximum_deletion, maximum deletion
//_negative_excess_deletion_max, maximum excess deletion, usually 3
//_v_allowed_errors, number maximum allowed errors in the alignment
//_error_cost, mismatch error cost, -4 or -5 (?)
//
//output:
// _V, the Alignment_obj pointer V holding the alignment details. The caller need to allocate the memory
//bool, indicating whether the alignment is successful.
bool match_V(const SequenceString& _seq,
	      const GenomicV* _genVs, const unsigned& _numOfVSegs, 
	      const unsigned& _V_minimum_alignment_length, const unsigned& _V_maximum_deletion, 
	      const unsigned& _negative_excess_deletion_max, const unsigned& _V_allowed_errors, 
	     const unsigned& _error_cost,
	     /*output*/ Alignment_Object& _V
 );


/* doing the alignment for V with no fixing on either end
//finds highest scoring alignment of between seq and target 
//doing the alignment forcing no sides on either end. 
//Allowing the alignment to start at any nt of target and seq
//
//input:
// _seq, the input sequence to align
// _target, the genomice template to align against.
// _maximium_errors, the maximum number of errors allow in this alignment
// _error_cost, the error cost, double
//
//output: Note, the caller need to initialize the memory
// _align_position, the array of size 2 to contain the alignment starting position 
//         for the best alignment. [0] for _seq, [1] for _target
// _n_errors, the pointer to the unsigned variable for holding the number of errors in the best alignment
// _error_positions, the array of size _maximum_errors. but it only have n_errors elements. Again, need to
//         initialized by the caller. 

//return: 
//   alignment_length, unsigned
*/
unsigned align_with_constraints_fast_no_fix
     (const string& _seq, const string& _target, 
      const unsigned& _maximum_errors, const unsigned& _minimum_align_length,  
      const double& _error_cost,
      /*output*/unsigned* _align_position, unsigned& n_errors, 
      unsigned* _error_positions);

/*the function wrapper for calculate the score of 
 *the alignment.
 *so far the score is simple, but we might want to change
 *later. so make it easy for now.
 * return: the score
 */ 
double CalculateScore
  (const unsigned& _total_aligned_len, 
   const unsigned* _error_positions,
   const unsigned& _n_errors,
   const double& _cost);

/*the function used to figure out the sub_alignment with the best
 *score along the full alignment. Here we add one more parameter
 *allow the subalignment could move freely on the left side, 
 *"_curr_error_position_index" (see below for explanation)
 *
 *input:
 *   _total_align_len:the full length of the aligment we 
 *            are working to get the best sub alignment 
 *   _all_error_positions: the error positoins of the full 
 *            alignment. The positions are relative to the
 *            beginning of the alignment.
 *   _curr_error_position_index: we add this one to make the 
 *            best alignment could start not necessarily with 
 *            the fixed left end. The outer caller can set up 
 *            to make the search from the left, ie varying
 *            left to make it removing left errors.
 *            THIS VALUE COULD BE -1, which
 *            means we started from very beginning of the full
 *            alignment and include all the errors. If it is 0, 
 *            it means we don't include first error and start
 *            to consider the sequence beginning from 
 *            error_position[0]+1 positoins.!!!!!
 *output:
 *   _best_n_errors: the n_errors give the best score. it starts
 *            right after the curr_error_position_index (DON'T
 *            include this one).
 *   NOTE: we don't include the error positoins, since we can
 *         figure it out by other information. 
 * return: the best aligned length 
 */
unsigned  GetBestScoreAmongTheAlignment
   (const unsigned& _total_align_len, 
    const unsigned* _all_error_positions,
    const unsigned& _all_n_errors,  
    const unsigned& _curr_error_position_index,
    const double& _cost,
    /*output*/ 
    unsigned& _best_n_errors  );


//this is exactly identical to GetBestScoreAmongTheAlignment(),
//see above for usage.
unsigned  align_with_constraints_after_left_remove_right_errors
   (const unsigned& _total_align_len, 
    const unsigned* _all_error_positions,
    const unsigned& _all_n_errors,  
    const unsigned& _curr_error_position_index,
    const double& _cost,
    /*output*/ 
    unsigned& _best_n_errors  );

/*align the two strings, but need to find the best scores for
 *the subalignments that could remove both left and right errors
 *
 *input:
 *     seq1 and seq2: two string to align fixed left (both start
 *           at position 0
 *     _cost, error cost
 *output:
 *     _n_errors: number of errors after return
 *     _error_positions: holding the position, need 
 *               to be initialized by the caller and 
 *               the size of _v_allowed_errors
 *     _align_position_start
 */
unsigned align_with_constraints_fixed_left_remove_both_errors
(const string& _seq1, const string& _seq2, 
 const unsigned& _v_allowed_errors, 
 const double& _cost,
 /*output*/ unsigned& _n_errors, unsigned* _error_positions, 
 unsigned& _align_position_start);
#endif
