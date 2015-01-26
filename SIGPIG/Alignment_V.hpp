#ifndef ALIGNMENT_V_HPP
#define ALIGNMENT_V_HPP

#include "Alignment.hpp"

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
	     /*output*/ Alignment_Object* _V
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
unsigned align_with_constraints_fast_no_fix(const string& _seq, const string& _target, 
					 const unsigned& _maximum_errors,  const double& _error_cost, 
					 /*output*/unsigned* _align_position, unsigned* n_errors, 
					 unsigned* _error_positions);

#endif
