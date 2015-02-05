#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include <vector>
#include "../SequenceString.hpp"
#include "../score.hpp"
#include "GenomicV.hpp"
#include "GenomicD.hpp"
#include "GenomicJ.hpp"

//for each object here, it contains possibly more then one alignment that passed the filter
//so, it is variable for each input sequences. that is why we need use vector


//will implement this later************
//the reason that we want to have a class instead of struct because we want to more control also more abstraction
//defining more accessor method so we can easily change the implementation detail later without changing the
//interface
//each object is initialized to be of Fixed elements, numOfAligned is determined by the type of 
//alignment, J or V. In the end, it is not containing all the entry. it is smaller, but the size is still longer
//we use numOfAligned to indicate how many are useful.
//we also need to delete in the constructor
//----
//finally decide to use (array)pointer to hold the data in the alignment object
//need to deconstruct DESTRUCTOR
//-----
//the outer caller will need to initialize and fill the data
//***************
class Alignment_Object
{
  public:
  //default constructor
  //Alignment_Object();
  Alignment_Object(const unsigned& numOfGenTemplates=0);

  //copy constructor
  Alignment_Object(const Alignment_Object& _ao);

  //assignment operator
  Alignment_Object& operator = (const Alignment_Object& _ao);
  
  //virtual destructor, not necessary 
  virtual ~Alignment_Object();
  
  /*
  //for initializing
  bool Initialize();

  //for populate the object
  bool Populate(const unsigned& _numOfAligned, vector<unsigned>& _align_length, vector<unsigned[2] >& _align_position,
		vector<unsigned>& _min_deletions, vector<unsigned> _n_errors, vector<vector<unsigned> >  

		private:*/
  unsigned numOfAligned;
  
  unsigned* align_length;//length of numOfAligned
  unsigned** align_position;//2D array(numOfAligend x 2), aligned position is a vector of 2 positions each alignment. 
                     //first one for the seq, and second one for the genomic sequence 
  unsigned* min_deletions; //len of numOfAligned
  unsigned* n_errors;//len of numOfAligned
  unsigned** error_positions; //this is a 2D array. 1D size is numOfAligned; 
                    //2D size is variable based on the n_errors[i]. in the original Matlab implementation, this is fixed and having 
                    //size of max_allowed_errors
  //for a J alignment, this error position is relative to the end of J chain. and we assume the alignment is 
  //forced to fixed on the end of J chain.!!!

  unsigned** p_region_max_length;//This one is storing the length of longest half-palindrome for each value of deletions
  //But not all number of deletions is possible for the specific alignment.
  //this is a 2 D vector, both dimensions are variable size. 1D size is numOfAligned;
  //2D size is fixed, having a size of max_deletion+1, 1 plus because it is also possible having zero deletion.  
  //again, 2nd Dimension might not be this long since if we don't have that long sequence to be deleted (unlikely), but keep the array the fixed size at this dimension always fixed -- as the maximum_numberOfdeletion allowed
  
  unsigned** excess_error_positions;//2D vectors, holding the negative_excess_error positions. it bases on the alignment
  //and counting one round of errors. Again, in this negative excess error, we assume the current deletion (min_deletion) is not the real 
  //deletion. it is because of sequence errors, so we want to figure out the error positions in this case.
  //2D array, 1D is of numberOfAligned long; 2D is max_nagative_excess_deletions long. it could be smaller (?), example because we don't have that many position, like we only have zero or 1 or 2 deletions. we will not have that many excess error. but in this case we allocate max_negative_error positions, just holding not that many errors.2D is fixed anyway.


  unsigned* alleles_all;//size of numOfAligned. contains the index to the position of genomic template
  unsigned* alleles_from_distinct_genes;//unique version Cof the above alleles_all

  string toString();
  
};



//see above for the definition (match_Vs)
bool match_J(const SequenceString& _seq,
	      const GenomicJ* _genJs, const unsigned& _numOfJSegs, 
	      const unsigned& _J_minimum_alignment_length, const unsigned& _J_maximum_deletion, 
	      const unsigned& _negative_excess_deletion_max, const unsigned& _J_allowed_errors, 
	     const unsigned& _error_cost,
	     /*output*/ Alignment_Object& _J
 );


//finds highest scoring alignment of between seq and target 
//doing the alignment forcing the alignment to begin from the first nucleotide on
//the left of target, but allowing the alignment to start at any nt of target
//note: we modified this code to force to start at first nt on target instead of seq!!!!

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
unsigned align_with_constraints_fast_left(const string& _seq, const string& _target, 
					 const unsigned& _maximum_errors,  const double& _error_cost, 
					 /*output*/unsigned* _align_position, unsigned* n_errors, 
					 unsigned* _error_positions);

//this is the code to finds highest scoring alignement of seq1 and seq2 that forces their
//left ends to match.
//this is the one in the original code named "align_with_constraints_fixed_left_remove_right_errors". 
//we named it so here in order to make the distinction. no further changes made.

//input:
// _seq1 and _seq2, the sequence strings to align, with forcing to begin from the left
// _maximum_errors and error_cost
//output: Note: caller need to initialize the memory for the output
// _n_errors, pointer to the variable for the output,
// _error_positions, array of size maximum_errors, but only the first _n_errors are used.
//return: the align_length.
unsigned align_with_constraints_fixed_left_remove_right_errors
       (const string& _seq1, const string& _seq2, const unsigned& _maximum_errors, 
	const double& _error_cost,
	/*output*/ unsigned* _n_errors, unsigned* _error_positions);

/*input:
 *
 */

void DeterminePalindromAndExcessError_J
( const SequenceString& _seq, const GenomicJ* _genJs, const unsigned* _ok_order,
  const unsigned* _min_deletions, 
  const unsigned& _negative_excess_deletions_max, const unsigned& _J_maximum_deletion,
  const unsigned* _align_length, const unsigned& _numOfAligned, unsigned** _align_positions,
  /*output*/ unsigned** _p_region_max_length, unsigned** _excess_error_position  
  );

void CleanUpMemory(unsigned** p_mem, unsigned size_first_dim);
#endif

