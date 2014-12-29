#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include <vector>
#include "../SquenceString.hpp"
#include "../score.hpp"

//for each object here, it contains possibly more then one alignment that passed the filter
//so, it is variable for each input sequences. that is why we need use vector

struct Alignment_D
{
  vector<unsigned> numOfAligned;
  vector<vector<unsigned> > align_length; //2D vector, 
  // 1D:numOfGenD alleles; 2D:numOfAligned;
  vector<vector<unsigned> > score; //same dimension as above
  vector<vector<unsigned> > n_errors; //same dimension as above
  vector<vector<vector<unsigned> > > error_positions;//3D,
  //1D,fixed, numOf alleles; 2D:numOfAligned;
  //3D,having a size of n_errors;
  vector<vector<vector<unsigned> > > excess_error_positions_left;// alleles by numOfAligned by max_negative_excess_errors
  vector<vector<vector<unsigned> > > excess_error_positions_right;//same as above
  
  vector<vector<unsigned> > align_position_left;//alleles by numOfAligned
  vector<vector<unsigned> > align_position_right;//same as above

  vector<vector<unsigned> > deletions_left; //alleles by numOfAligned
  vector<vector<unsigned> > deletions_right; //alleles by numOfAligned

  vector<vector<vector<unsigned> > >p_region_max_length_left;//allele by numOfAligned by  max_number of deletions
  vector<vector<vector<unsigned> > > p_region_max_length_right; //above

  unsigned[] allele_order={ 0,
                            1, 2, 3, 4, 5, 6, 7, 8, 9,10,
			   11,12,13,14,15,16,17,18,19,20,
			    21,22,23,24,25,26,27,28,29,30,
			    31,32,33};
  
};

//will implement this later************
//the reason that we want to have a class instead of struct because we want to more control also more abstraction
//defining more accessor method so we can easily change the implementation detail later without changing the
//interface
//***************
struct Alignmet_Object
{
  /*public:
  Alignment_Object();

  ~Alignment_Object();

  //for initializing
  bool Initialize();

  //for populate the object
  bool Populate(const unsigned& _numOfAligned, vector<unsigned>& _align_length, vector<unsigned[2] >& _align_position,
		vector<unsigned>& _min_deletions, vector<unsigned> _n_errors, vector<vector<unsigned> >  

		private:*/
  unsigned numOfAligned;
  vector<unsigned> align_length;
  vector<unsigned [2] > align_position;//aligned position is a vector of 2 positions each alignment. 
                     //first one for the seq, and second one for the genomic sequence 
  vector<unsigned> min_deletions;
  vector<unsigned> n_errors;
  vector< vector<unsigned> > error_positions; //this is a 2D vector, both dimensions are variable size. 1D size is numOfAligned; 
                    //2D size is variable based on the n_errors[i]. in the original Matlab implementation, this is fixed and having 
                    //size of max_allowed_errors

  vector<vector<unsigned> > p_region_max_length;//This one is storing the length of longest half-palindrome for each value of deletions
  //But not all number of deletions is possible for the specific alignment.
  //this is a 2 D vector, both dimensions are variable size. 1D size is numOfAligned;
  //2D size is fixed, having a size of max_deletion+1, 1 plus because it is also possible having zero deletion.  
  
  vector<vector <unsigned> > excess_error_positions;//2D vectors, holding the negative_excess_error positions. it bases on the alignment
  //and counting one round of errors. Again, in this negative excess error, we assume the current deletion (min_deletion) is not the real 
  //deletion. it is because of sequence errors, so we want to figure out the error positions in this case.

  vector<unsigned> alleles_all;//size of numOfAligned. contains the index to the position of genomic template
  vector<unsigned> alleles_from_distinct_genes;//unique version Cof the above alleles_all
  
};

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
	      const GenomicV* _genVs, unsigned _numOfVSegs, 
	      unsigned _V_minimum_alignment_length, unsigned _V_maximum_deletion, 
	      unsigned _nagative_excess_deletion_max, unsigned _V_allowed_errors, unsigned _error_cost,
	      Alignment_Obj* _V
 );


//see above for the definition (match_Vs)
bool match_J(const SequenceString& _seq,
	      const GenomicV* _genJs, unsigned _numOfJSegs, 
	      unsigned _J_minimum_alignment_length, unsigned _J_maximum_deletion, 
	      unsigned _nagative_excess_deletion_max, unsigned _J_allowed_errors, unsigned _error_cost,
	      Alignment_Obj* _J
 );

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
	      const GenomicV* _genDs, unsigned _numOfDSegs,
	      unsigned _V_end, unsigned _J_start,
	      ScoreMatrix* _ScoreMatrix, unsigned _D_maximum_deletion, 
	     unsigned _nagative_excess_deletion_max, unsigned _max_align,
	      Alignment_D*, _D
 );





#endif
