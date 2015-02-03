#ifndef ALIGNMENT_D_HPP
#define ALIGNMENT_D_HPP

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

  string toString();
  
  static unsigned n_D_alleles;//total number of D alleles
  unsigned* numOfAligned;//of length n_D_alleles, indicating the length of align
                         //for each D alleles
    
  unsigned** align_length; //2D vector, n_D_alleles x numOfAligned 
  // 1D:numOfGenD alleles; 2D:numOfAligned;
  
  unsigned** score; //same dimension as above, holding the alignment score
  unsigned** n_errors; //same dimension as above
  unsigned*** error_positions;//3D,
  //1D,fixed, numOf alleles; 2D:numOfAligned;
  //3D,having a size of n_errors[i][j];
  unsigned*** excess_error_positions_left;//3D,
                    //alleles x numOfAligned x max_negative_excess_errors

  unsigned*** excess_error_positions_right;//same as above
  
  unsigned** align_position_left;//alleles by numOfAligned
  unsigned** align_position_right;//same as above

  unsigned** deletions_left; //alleles by numOfAligned
  unsigned** deletions_right; //alleles by numOfAligned

  unsigned*** p_region_max_length_left;//allele x numOfAligned x  max_number of deletions
  unsigned*** p_region_max_length_right; //above

  static unsigned allele_order [];//={ 0,
  //1, 2, 3, 4, 5, 6, 7, 8, 9,10,
  //			   11,12,13,14,15,16,17,18,19,20,
  //			    21,22,23,24,25,26,27,28,29,30,
  //			    31,32,33};
  
};




#endif
