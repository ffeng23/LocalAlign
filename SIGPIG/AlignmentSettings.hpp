#ifndef ALIGNMENT_SETTINGS_HPP
#define ALIGNMENT_SETTINGS_HPP
#include <string>

//this is the module defining the aligment parameters
namespace AlignmentSettings
{
  extern unsigned N_per_file;// = 200;
  //%% Settings and Initializations.

  extern unsigned J_allowed_errors;// = 10; //% Maximum number of errors allowed in J match
  extern unsigned V_allowed_errors;// = 20; //% Maximum number of errors allowed in V match

  //% V and J matches are sorted by score.
  //% score = match_length - error_cost*n_errors;
  //% I am being lenient with errors at this alignment step.
  extern double  J_error_cost;// = 4;
  extern double  V_error_cost;// = 4;

  //% Alignment Score Matrix. USED ONLY FOR D MATCHING
  //====>>>>>>>>>we need to take care of this later
  extern std::string scorematrix;//("nuc44");
  //scorematrix=scorematrix(1:4,1:4);
  //scorematrix(scorematrix==-4)=-14;
  //*/

  extern unsigned max_D_aligns;// = 200;
  //% Maximum number of deletions from genomic sequence.
  extern unsigned J_maximum_deletion;// = 18;
  extern unsigned V_maximum_deletion;// = 20;
  extern unsigned D_maximum_deletion;// = 16;

  //% Minimum match length for alignment step. This value can be different in
  //% the assignments step.
  extern unsigned J_minimum_alignment_length;// = 2;
  //if Read_Length==101
  extern unsigned V_minimum_alignment_length;// = 30; //% For 101 bp reads, this is 25, for 60 bp reads, this is 15.
  //else
  //  V_minimum_alignment_length = 15; % For 101 bp reads, this is 25, for 60 bp reads, this is 15.
  //end

  //% To store error and palindrome info for the cases where number of
  //% deletions is less than the 'min_deletions' (from alignments), how far do
  //% you go:
  extern unsigned negative_excess_deletions_max;// = 3;

  //to define the length of the flank sequences used to run D match
  extern unsigned flank_length;

  //extern unsigned n_D_alleles; we don't set it in advance, we get its value by reading the file of D seg

  extern unsigned max_length_D_genes;//used by match D to determine the max number of errors 
}


#endif
