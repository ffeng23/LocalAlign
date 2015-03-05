


#include <string>
#include "AlignmentSettings.hpp"
using namespace std;

//this is the module defining the aligment parameters
namespace AlignmentSettings
{
  unsigned N_per_file = 200;
  //%% Settings and Initializations.

  unsigned J_allowed_errors= 10; //% Maximum number of errors allowed in J match
   unsigned V_allowed_errors= 20; //% Maximum number of errors allowed in V match

  //% V and J matches are sorted by score.
  //% score = match_length - error_cost*n_errors;
  //% I am being lenient with errors at this alignment step.
   double  J_error_cost= 4;
   double  V_error_cost= 4;

  //% Alignment Score Matrix. USED ONLY FOR D MATCHING
  //====>>>>>>>>>we need to take care of this later
   string scorematrix("nuc44");
  //scorematrix=scorematrix(1:4,1:4);
  //scorematrix(scorematrix==-4)=-14;
  //*/

   unsigned max_D_aligns= 200;
  //% Maximum number of deletions from genomic sequence.
   unsigned J_maximum_deletion= 18;
   unsigned V_maximum_deletion= 20;
   unsigned D_maximum_deletion= 16;

  //% Minimum match length for alignment step. This value can be different in
  //% the assignments step.
   unsigned J_minimum_alignment_length= 2;
  //if Read_Length==101
   unsigned V_minimum_alignment_length= 30; //% For 101 bp reads, this is 25, for 60 bp reads, this is 15.
  //else
  //  V_minimum_alignment_length = 15; % For 101 bp reads, this is 25, for 60 bp reads, this is 15.
  //end

  //% To store error and palindrome info for the cases where number of
  //% deletions is less than the 'min_deletions' (from alignments), how far do
  //% you go:
   unsigned negative_excess_deletions_max= 3;
  
  unsigned flank_length=10;
  //unsigned n_D_alleles=34;
  
  unsigned max_length_D_genes=37;

  unsigned max_V_length=270;
  unsigned max_J_length=100;
  unsigned max_D_length=66;
  AlignmentSettings::AlignmentSettings()
  {
    //default one
  }
  AlignmentSettings::~AlignmentSettings()
  {
    //empty one
  }
}



