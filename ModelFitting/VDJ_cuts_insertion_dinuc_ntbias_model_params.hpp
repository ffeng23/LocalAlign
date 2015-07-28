#ifndef VDJ_CUTS_INSERTION_DINUC_NTBIAS_MODEL_PARAMS_HPP
#define VDJ_CUTS_INSERTION_DINUC_NTBIAS_MODEL_PARAMS_HPP

#include "../SIGPIG/genomicSegments.hpp"
//this is defining the global variables
struct VDJ_cuts_insertion_dinuc_ntbias_model_params
{
  VDJ_cuts_insertion_dinuc_ntbias_model_params
  (
   const GenomicV* _genV, const unsigned& _numV, 
   const GenomicD* _genD, const unsigned& _numD,
   const GenomicJ* _genJ, const unsigned& _numJ);
  //define the members
  unsigned max_assignments; //max no. of assignments to explore during Expectation step of EM algorithm
  unsigned max_insertions;

  unsigned max_V_deletions;
  unsigned max_D_deletions;
  unsigned max_J_deletions;

  //here, the _numV/D/J are the number of distinct genes,
  //not toal number of alleles of genes
  unsigned number_V_genes;
  unsigned number_D_genes;
  unsigned number_J_genes;

  //this is for the maxi number of alleles across all the genes
  unsigned max_V_n_alleles;
  unsigned max_D_n_alleles;
  unsigned max_J_n_alleles;

  unsigned max_excess_V_deletions; //maximum excess deletions to consider, this is the extra number of deletions on top of observed deletion/min_deletion/cut(?)
  unsigned max_excess_D_deletions;
  unsigned max_excess_J_deletions;

  unsigned max_palindrome; //maximum half-palindrome length

  int min_V_cut; //cut is the observed deletion==real deletion+insertion+excess error etc.
  int min_D_cut;
  int min_J_cut;
  int max_V_cut;
  int max_D_cut;
  int max_J_cut;
  
  unsigned negative_excess_deletions_max;// (0), /*not sure why we set it up as 0 in here*/
  unsigned min_J_align_length;//(2), 
  unsigned min_J_assign_length;//(1) 
  unsigned min_V_length;//(20)/*originally in matlab code is 15*/,
  unsigned high_error_region;//(15) /*we PROBABLY will NOT use this one*/,
  bool use_no_D_match_seqs;//(true),

  //these below are the maximum possible read length of the 
  unsigned max_J_depth;
  unsigned max_V_depth;
  unsigned max_D_depth;

  //unsigned read_length;
  unsigned maximum_read_length;
};
//extern VDJ_cuts_insertion_dinuc_ntbias_model_params vdj_mps;
#endif
