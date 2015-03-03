

#include "VDJ_cuts_insertion_dinuc_ntbias_model.hpp"

VDJ_cuts_insertion_dinuc_ntbias_model::VDJ_cuts_insertion_dinuc_ntbias_model
(const GenomicV* _genV, const unsigned& _numV, 
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ):
  /*initilization list*/
  max_assignments(6000), max_insertions(30),
  max_V_deletions(16), max_D_deletions(16), max_J_deletions(18),
  max_excess_V_deletions(12),  max_excess_D_deletions(15), max_excess_J_deletions(10),
  max_palindrome(6), 
  min_V_cut(0) /*cut will be set later*/, min_D_cut(0)/*later*/, min_J_cut(0)/*later*/,
  max_V_cut(0), max_D_cut(0),  max_J_cut(0), /*later*/

  /*entropies*/
  S_total(0), S_gene(0), S_V(0), S_DJ(0), S_D(0), S_J(0), S_insVD(0), S_insDJ(0),
  S_insVD_length(0), S_insDJ_length(0), S_insVD_nt(0), S_insDJ_nt(0), S_delV(0),
  S_delD(0), S_delJ(0), S_delDJ(0),

  negative_excess_deletions_max (0), /*not sure why we set it up as 0 in here*/
  min_J_align_length(2), min_J_assign_length(1), 
  min_V_length(20)/*originally in matlab code is 15*/,
  high_error_region(15) /*we PROBABLY will NOT use this one*/,
  use_no_D_match_seqs(true),
  read_length(101)/*the minmum read length has to be 101nts, is this good*/,

/*model parameters*/
  PinsVD(), PinsDJ(),
  RnucleotideVD_per_nucleotideVD_5prime(), RnucleotideDJ_per_nucleotideDJ_3prime(),

  PcutV_given_v(), PcutJ_give_J(), PcutDlcutDr_give_D(),
  PV(), PDJ(), PVallele_given_gen(), PDallele_given_gen(),

  Rerror_per_sequenced_nucleotde()
{
  //need to do things to initialize the parameters and arrays
  min_V_cut=-1*max_palindrome;
  min_D_cut=-1*max_palindrome;
  min_J_cut=-1*max_palindrome;
  
  max_V_cut=max_V_deletions;
  max_D_cut=max_D_deletions;
  max_J_cut=max_J_deletions;

  //start intialize the various matrix
  unsigned dim_size[1]={max_insertions+1};
  
  PinsVD.initialize(1, dim_size, 1);
  PinsDJ.initialize(1, dim_size, 1);

  PV.initialize(1, max_gene_index(genV, _numV), 1);
  unsigned dim_size2[2]={max_gene_index(genD, _numD), max_gene_index(genJ, _numJ)};
  PDJ.initialize(2, dim_size2, 1);

  unsigned dim_size3[2]={max_V_cut-min_V_cut+1, max_gene_index(genV, _numV)};
  PcutV_given_V.initialize(2, dim_size3,1);
}

VDJ_cuts_insertion_dinuc_ntbias_model::~VDJ_cuts_insertion_dinuc_ntbias_model
()
{
  //empty so far
}

bool VDJ_cuts_insertion_dinuc_ntbias_model::ValidateModel()
{
  return true;
}

bool VDJ_cuts_insertion_dinuc_ntbias_model::Normalize()//normalize the model
{

  return true;
}

bool VDJ_cuts_insertion_dinuc_ntbias_model::InitializeCounter(Counter& _c) const
{
  return true;
}

void VDJ_cuts_insertion_dinuc_ntbias_model::GetModelFromCounter(const Counter& _c) 
{
}

void VDJ_cuts_insertion_dinuc_ntbias_model::InitializeAssign(Assigns& _a) const
{
}

void VDJ_cuts_insertion_dinuc_ntbias_model::UpdateCounter(const Assigns& _a, Counter& _c) const
{
}

void VDJ_cuts_insertion_dinuc_ntbias_model::SumCounter(const Counter& _c1, const Counter& _c2, Counter& _retC) const
{
}
