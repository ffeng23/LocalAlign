#ifndef VDJ_CUTS_INSERTION_DINUC_NTBIAS_MODEL_HPP
#define VDJ_CUTS_INSERTION_DINUC_NTBIAS_MODEL_HPP

#include "BaseModel.hpp"

//inherited class
class VDJ_cuts_insertion_dinuc_ntbias_model:public BaseModel
{
public:
  VDJ_cuts_insertion_dinuc_ntbias_model();

  virtual ~VDJ_cuts_insertion_dinuc_ntbias_model();

  virtual bool ValidateModel();//done know what to validate.
  //check for some "important" model parameters and make sure they are set!!!
  //who know what I will be check.
  //the inherited class will know what to do.
  
  virtual bool Normalize();//normalize the model
  //the inherited class will decide what to do.

  //initialize the input counter, set up its parameter and make it ready to
  //do the model counting.
  //the outer caller need to make the couter available as well as
  //the counter is base class and will be determine by the 
  //inherited call to decide which one to use
  virtual bool InitializeCounter(Counter& _c) const;
  
  //the user will supply the model to be populated
  //polymorphism here!! 
  virtual void GetModelFromCounter(const Counter& _c); 
  
  //initialize the assign and set up the parameter
  //the caller needs to make the assigns available
  //inherited class will decide which one to work on
  //polymorphism
  virtual void InitializeAssign(Assigns& _a) const;
  
  //_c is output
  virtual void UpdateCounter(const Assigns& _a, Counter& _c) const;

  //sum counter
  virtual void SumCounter(const Counter& _c1, const Counter& _c2, Counter& _retC) const;

  //==================================
  //define the members
  unsigned max_assignments; //max no. of assignments to explore during Expectation step of EM algorithm
  unsigned max_insertions;

  unsigned max_V_deletions;
  unsigned max_D_deletions;
  unsigned max_J_deletions;

  unsigned max_excess_V_deletions; //maximum excess deletions to consider, this is the extra number of deletions on top of observed deletion/min_deletion/cut(?)
  unsigned max_excess_D_deletions;
  unsigned max_excess_J_deletions;

  unsigned max_plindrome; //maximum half-palindrome length

  unsigned min_V_cut; //cut is the observed deletion==real deletion+insertion+excess error etc.
  unsigned min_D_cut;
  unsigned min_J_cut;
  unsigned max_V_cut;
  unsigned max_D_cut;
  unsigned max_J_cut;

  //various entropies of the features from the model
  double S_total;
  double S_gene;
  double S_V;
  double S_DJ;
  double S_D;
  double S_J;
  double S_insVD;
  double S_insDJ;
  double S_insVD_length;
  double S_insDJ_length;
  double S_insVD_nt;
  double S_insDJ_nt;
  double S_delV;
  double S_delD;
  double S_delJ;
  double S_delDJ;

  unsigned negative_excess_deletions_max;

  //when considering more deletions than implied by alignment,
  //what are the minimum V, D and J lengths to maintain.
  unsigned min_J_align_length;
  unsigned min_J_assign_length;
  unsigned min_V_length;

  //region on left of sequence that is discarded
  //because of high position dependent error rate
  unsigned high_error_region; //we PROBABLY will NOT use this one.

  //where sequences with too short D alignments are even considered
  unsigned use_no_D_match_seqs;

  //Read length in data set, usually set by main model fitting script
  //we don't use this as a fixed param, but instead, we MIGHT use
  //it as a cutoff value
  unsigned read_length;

  //model parameters
  double PinsVD; //P(insertion)
  double PinsDJ;

  double RnucleotideVD_per_nucleotideVD_5prime;//nucleotide distr's
  double RnucleotideDJ_per_nucleotideDJ_3prime;

  double PcutV_given_v;
  double PcutJ_give_J;
  double PcutDlcutDr_give_D;

  double PV;
  double PDJ; //Joint P(V, D, J gene choices)
  double PVallele_given_gen; //Probabilities of alleles given gene for each gene
  double PDallele_given_gen;

  double Rerror_per_sequenced_nucleotde ;//error rate
  
};


#endif

