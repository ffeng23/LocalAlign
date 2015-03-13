#ifndef VDJ_CUTS_INSERTION_DINUC_NTBIAS_ASSIGNS_HPP
#define VDJ_CUTS_INSERTION_DINUC_NTBIAS_ASSIGNS_HPP

#include "Assigns.hpp"
#include "VDJ_cuts_insertion_dinuc_ntbias_model.hpp"
#include "VDJ_cuts_insertion_dinuc_ntbias_counter.hpp"
class VDJ_cuts_insertion_dinuc_ntbias_assigns: public Assigns
{
public: 
  VDJ_cuts_insertion_dinuc_ntbias_assigns(const VDJ_cuts_insertion_dinuc_ntbias_model& _model, const VDJ_cuts_insertion_dinuc_ntbias_counter& _counter);
  virtual ~VDJ_cuts_insertion_dinuc_ntbias_assigns();

  //==========================>member declaration from counter
  Matrix<double> proba; //zeros(model.max_assignments, 1);
  Matrix<double> event_probability; //zeros(model.max_assignments,1);
  unsigned n_assignments;//0, number of valid assignments
  unsigned skips;//0, skipped sequences
  double likelihood;//0, total probability of generating the sequence (sum over ssignments
  double generation_probability;//0, total probability of generating the sequence (assuming zero error rate);
  unsigned max_proba_index;//0, index of the best one
  
  //-------nP****
  //% If counter variable begins with 'nP' i.e. the model variable is a probability distribution P(x),
  //% then the assigns variable should store the index of the value of x, for each assignment for a sequence.
  //% eg. counter.nPpVdelV -> assigns.pVdelV which stores the indices of palindrome half-length and deletions for V as [1 + npV, 1 + ndV]

  //-------nM****
  //% If counter variable begins with 'nM' i.e. the model variable is a mean value of x,
  //% then the assigns variable should store the value of x, for each assignment for a sequence.
  //% eg. counter.nMerror_vs_position -> assigns.error_vs_position which stores number of errors at each position for each assignment for a sequence.
  
  Matrix<unsigned>  insVD;
  Matrix<unsigned>  insDJ;

  Matrix<unsigned>  cutV_given_V;
  Matrix<unsigned>  cutJ_given_J;
  Matrix<unsigned>  cutDlcutDr_given_D;

  Matrix<unsigned>  V;
  Matrix<unsigned>  DJ; //Joint P(V, D, J gene choices)
  Matrix<unsigned>  Vallele_given_gene; //Probabilities of alleles given gene for each gene
  Matrix<unsigned>  Dallele_given_gene;
  Matrix<unsigned>  Jallele_given_gene;


  Matrix<unsigned>  VD_left_edge_dinucleotide;// = zeros(4,4);
  Matrix<unsigned>  VD_right_edge_dinucleotide;// = zeros(4,4);
  Matrix<unsigned>  DJ_left_edge_dinucleotide;// = zeros(4,4);
  Matrix<unsigned>  DJ_right_edge_dinucleotide;// = zeros(4,4);
  
  Matrix<unsigned>  pVmax_delV_V;// =  zeros(model.max_palindrome + 1, model.max_V_deletions + 1, size(model.PV,1));
  Matrix<unsigned>  pJmax_delJ_J;// = zeros(model.max_palindrome + 1, model.max_J_deletions + 1, size(model.PDJ,2));

  Matrix<unsigned>  pDlmax_delDl_D;// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));
  Matrix<unsigned>  pDrmax_delDr_D;// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));

  Matrix<unsigned>  VDJ;// = zeros(size(model.PV,1), size(model.PDJ,1), size(model.PDJ,2));
  Matrix<unsigned>  pVdelV;// = zeros( model.max_palindrome + 1, model.max_V_deletions + 1);
  Matrix<unsigned>  pDldelDl;// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
  Matrix<unsigned>  pDrdelDr;// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
  Matrix<unsigned>  pJdelJ;// = zeros( model.max_palindrome + 1, model.max_J_deletions + 1);


  Matrix<unsigned>  VcutV;// = zeros(size(model.PV,1),model.max_palindrome + model.max_V_deletions + 1);
  Matrix<unsigned>  DcutDl;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<unsigned>  DcutDr;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<unsigned>  JcutJ;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_J_deletions + 1);

  Matrix<unsigned>  DcutV;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_V_deletions + 1);
  Matrix<unsigned>  DcutJ;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_J_deletions + 1);

  Matrix<unsigned>  VcutJ;// = zeros(size(model.PV,1),model.max_palindrome + model.max_J_deletions + 1);
  Matrix<unsigned>  VcutDl;// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<unsigned>  VcutDr;// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);

  Matrix<unsigned>  JcutDl;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<unsigned>  JcutDr;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<unsigned>  JcutV;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_V_deletions + 1);

  Matrix<unsigned>  insVDcutV;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);
  Matrix<unsigned>  insDJcutV;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);

  Matrix<unsigned>  insVDcutDl;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  Matrix<unsigned>  insDJcutDl;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

  Matrix<unsigned>  insVDcutDr;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  Matrix<unsigned>  insDJcutDr;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

  Matrix<unsigned>  insVDcutJ;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);
  Matrix<unsigned>  insDJcutJ;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);

  /*
  unsigned max_V_length;// = model.read_length;
  unsigned max_D_length;// = 16;
  unsigned max_J_length;// = 18;
  */
  Matrix<unsigned>  V_align_length;// = zeros(max_V_length + 1,1);
  Matrix<unsigned>  D_align_length;// = zeros(max_D_length + 1,1);
  Matrix<unsigned>  J_align_length;// = zeros(max_J_length + 1,1);

  Matrix<unsigned>  JJ_align_length;// = zeros(size(model.PDJ,2), 1 + max_J_length);
  Matrix<unsigned>  VV_align_length;// = zeros(size(model.PV,1), 1 + max_V_length);


  Matrix<unsigned>  insDJ_D_align_length;// = zeros(model.max_insertions +1, max_D_length + 1);
  Matrix<unsigned>  insVD_D_align_length;// = zeros(model.max_insertions +1, max_D_length + 1);

  Matrix<unsigned>  insDJ_J_align_length;// = zeros(model.max_insertions +1, max_J_length + 1);
  Matrix<unsigned>  insVD_J_align_length;// = zeros(model.max_insertions +1, max_J_length + 1);

  Matrix<unsigned>  insDJ_V_align_length;// = zeros(model.max_insertions +1, max_V_length + 1);
  Matrix<unsigned>  insVD_V_align_length;// = zeros(model.max_insertions +1, max_V_length + 1);

  Matrix<unsigned>  Dallele_D_align_length;// = zeros(3, max_D_length + 1);


  Matrix<unsigned>  delVinsVD;// = zeros(model.max_V_deletions +1, model.max_insertions +1);
  Matrix<unsigned>  delVinsDJ;// = zeros(model.max_V_deletions +1, model.max_insertions +1);
  Matrix<unsigned>  delVdelDl;// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
  Matrix<unsigned>  delVdelDr;// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
  Matrix<unsigned>  delVdelJ;// = zeros(model.max_V_deletions +1, model.max_J_deletions +1);

  Matrix<unsigned>  delJinsVD;// = zeros(model.max_J_deletions +1, model.max_insertions +1);
  Matrix<unsigned>  delJinsDJ;//zeros(model.max_J_deletions +1, model.max_insertions +1);
  Matrix<unsigned>  delJdelDl;//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
  Matrix<unsigned>  delJdelDr;//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
  
  Matrix<unsigned>  delDlinsVD;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  Matrix<unsigned>  delDlinsDJ;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  Matrix<unsigned>  delDldelDr;//zeros(model.max_D_deletions +1, model.max_D_deletions +1);
  
  Matrix<unsigned>  delDrinsVD;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  Matrix<unsigned>  delDrinsDJ;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  
  Matrix<unsigned>  insVDinsDJ;//zeros(model.max_insertions +1, model.max_insertions +1);
  
  Matrix<unsigned>  VdelV;//zeros(size(model.PV,1), model.max_V_deletions+1 );
  Matrix<unsigned>  DdelDl;//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
  Matrix<unsigned>  DdelDr;//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
  Matrix<unsigned>  JdelJ;//zeros(size(model.PDJ,2), model.max_J_deletions+1 );
  
  Matrix<unsigned>  VinsVD;//zeros(size(model.PV,1), model.max_insertions+1 );
  Matrix<unsigned>  DinsVD;//zeros(size(model.PDJ,1), model.max_insertions+1 );
  Matrix<unsigned>  DinsDJ;//zeros(size(model.PDJ,1), model.max_insertions+1 );
  Matrix<unsigned>  JinsDJ;//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
  Matrix<unsigned>  VdelDl;//zeros(size(model.PV,1), model.max_D_deletions+1 );
  Matrix<unsigned>  VdelDr;//zeros(size(model.PV,1), model.max_D_deletions+1 );
  Matrix<unsigned>  VdelJ;//zeros(size(model.PV,1), model.max_J_deletions+1 );
  
  Matrix<unsigned>  JdelV;//zeros(size(model.PDJ,2), model.max_V_deletions+1 );
  Matrix<unsigned>  JdelDl;//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
  Matrix<unsigned>  JdelDr;//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
  
  Matrix<unsigned>  DdelV;//zeros(size(model.PDJ,1), model.max_V_deletions+1 );
  Matrix<unsigned>  DdelJ;//zeros(size(model.PDJ,1), model.max_J_deletions+1 );
  
  Matrix<unsigned>  VinsDJ;//zeros(size(model.PV,1), model.max_insertions+1 );
  Matrix<unsigned>  JinsVD;//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
  Matrix<unsigned>  pVinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<unsigned>  pVinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<unsigned>  pVdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<unsigned>  pVdelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<unsigned>  pVdelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  Matrix<unsigned>  VpV;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<unsigned>  JpV;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<unsigned>  DpV;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<unsigned>  pJinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<unsigned>  pJinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<unsigned>  pJdelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  Matrix<unsigned>  pJdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<unsigned>  pJdelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<unsigned>  VpJ;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<unsigned>  JpJ;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<unsigned>  DpJ;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<unsigned>  pDlinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<unsigned>  pDlinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<unsigned>  pDldelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  Matrix<unsigned>  pDldelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  Matrix<unsigned>  pDldelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<unsigned>  VpDl;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<unsigned>  JpDl;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<unsigned>  DpDl;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<unsigned>  pDrinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<unsigned>  pDrinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<unsigned>  pDrdelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  Matrix<unsigned>  pDrdelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  Matrix<unsigned>  pDrdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<unsigned>  VpDr;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<unsigned>  JpDr;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<unsigned>  DpDr;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<unsigned>  pVpDl;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<unsigned>  pVpDr;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<unsigned>  pVpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<unsigned>  pDlpDr;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<unsigned>  pDlpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<unsigned>  pDrpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  
  //%------nM*** fields
    
  Matrix<double> mononucleotideVD;// zeros(4,1);
  Matrix<double> mononucleotideDJ;// = zeros(4,1);
  Matrix<double> insertionVD;// = zeros(4,1);
  Matrix<double> insertionDJ;// = zeros(4,1);

  
  Matrix<double> trinucleotideVD;// = zeros(4,4,4);
  Matrix<double> trinucleotideDJ;// = zeros(4,4,4);

  Matrix<double> VV_err_pos;// = zeros(size(model.PV,1), max_V_length + model.max_V_deletions);
  Matrix<double> JJ_err_pos;// = zeros(size(model.PDJ,2), max_J_length);

  Matrix<double> zeroD ;
  Matrix<double> nucleotideVD; Matrix<double> nucleotideVD_5prime;//nucleotide distr's
  Matrix<double> nucleotideDJ; Matrix<double> nucleotideDJ_3prime;

  Matrix<double> error; Matrix<double> sequenced_nucleotide ;//error rate, will be 


  Matrix<double> error_vs_position;//!!this is nM zeros(model.read_length,1);
  Matrix<double> coverage;//zeros(model.read_length,1);
  
};

#endif
