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
  unsigned skip;//0, skipped sequences
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
  
  Matrix<double>  insVD;
  Matrix<double>  insDJ;

  Matrix<double>  cutV_given_V;
  Matrix<double>  cutJ_given_J;
  Matrix<double>  cutDlcutDr_given_D;

  Matrix<double>  V;
  Matrix<double>  DJ; //Joint P(V, D, J gene choices)
  Matrix<double>  Vallele_given_gene; //Probabilities of alleles given gene for each gene
  Matrix<double>  Dallele_given_gene;


  Matrix<double>  VD_left_edge_dinucleotide;// = zeros(4,4);
  Matrix<double>  VD_right_edge_dinucleotide;// = zeros(4,4);
  Matrix<double>  DJ_left_edge_dinucleotide;// = zeros(4,4);
  Matrix<double>  DJ_right_edge_dinucleotide;// = zeros(4,4);
  
  Matrix<double>  pVmax_delV_V;// =  zeros(model.max_palindrome + 1, model.max_V_deletions + 1, size(model.PV,1));
  Matrix<double>  pJmax_delJ_J;// = zeros(model.max_palindrome + 1, model.max_J_deletions + 1, size(model.PDJ,2));

  Matrix<double>  pDlmax_delDl_D;// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));
  Matrix<double>  pDrmax_delDr_D;// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));

  Matrix<double>  VDJ;// = zeros(size(model.PV,1), size(model.PDJ,1), size(model.PDJ,2));
  Matrix<double>  pVdelV;// = zeros( model.max_palindrome + 1, model.max_V_deletions + 1);
  Matrix<double>  pDldelDl;// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
  Matrix<double>  pDrdelDr;// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
  Matrix<double>  pJdelJ;// = zeros( model.max_palindrome + 1, model.max_J_deletions + 1);


  Matrix<double>  VcutV;// = zeros(size(model.PV,1),model.max_palindrome + model.max_V_deletions + 1);
  Matrix<double>  DcutDl;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double>  DcutDr;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double>  JcutJ;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_J_deletions + 1);

  Matrix<double>  DcutV;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_V_deletions + 1);
  Matrix<double>  DcutJ;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_J_deletions + 1);

  Matrix<double>  VcutJ;// = zeros(size(model.PV,1),model.max_palindrome + model.max_J_deletions + 1);
  Matrix<double>  VcutDl;// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double>  VcutDr;// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);

  Matrix<double>  JcutDl;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double>  JcutDr;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double>  JcutV;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_V_deletions + 1);

  Matrix<double>  insVDcutV;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);
  Matrix<double>  insDJcutV;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);

  Matrix<double>  insVDcutDl;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double>  insDJcutDl;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

  Matrix<double>  insVDcutDr;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double>  insDJcutDr;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

  Matrix<double>  insVDcutJ;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);
  Matrix<double>  insDJcutJ;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);

  /*
  unsigned max_V_length;// = model.read_length;
  unsigned max_D_length;// = 16;
  unsigned max_J_length;// = 18;
  */
  Matrix<double>  V_align_length;// = zeros(max_V_length + 1,1);
  Matrix<double>  D_align_length;// = zeros(max_D_length + 1,1);
  Matrix<double>  J_align_length;// = zeros(max_J_length + 1,1);

  Matrix<double>  JJ_align_length;// = zeros(size(model.PDJ,2), 1 + max_J_length);
  Matrix<double>  VV_align_length;// = zeros(size(model.PV,1), 1 + max_V_length);


  Matrix<double>  insDJ_D_align_length;// = zeros(model.max_insertions +1, max_D_length + 1);
  Matrix<double>  insVD_D_align_length;// = zeros(model.max_insertions +1, max_D_length + 1);

  Matrix<double>  insDJ_J_align_length;// = zeros(model.max_insertions +1, max_J_length + 1);
  Matrix<double>  insVD_J_align_length;// = zeros(model.max_insertions +1, max_J_length + 1);

  Matrix<double>  insDJ_V_align_length;// = zeros(model.max_insertions +1, max_V_length + 1);
  Matrix<double>  insVD_V_align_length;// = zeros(model.max_insertions +1, max_V_length + 1);

  Matrix<double>  Dallele_D_align_length;// = zeros(3, max_D_length + 1);


  Matrix<double>  delVinsVD;// = zeros(model.max_V_deletions +1, model.max_insertions +1);
  Matrix<double>  delVinsDJ;// = zeros(model.max_V_deletions +1, model.max_insertions +1);
  Matrix<double>  delVdelDl;// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
  Matrix<double>  delVdelDr;// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
  Matrix<double>  delVdelJ;// = zeros(model.max_V_deletions +1, model.max_J_deletions +1);

  Matrix<double>  delJinsVD;// = zeros(model.max_J_deletions +1, model.max_insertions +1);
  Matrix<double>  delJinsDJ;//zeros(model.max_J_deletions +1, model.max_insertions +1);
  Matrix<double>  delJdelDl;//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
  Matrix<double>  delJdelDr;//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
  
  Matrix<double>  delDlinsVD;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  Matrix<double>  delDlinsDJ;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  Matrix<double>  delDldelDr;//zeros(model.max_D_deletions +1, model.max_D_deletions +1);
  
  Matrix<double>  delDrinsVD;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  Matrix<double>  delDrinsDJ;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  
  Matrix<double>  insVDinsDJ;//zeros(model.max_insertions +1, model.max_insertions +1);
  
  Matrix<double>  VdelV;//zeros(size(model.PV,1), model.max_V_deletions+1 );
  Matrix<double>  DdelDl;//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
  Matrix<double>  DdelDr;//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
  Matrix<double>  JdelJ;//zeros(size(model.PDJ,2), model.max_J_deletions+1 );
  
  Matrix<double>  VinsVD;//zeros(size(model.PV,1), model.max_insertions+1 );
  Matrix<double>  DinsVD;//zeros(size(model.PDJ,1), model.max_insertions+1 );
  Matrix<double>  DinsDJ;//zeros(size(model.PDJ,1), model.max_insertions+1 );
  Matrix<double>  JinsDJ;//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
  Matrix<double>  VdelDl;//zeros(size(model.PV,1), model.max_D_deletions+1 );
  Matrix<double>  VdelDr;//zeros(size(model.PV,1), model.max_D_deletions+1 );
  Matrix<double>  VdelJ;//zeros(size(model.PV,1), model.max_J_deletions+1 );
  
  Matrix<double>  JdelV;//zeros(size(model.PDJ,2), model.max_V_deletions+1 );
  Matrix<double>  JdelDl;//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
  Matrix<double>  JdelDr;//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
  
  Matrix<double>  DdelV;//zeros(size(model.PDJ,1), model.max_V_deletions+1 );
  Matrix<double>  DdelJ;//zeros(size(model.PDJ,1), model.max_J_deletions+1 );
  
  Matrix<double>  VinsDJ;//zeros(size(model.PV,1), model.max_insertions+1 );
  Matrix<double>  JinsVD;//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
  Matrix<double>  pVinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double>  pVinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double>  pVdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double>  pVdelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double>  pVdelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  Matrix<double>  VpV;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<double>  JpV;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<double>  DpV;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<double>  pJinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double>  pJinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double>  pJdelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  Matrix<double>  pJdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double>  pJdelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double>  VpJ;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<double>  JpJ;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<double>  DpJ;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<double>  pDlinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double>  pDlinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double>  pDldelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  Matrix<double>  pDldelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  Matrix<double>  pDldelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double>  VpDl;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<double>  JpDl;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<double>  DpDl;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<double>  pDrinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double>  pDrinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double>  pDrdelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  Matrix<double>  pDrdelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  Matrix<double>  pDrdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double>  VpDr;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<double>  JpDr;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<double>  DpDr;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<double>  pVpDl;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<double>  pVpDr;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<double>  pVpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<double>  pDlpDr;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<double>  pDlpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<double>  pDrpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  
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
