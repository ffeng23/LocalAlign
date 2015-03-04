#ifndef VDJ_CUTS_INSERTION_DINUC_NTBIAS_COUNTER_HPP
#define VDJ_CUTS_INSERTION_DINUC_NTBIAS_COUNTER_HPP

#include "Counter.hpp"
#include "VDJ_cuts_insertion_dinuc_ntbias_model.hpp"

class VDJ_cuts_insertion_dinuc_ntbias_counter:public Counter
{
  VDJ_cuts_insertion_dinuc_ntbias_counter(VDJ_cuts_insertion_dinuc_ntbias_model, _model);

  virtual ~VDJ_cuts_insertion_dinuc_ntbias_counter();

  //==========================>member declaration, from the model
  double logLikelihood;//This adds up the log likelihood of generating a sequence given the model
  unsigned N_processed;

  //% For all parameters of model that begin with a 'P' or an 'M' AND don't have the word 'specific' in them,
  //% make a corresponding variable in counter, with an 'n' prefixed to it.
  Matrix<double> nPinsVD;
  Matrix<double> nPinsDJ;

  Matrix<double> nPcutV_given_V;
  Matrix<double> nPcutJ_given_J;
  Matrix<double> nPcutDlcutDr_given_D;

  Matrix<double> nPV;
  Matrix<double> nPDJ; //Joint P(V, D, J gene choices)
  Matrix<double> nPVallele_given_gene; //Probabilities of alleles given gene for each gene
  Matrix<double> nPDallele_given_gene;

  //_per_
  //% Rate parameter. So make two 'nM' variables in counter, one for numerator and one for denominator.
  //% Names of these two should be separated by '_per_'.
  //% The ratio of these will be estimate of the model parameter.
  Matrix<double> nMnucleotideVD; Matrix<double> nMnucleotideVD_5prime;//nucleotide distr's
  Matrix<double> nMnucleotideDJ; Matrix<double> nMnucleotideDJ_3prime;

  double nMerror; double nMsequenced_nucleotide ;//error rate, will be 
  //============the above members are from model

  //==========by the counter
  //add any other fields you might want to keep track of, but are not in model,
  //below, begin the name with an 'nP' or 'nM',
  //if it is a distribution or mean value, respectively;
  Matrix<double> nMmononucleotideVD;// zeros(4,1);
  Matrix<double> nMmononucleotideDJ;// = zeros(4,1);
  Matrix<double> nMinsertionVD;// = zeros(4,1);
  Matrix<double> nMinsertionDJ;// = zeros(4,1);

  Matrix<double> nPVD_left_edge_dinucleotide;// = zeros(4,4);
  Matrix<double> nPVD_right_edge_dinucleotide;// = zeros(4,4);
  Matrix<double> nPDJ_left_edge_dinucleotide;// = zeros(4,4);
  Matrix<double> nPDJ_right_edge_dinucleotide;// = zeros(4,4);

  Matrix<double> nMtrinucleotideVD;// = zeros(4,4,4);
  Matrix<double> nMtrinucleotideDJ;// = zeros(4,4,4);
  
  Matrix<double> nPpVmax_delV_V;// =  zeros(model.max_palindrome + 1, model.max_V_deletions + 1, size(model.PV,1));
  Matrix<double> nPpJmax_delJ_J;// = zeros(model.max_palindrome + 1, model.max_J_deletions + 1, size(model.PDJ,2));

  Matrix<double> nPpDlmax_delDl_D;// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));
  Matrix<double> nPpDrmax_delDr_D;// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));

  double nMzeroD = 0;
  Matrix<double> nPVDJ;// = zeros(size(model.PV,1), size(model.PDJ,1), size(model.PDJ,2));
  Matrix<double> nPpVdelV;// = zeros( model.max_palindrome + 1, model.max_V_deletions + 1);
  Matrix<double> nPpDldelDl;// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
  Matrix<double> nPpDrdelDr;// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
  Matrix<double> nPpJdelJ;// = zeros( model.max_palindrome + 1, model.max_J_deletions + 1);


  Matrix<double> nPVcutV;// = zeros(size(model.PV,1),model.max_palindrome + model.max_V_deletions + 1);
  Matrix<double> nPDcutDl;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double> nPDcutDr;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double> nPJcutJ;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_J_deletions + 1);

  Matrix<double> nPDcutV;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_V_deletions + 1);
  Matrix<double> nPDcutJ;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_J_deletions + 1);

  Matrix<double> nPVcutJ;// = zeros(size(model.PV,1),model.max_palindrome + model.max_J_deletions + 1);
  Matrix<double> nPVcutDl;// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double> nPVcutDr;// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);

  Matrix<double> nPJcutDl;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double> nPJcutDr;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double> nPJcutV;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_V_deletions + 1);

  Matrix<double> nPinsVDcutV;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);
  Matrix<double> nPinsDJcutV;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);

  Matrix<double> nPinsVDcutDl;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double> nPinsDJcutDl;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

  Matrix<double> nPinsVDcutDr;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  Matrix<double> nPinsDJcutDr;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

  Matrix<double> nPinsVDcutJ;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);
  Matrix<double> nPinsDJcutJ;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);


  unsigned max_V_length;// = model.read_length;
  unsigned max_D_length;// = 16;
  unsigned max_J_length;// = 18;

  Matrix<double> nPV_align_length;// = zeros(max_V_length + 1,1);
  Matrix<double> nPD_align_length;// = zeros(max_D_length + 1,1);
  Matrix<double> nPJ_align_length;// = zeros(max_J_length + 1,1);

  Matrix<double> nMVV_err_pos;// = zeros(size(model.PV,1), max_V_length + model.max_V_deletions);
  Matrix<double> nMJJ_err_pos;// = zeros(size(model.PDJ,2), max_J_length);

  Matrix<double> nPJJ_align_length;// = zeros(size(model.PDJ,2), 1 + max_J_length);
  Matrix<double> nPVV_align_length;// = zeros(size(model.PV,1), 1 + max_V_length);


  Matrix<double> nPinsDJ_D_align_length;// = zeros(model.max_insertions +1, max_D_length + 1);
  Matrix<double> nPinsVD_D_align_length;// = zeros(model.max_insertions +1, max_D_length + 1);

  Matrix<double> nPinsDJ_J_align_length;// = zeros(model.max_insertions +1, max_J_length + 1);
  Matrix<double> nPinsVD_J_align_length;// = zeros(model.max_insertions +1, max_J_length + 1);

  Matrix<double> nPinsDJ_V_align_length;// = zeros(model.max_insertions +1, max_V_length + 1);
  Matrix<double> nPinsVD_V_align_length;// = zeros(model.max_insertions +1, max_V_length + 1);

  Matrix<double> nPDallele_D_align_length;// = zeros(3, max_D_length + 1);


  Matrix<double> nPdelVinsVD;// = zeros(model.max_V_deletions +1, model.max_insertions +1);
  Matrix<double> nPdelVinsDJ;// = zeros(model.max_V_deletions +1, model.max_insertions +1);
  Matrix<double> nPdelVdelDl;// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
  Matrix<double> nPdelVdelDr;// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
  Matrix<double> nPdelVdelJ;// = zeros(model.max_V_deletions +1, model.max_J_deletions +1);

  Matrix<double> nPdelJinsVD;// = zeros(model.max_J_deletions +1, model.max_insertions +1);
  Matrix<double> nPdelJinsDJ;//zeros(model.max_J_deletions +1, model.max_insertions +1);
  Matrix<double> nPdelJdelDl;//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
  Matrix<double> nPdelJdelDr;//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
  
  Matrix<double> nPdelDlinsVD;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  Matrix<double> nPdelDlinsDJ;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  Matrix<double> nPdelDldelDr;//zeros(model.max_D_deletions +1, model.max_D_deletions +1);
  
  Matrix<double> nPdelDrinsVD;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  Matrix<double> nPdelDrinsDJ;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  
  Matrix<double> nPinsVDinsDJ;//zeros(model.max_insertions +1, model.max_insertions +1);
  
  Matrix<double> nPVdelV;//zeros(size(model.PV,1), model.max_V_deletions+1 );
  Matrix<double> nPDdelDl;//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
  Matrix<double> nPDdelDr;//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
  Matrix<double> nPJdelJ;//zeros(size(model.PDJ,2), model.max_J_deletions+1 );
  
  Matrix<double> nPVinsVD;//zeros(size(model.PV,1), model.max_insertions+1 );
  Matrix<double> nPDinsVD;//zeros(size(model.PDJ,1), model.max_insertions+1 );
  Matrix<double> nPDinsDJ;//zeros(size(model.PDJ,1), model.max_insertions+1 );
  Matrix<double> nPJinsDJ;//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
  Matrix<double> nPVdelDl;//zeros(size(model.PV,1), model.max_D_deletions+1 );
  Matrix<double> nPVdelDr;//zeros(size(model.PV,1), model.max_D_deletions+1 );
  Matrix<double> nPVdelJ;//zeros(size(model.PV,1), model.max_J_deletions+1 );
  
  Matrix<double> nPJdelV;//zeros(size(model.PDJ,2), model.max_V_deletions+1 );
  Matrix<double> nPJdelDl;//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
  Matrix<double> nPJdelDr;//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
  
  Matrix<double> nPDdelV;//zeros(size(model.PDJ,1), model.max_V_deletions+1 );
  Matrix<double> nPDdelJ;//zeros(size(model.PDJ,1), model.max_J_deletions+1 );
  
  Matrix<double> nPVinsDJ;//zeros(size(model.PV,1), model.max_insertions+1 );
  Matrix<double> nPJinsVD;//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
  Matrix<double> nPpVinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double> nPpVinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double> nPpVdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double> nPpVdelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double> nPpVdelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  Matrix<double> nPVpV;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<double> nPJpV;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<double> nPDpV;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<double> nPpJinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double> nPpJinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double> nPpJdelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  Matrix<double> nPpJdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double> nPpJdelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double> nPVpJ;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<double> nPJpJ;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<double> nPDpJ;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<double> nPpDlinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double> nPpDlinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double> nPpDldelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  Matrix<double> nPpDldelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  Matrix<double> nPpDldelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double> nPVpDl;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<double> nPJpDl;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<double> nPDpDl;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<double> nPpDrinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double> nPpDrinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  Matrix<double> nPpDrdelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  Matrix<double> nPpDrdelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  Matrix<double> nPpDrdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  Matrix<double> nPVpDr;//zeros(size(model.PV,1), model.max_palindrome+1 );
  Matrix<double> nPJpDr;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  Matrix<double> nPDpDr;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  Matrix<double> nPpVpDl;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<double> nPpVpDr;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<double> nPpVpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<double> nPpDlpDr;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<double> nPpDlpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  Matrix<double> nPpDrpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  
  //%
  Matrix<double> nMerror_vs_position;//zeros(model.read_length,1);
  Matrix<double> nMcoverage;//zeros(model.read_length,1);


};

#endif
