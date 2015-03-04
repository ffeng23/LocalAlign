
#include "VDJ_cuts_insertion_dinuc_ntbias_counter.hpp"

VDJ_cuts_insertion_dinuc_ntbias_counter::VDJ_cuts_insertion_dinuc_ntbias_counter():
  /*initialization list, intialize default value and do matrices too*/
  logLikelihood,//This adds up the log likelihood of generating a sequence given the model
  unsigned N_processed,

  //% For all parameters of model that begin with a 'P' or an 'M' AND don't have the word 'specific' in them,
  //% make a corresponding variable in counter, with an 'n' prefixed to it.
   nPinsVD(),
   nPinsDJ(),

   nPcutV_given_V(),
   nPcutJ_given_J(),
   nPcutDlcutDr_given_D(),

   nPV(),
   nPDJ(), //Joint P(V, D, J gene choices)
   nPVallele_given_gene(), //Probabilities of alleles given gene for each gene
   nPDallele_given_gene(),

  //_per_
  //% Rate parameter. So make two 'nM' variables in counter, one for numerator and one for denominator.
  //% Names of these two should be separated by '_per_'.
  //% The ratio of these will be estimate of the model parameter.
   nMnucleotideVD(),  nMnucleotideVD_5prime(),//nucleotide distr's
   nMnucleotideDJ(),  nMnucleotideDJ_3prime(),

  double nMerror, double nMsequenced_nucleotide ,//error rate, will be 
  //============the above members are from model

  //==========by the counter
  //add any other fields you might want to keep track of, but are not in model,
  //below, begin the name with an 'nP' or 'nM',
  //if it is a distribution or mean value, respectively,
   nMmononucleotideVD(),// zeros(4,1),
   nMmononucleotideDJ(),// = zeros(4,1),
   nMinsertionVD(),// = zeros(4,1),
   nMinsertionDJ(),// = zeros(4,1),

   nPVD_left_edge_dinucleotide(),// = zeros(4,4),
   nPVD_right_edge_dinucleotide(),// = zeros(4,4),
  nPDJ_left_edge_dinucleotide(),// = zeros(4,4),
  nPDJ_right_edge_dinucleotide(),// = zeros(4,4),

  nMtrinucleotideVD(),// = zeros(4,4,4),
  nMtrinucleotideDJ(),// = zeros(4,4,4),
  
  nPpVmax_delV_V(),// =  zeros(model.max_palindrome + 1, model.max_V_deletions + 1, size(model.PV,1)),
  nPpJmax_delJ_J(),// = zeros(model.max_palindrome + 1, model.max_J_deletions + 1, size(model.PDJ,2)),

  nPpDlmax_delDl_D(),// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1)),
  nPpDrmax_delDr_D(),// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1)),

  double nMzeroD = 0,
  nPVDJ(),// = zeros(size(model.PV,1), size(model.PDJ,1), size(model.PDJ,2)),
  nPpVdelV(),// = zeros( model.max_palindrome + 1, model.max_V_deletions + 1),
   nPpDldelDl(),// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1),
   nPpDrdelDr(),// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1),
   nPpJdelJ(),// = zeros( model.max_palindrome + 1, model.max_J_deletions + 1),


   nPVcutV(),// = zeros(size(model.PV,1),model.max_palindrome + model.max_V_deletions + 1),
   nPDcutDl(),// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1),
   nPDcutDr(),// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1),
   nPJcutJ(),// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_J_deletions + 1),

   nPDcutV(),// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_V_deletions + 1),
   nPDcutJ(),// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_J_deletions + 1),

   nPVcutJ(),// = zeros(size(model.PV,1),model.max_palindrome + model.max_J_deletions + 1),
   nPVcutDl(),// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1),
   nPVcutDr(),// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1),

   nPJcutDl(),// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1),
   nPJcutDr(),// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1),
   nPJcutV(),// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_V_deletions + 1),

   nPinsVDcutV(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1),
   nPinsDJcutV(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1),

   nPinsVDcutDl(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1),
   nPinsDJcutDl(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1),

   nPinsVDcutDr(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1),
   nPinsDJcutDr(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1),

   nPinsVDcutJ(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1),
   nPinsDJcutJ(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1),


  unsigned max_V_length,// = model.read_length,
  unsigned max_D_length,// = 16,
  unsigned max_J_length,// = 18,

   nPV_align_length(),// = zeros(max_V_length + 1,1),
   nPD_align_length(),// = zeros(max_D_length + 1,1),
   nPJ_align_length(),// = zeros(max_J_length + 1,1),

   nMVV_err_pos(),// = zeros(size(model.PV,1), max_V_length + model.max_V_deletions),
   nMJJ_err_pos(),// = zeros(size(model.PDJ,2), max_J_length),

   nPJJ_align_length(),// = zeros(size(model.PDJ,2), 1 + max_J_length),
   nPVV_align_length(),// = zeros(size(model.PV,1), 1 + max_V_length),


   nPinsDJ_D_align_length(),// = zeros(model.max_insertions +1, max_D_length + 1),
   nPinsVD_D_align_length(),// = zeros(model.max_insertions +1, max_D_length + 1),

   nPinsDJ_J_align_length(),// = zeros(model.max_insertions +1, max_J_length + 1),
   nPinsVD_J_align_length(),// = zeros(model.max_insertions +1, max_J_length + 1),

   nPinsDJ_V_align_length(),// = zeros(model.max_insertions +1, max_V_length + 1),
   nPinsVD_V_align_length(),// = zeros(model.max_insertions +1, max_V_length + 1),

   nPDallele_D_align_length(),// = zeros(3, max_D_length + 1),


   nPdelVinsVD(),// = zeros(model.max_V_deletions +1, model.max_insertions +1),
   nPdelVinsDJ(),// = zeros(model.max_V_deletions +1, model.max_insertions +1),
   nPdelVdelDl(),// = zeros(model.max_V_deletions +1, model.max_D_deletions +1),
   nPdelVdelDr(),// = zeros(model.max_V_deletions +1, model.max_D_deletions +1),
   nPdelVdelJ(),// = zeros(model.max_V_deletions +1, model.max_J_deletions +1),

   nPdelJinsVD(),// = zeros(model.max_J_deletions +1, model.max_insertions +1),
   nPdelJinsDJ(),//zeros(model.max_J_deletions +1, model.max_insertions +1),
   nPdelJdelDl(),//zeros(model.max_J_deletions +1, model.max_D_deletions +1),
   nPdelJdelDr(),//zeros(model.max_J_deletions +1, model.max_D_deletions +1),
  
   nPdelDlinsVD(),//zeros(model.max_D_deletions +1, model.max_insertions +1),
   nPdelDlinsDJ(),//zeros(model.max_D_deletions +1, model.max_insertions +1),
   nPdelDldelDr(),//zeros(model.max_D_deletions +1, model.max_D_deletions +1),
  
   nPdelDrinsVD(),//zeros(model.max_D_deletions +1, model.max_insertions +1),
   nPdelDrinsDJ(),//zeros(model.max_D_deletions +1, model.max_insertions +1),
  
   nPinsVDinsDJ(),//zeros(model.max_insertions +1, model.max_insertions +1),
  
   nPVdelV(),//zeros(size(model.PV,1), model.max_V_deletions+1 ),
   nPDdelDl(),//zeros(size(model.PDJ,1), model.max_D_deletions+1 ),
   nPDdelDr(),//zeros(size(model.PDJ,1), model.max_D_deletions+1 ),
   nPJdelJ(),//zeros(size(model.PDJ,2), model.max_J_deletions+1 ),
  
   nPVinsVD(),//zeros(size(model.PV,1), model.max_insertions+1 ),
   nPDinsVD(),//zeros(size(model.PDJ,1), model.max_insertions+1 ),
   nPDinsDJ(),//zeros(size(model.PDJ,1), model.max_insertions+1 ),
   nPJinsDJ(),//zeros(size(model.PDJ,2), model.max_insertions+1 ),
  
   nPVdelDl(),//zeros(size(model.PV,1), model.max_D_deletions+1 ),
   nPVdelDr(),//zeros(size(model.PV,1), model.max_D_deletions+1 ),
   nPVdelJ(),//zeros(size(model.PV,1), model.max_J_deletions+1 ),
  
   nPJdelV(),//zeros(size(model.PDJ,2), model.max_V_deletions+1 ),
   nPJdelDl(),//zeros(size(model.PDJ,2), model.max_D_deletions+1 ),
   nPJdelDr(),//zeros(size(model.PDJ,2), model.max_D_deletions+1 ),
  
   nPDdelV(),//zeros(size(model.PDJ,1), model.max_V_deletions+1 ),
   nPDdelJ(),//zeros(size(model.PDJ,1), model.max_J_deletions+1 ),
  
   nPVinsDJ(),//zeros(size(model.PV,1), model.max_insertions+1 ),
   nPJinsVD(),//zeros(size(model.PDJ,2), model.max_insertions+1 ),
  
   nPpVinsVD(),//zeros(model.max_palindrome +1, model.max_insertions +1),
   nPpVinsDJ(),//zeros(model.max_palindrome +1, model.max_insertions +1),
   nPpVdelDl(),//zeros(model.max_palindrome +1, model.max_D_deletions +1),
   nPpVdelDr(),//zeros(model.max_palindrome +1, model.max_D_deletions +1),
   nPpVdelJ(),//zeros(model.max_palindrome +1, model.max_J_deletions +1),
   nPVpV(),//zeros(size(model.PV,1), model.max_palindrome+1 ),
   nPJpV(),//zeros(size(model.PDJ,2), model.max_palindrome+1 ),
   nPDpV(),//zeros(size(model.PDJ,1), model.max_palindrome+1 ),
  
   nPpJinsVD(),//zeros(model.max_palindrome +1, model.max_insertions +1),
   nPpJinsDJ(),//zeros(model.max_palindrome +1, model.max_insertions +1),
   nPpJdelV(),//zeros(model.max_palindrome +1, model.max_V_deletions +1),
   nPpJdelDl(),//zeros(model.max_palindrome +1, model.max_D_deletions +1),
   nPpJdelDr(),//zeros(model.max_palindrome +1, model.max_D_deletions +1),
   nPVpJ(),//zeros(size(model.PV,1), model.max_palindrome+1 ),
   nPJpJ(),//zeros(size(model.PDJ,2), model.max_palindrome+1 ),
   nPDpJ(),//zeros(size(model.PDJ,1), model.max_palindrome+1 ),
  
   nPpDlinsVD(),//zeros(model.max_palindrome +1, model.max_insertions +1),
   nPpDlinsDJ(),//zeros(model.max_palindrome +1, model.max_insertions +1),
   nPpDldelV(),//zeros(model.max_palindrome +1, model.max_V_deletions +1),
   nPpDldelJ(),//zeros(model.max_palindrome +1, model.max_J_deletions +1),
   nPpDldelDr(),//zeros(model.max_palindrome +1, model.max_D_deletions +1),
   nPVpDl(),//zeros(size(model.PV,1), model.max_palindrome+1 ),
   nPJpDl(),//zeros(size(model.PDJ,2), model.max_palindrome+1 ),
   nPDpDl(),//zeros(size(model.PDJ,1), model.max_palindrome+1 ),
  
   nPpDrinsVD(),//zeros(model.max_palindrome +1, model.max_insertions +1),
   nPpDrinsDJ(),//zeros(model.max_palindrome +1, model.max_insertions +1),
   nPpDrdelV(),//zeros(model.max_palindrome +1, model.max_V_deletions +1),
   nPpDrdelJ(),//zeros(model.max_palindrome +1, model.max_J_deletions +1),
   nPpDrdelDl(),//zeros(model.max_palindrome +1, model.max_D_deletions +1),
   nPVpDr(),//zeros(size(model.PV,1), model.max_palindrome+1 ),
   nPJpDr(),//zeros(size(model.PDJ,2), model.max_palindrome+1 ),
   nPDpDr(),//zeros(size(model.PDJ,1), model.max_palindrome+1 ),
  
   nPpVpDl(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
   nPpVpDr(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  nPpVpJ(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  nPpDlpDr(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  nPpDlpJ(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  nPpDrpJ(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  
  nMerror_vs_position(1,_model.read_length),//zeros(model.read_length,1),
  nMcoverage(1, _model.read_length)//zeros(model.read_length,1);
{
    
}

VDJ_cuts_insertion_dinuc_ntbias_counter::~VDJ_cuts_insertion_dinuc_ntbias_counter()
{
  //empty so far
}
