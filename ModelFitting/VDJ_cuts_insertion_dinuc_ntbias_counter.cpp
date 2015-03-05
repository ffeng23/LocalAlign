
#include "VDJ_cuts_insertion_dinuc_ntbias_counter.hpp"
#include "../SIGPIG/AlignmentSettings.hpp"
//NOTE: to do:
//model read length ??? what is the value, what is the purpose
//max_V/D/J_length in the alignmentSettings???

VDJ_cuts_insertion_dinuc_ntbias_counter::VDJ_cuts_insertion_dinuc_ntbias_counter(const 
VDJ_cuts_insertion_dinuc_ntbias_model& _model):
  /*initialization list, intialize default value and do matrices too*/
  logLikelihood(0),//This adds up the log likelihood of generating a sequence given the model
  N_processed(0),

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
  
  nMerror(0), nMsequenced_nucleotide(0) ,//error rate, will be 
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

  nMzeroD(0),// = 0,
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


  max_V_length(AlignmentSettings::max_V_length),//max v_region length = model.read_length,
  max_D_length(AlignmentSettings::max_D_length),// = 16, max d region length
  max_J_length(AlignmentSettings::max_J_length),// = 18, max j region length

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
  
  nMerror_vs_position(),//1,_model.read_length),//zeros(model.read_length,1),
  nMcoverage()//1, _model.read_length)//zeros(model.read_length,1);
{
  //here we will correctly initialize the matrices
  unsigned dim_size1[]={_model.max_insertions+1};
  nPinsVD.initialize(1, dim_size1, 0.0);
  nPinsDJ.initialize(1, dim_size1, 0.0);
  
  unsigned dim_size2[]={_model.number_V_genes, (unsigned)(_model.max_V_cut-_model.min_V_cut+1)};
  nPcutV_given_V.initialize(2, dim_size2, 0.0);
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=(unsigned)(_model.max_J_cut-_model.min_J_cut+1);
  nPcutJ_given_J.initialize(2,dim_size2, 0.0);
  
  unsigned dim_size3[3]={_model.number_D_genes, (unsigned)(_model.max_D_cut-_model.min_D_cut+1), (unsigned)(_model.max_D_cut-_model.min_D_cut+1) };
  nPcutDlcutDr_given_D.initialize(3, dim_size3, 0.0);
  
  unsigned dim_size[1]={_model.number_V_genes};
  nPV.initialize(1, dim_size, 0.0);
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.number_J_genes;
  nPDJ.initialize(2, dim_size2, 0.0);
  
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_V_n_alleles;
  nPVallele_given_gene.initialize(2, dim_size2, 0.0);
  
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_D_n_alleles;
  nPDallele_given_gene.initialize(2, dim_size2,0.0);// zeros(max_D_alleles , D_genes);
  
  dim_size2[0]=4;
  dim_size2[1]=4;
  nMnucleotideVD.initialize(2, dim_size2, 0.0); nMnucleotideVD_5prime.initialize(2, dim_size2,0.0);
  
  //RnucleotideVD_per_nucleotideVD_5prime=RnucleotideVD_per_nucleotideVD_5prime/4;
  nMnucleotideDJ.initialize(2, dim_size2, 0.0);nMnucleotideDJ_3prime.initialize(2, dim_size2,0.0);
  
  //===================>>>>>>>>>>>
  //start initializing the counter specified fields
  dim_size[0]=4;
  nMmononucleotideVD.initialize(1,dim_size,0.0);// zeros(4,1),
  nMmononucleotideDJ.initialize(1, dim_size, 0.0);// = zeros(4,1),
  nMinsertionVD.initialize(1, dim_size, 0.0);// = zeros(4,1),
  nMinsertionDJ.initialize(1, dim_size, 0.0);// = zeros(4,1),

  dim_size2[0]=4;
  dim_size2[1]=4;
  nPVD_left_edge_dinucleotide.initialize(2, dim_size2,0.0);// = zeros(4,4),
  nPVD_right_edge_dinucleotide.initialize(2, dim_size2,0.0);// = zeros(4,4),
  nPDJ_left_edge_dinucleotide.initialize(2, dim_size2,0.0);// = zeros(4,4),
  nPDJ_right_edge_dinucleotide.initialize(2, dim_size2,0.0);// = zeros(4,4),

  nMtrinucleotideVD.initialize(2, dim_size2,0.0);// = zeros(4,4,4),
  nMtrinucleotideDJ.initialize(2, dim_size2,0.0);// = zeros(4,4,4),
  
  dim_size3[0]=_model.number_V_genes;
  dim_size3[1]=_model.max_palindrome+1;
  dim_size3[2]=_model.max_V_deletions+1;
  nPpVmax_delV_V.initialize(3, dim_size3,0.0);// =  zeros(model.max_palindrome + 1, model.max_V_deletions + 1, size(model.PV,1)),
  dim_size3[0]=_model.number_J_genes;
  dim_size3[1]=_model.max_J_deletions+1;
  nPpJmax_delJ_J.initialize(3, dim_size3,0.0);// = zeros(model.max_palindrome + 1, model.max_J_deletions + 1, size(model.PDJ,2)),
  dim_size3[0]=_model.number_D_genes;
  dim_size3[1]=_model.max_D_deletions+1;
  nPpDlmax_delDl_D.initialize(3, dim_size3,0.0);// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1)),
  nPpDrmax_delDr_D.initialize(3, dim_size3,0.0);// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1)),

  //nMzeroD(0),// = 0,
  dim_size3[0]=_model.number_V_genes;
  dim_size3[1]=_model.number_D_genes;
  dim_size3[2]=_model.number_J_genes;
  nPVDJ.initialize(3, dim_size3,0.0);// = zeros(size(model.PV,1), size(model.PDJ,1), size(model.PDJ,2)),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_V_deletions+1;
  nPpVdelV.initialize(2, dim_size2,0.0);// = zeros( model.max_palindrome + 1, model.max_V_deletions + 1),
  dim_size2[1]=_model.max_D_deletions+1;
  nPpDldelDl.initialize(2, dim_size2,0.0);// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1),
  nPpDrdelDr.initialize(2, dim_size2,0.0);// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1),
  dim_size2[1]=_model.max_J_deletions+1;
  nPpJdelJ.initialize(2, dim_size2,0.0);// = zeros( model.max_palindrome + 1, model.max_J_deletions + 1),
  
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_V_deletions+1+_model.max_palindrome;
  nPVcutV.initialize(2, dim_size2,0.0);// = zeros(size(model.PV,1),model.max_palindrome + model.max_V_deletions + 1),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_D_deletions+1+_model.max_palindrome;
  nPDcutDl.initialize(2, dim_size2,0.0);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_D_deletions+1+_model.max_palindrome;
  nPDcutDr.initialize(2, dim_size2,0.0);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1),
  
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_J_deletions+1+_model.max_palindrome;
  nPJcutJ.initialize(2, dim_size2,0.0);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_J_deletions + 1),
  
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_V_deletions+1+_model.max_palindrome;
  nPDcutV.initialize(2, dim_size2,0.0);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_V_deletions + 1),
  
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_J_deletions+1+_model.max_palindrome;
  nPDcutJ.initialize(2, dim_size2,0.0);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_J_deletions + 1),
  
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_J_deletions+1+_model.max_palindrome;
  nPVcutJ.initialize(2, dim_size2,0.0);// = zeros(size(model.PV,1),model.max_palindrome + model.max_J_deletions + 1),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_D_deletions+1+_model.max_palindrome;
  nPVcutDl.initialize(2, dim_size2,0.0);// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_D_deletions+1+_model.max_palindrome;
  nPVcutDr.initialize(2, dim_size2,0.0);// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1),
  
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_D_deletions+1+_model.max_palindrome;
  nPJcutDl.initialize(2, dim_size2,0.0);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_D_deletions+1+_model.max_palindrome;
  nPJcutDr.initialize(2, dim_size2,0.0);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_V_deletions+1+_model.max_palindrome;
  nPJcutV.initialize(2, dim_size2,0.0);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_V_deletions + 1),
  
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=_model.max_V_deletions+1+_model.max_palindrome;
  nPinsVDcutV.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1),
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=_model.max_V_deletions+1+_model.max_palindrome;
  nPinsDJcutV.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1),
  
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=_model.max_D_deletions+1+_model.max_palindrome;
  nPinsVDcutDl.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1),
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=_model.max_D_deletions+1+_model.max_palindrome;
  nPinsDJcutDl.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1),
  
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=_model.max_D_deletions+1+_model.max_palindrome;
  nPinsVDcutDr.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1),
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=_model.max_D_deletions+1+_model.max_palindrome;
  nPinsDJcutDr.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1),
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=_model.max_J_deletions+1+_model.max_palindrome;
  nPinsVDcutJ.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1),
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=_model.max_J_deletions+1+_model.max_palindrome;
  nPinsDJcutJ.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1),

  //???????????????????????????/to do
  dim_size[0]=max_V_length+1;
  nPV_align_length.initialize(1, dim_size,0.0);// = zeros(max_V_length + 1,1),
  dim_size[0]=max_D_length+1;
  nPD_align_length.initialize(1, dim_size,0.0);// = zeros(max_D_length + 1,1),
  dim_size[0]=max_J_length+1;
  nPJ_align_length.initialize(1, dim_size,0.0);// = zeros(max_J_length + 1,1),

  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=max_V_length+_model.max_V_deletions;
  nMVV_err_pos.initialize(2, dim_size2,0.0);// = zeros(size(model.PV,1), max_V_length + model.max_V_deletions),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=max_J_length;
  nMJJ_err_pos.initialize(2, dim_size2,0.0);// = zeros(size(model.PDJ,2), max_J_length),
  //Note:: WHY this is different from above V case, not adding max_V_deletions

  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=max_J_length+1;
  nPJJ_align_length.initialize(2, dim_size2,0.0);// = zeros(size(model.PDJ,2), 1 + max_J_length),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=max_V_length+1;  
  nPVV_align_length.initialize(2, dim_size2,0.0);// = zeros(size(model.PV,1), 1 + max_V_length),
  
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=max_D_length+1;  
  nPinsDJ_D_align_length.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions +1, max_D_length + 1),
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=max_D_length+1;  
  nPinsVD_D_align_length.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions +1, max_D_length + 1),

  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=max_J_length+1;  
  nPinsDJ_J_align_length.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions +1, max_J_length + 1),
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=max_J_length+1;  
  nPinsVD_J_align_length.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions +1, max_J_length + 1),

  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=max_V_length+1;  
  nPinsDJ_V_align_length.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions +1, max_V_length + 1),
  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=max_V_length+1;
  nPinsVD_V_align_length.initialize(2, dim_size2,0.0);// = zeros(model.max_insertions +1, max_V_length + 1),

  dim_size2[0]=_model.max_D_n_alleles;
  dim_size2[1]=max_D_length+1;
  nPDallele_D_align_length.initialize(2, dim_size2,0.0);// = zeros(3, max_D_length + 1),

  dim_size2[0]=_model.max_V_deletions+1;
  dim_size2[1]=_model.max_insertions+1;
  nPdelVinsVD.initialize(2, dim_size2,0.0);// = zeros(model.max_V_deletions +1, model.max_insertions +1),
  dim_size2[0]=_model.max_V_deletions+1;
  dim_size2[1]=_model.max_insertions+1;
  nPdelVinsDJ.initialize(2, dim_size2,0.0);// = zeros(model.max_V_deletions +1, model.max_insertions +1),
  dim_size2[0]=_model.max_V_deletions+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPdelVdelDl.initialize(2, dim_size2,0.0);// = zeros(model.max_V_deletions +1, model.max_D_deletions +1),
  
  dim_size2[0]=_model.max_V_deletions+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPdelVdelDr.initialize(2, dim_size2,0.0);// = zeros(model.max_V_deletions +1, model.max_D_deletions +1),
  
  dim_size2[0]=_model.max_V_deletions+1;
  dim_size2[1]=_model.max_J_deletions+1;
  nPdelVdelJ.initialize(2, dim_size2,0.0);// = zeros(model.max_V_deletions +1, model.max_J_deletions +1),
  
  dim_size2[0]=_model.max_J_deletions+1;
  dim_size2[1]=_model.max_insertions+1;
  nPdelJinsVD.initialize(2, dim_size2,0.0);// = zeros(model.max_J_deletions +1, model.max_insertions +1),
  dim_size2[0]=_model.max_J_deletions+1;
  dim_size2[1]=_model.max_insertions+1;
  nPdelJinsDJ.initialize(2, dim_size2,0.0);//zeros(model.max_J_deletions +1, model.max_insertions +1),
  dim_size2[0]=_model.max_J_deletions+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPdelJdelDl.initialize(2, dim_size2,0.0);//zeros(model.max_J_deletions +1, model.max_D_deletions +1),
  dim_size2[0]=_model.max_J_deletions+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPdelJdelDr.initialize(2, dim_size2,0.0);//zeros(model.max_J_deletions +1, model.max_D_deletions +1),
  dim_size2[0]=_model.max_D_deletions+1;
  dim_size2[1]=_model.max_insertions+1;
  nPdelDlinsVD.initialize(2, dim_size2,0.0);//zeros(model.max_D_deletions +1, model.max_insertions +1),
  dim_size2[0]=_model.max_J_deletions+1;
  dim_size2[1]=_model.max_insertions+1;
  nPdelDlinsDJ.initialize(2, dim_size2,0.0);//zeros(model.max_D_deletions +1, model.max_insertions +1),
  dim_size2[0]=_model.max_D_deletions+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPdelDldelDr.initialize(2, dim_size2,0.0);//zeros(model.max_D_deletions +1, model.max_D_deletions +1),

  dim_size2[0]=_model.max_D_deletions+1;
  dim_size2[1]=_model.max_insertions+1;
  nPdelDrinsVD.initialize(2, dim_size2,0.0);//zeros(model.max_D_deletions +1, model.max_insertions +1),
  dim_size2[0]=_model.max_D_deletions+1;
  dim_size2[1]=_model.max_insertions+1;
  nPdelDrinsDJ.initialize(2, dim_size2,0.0);//zeros(model.max_D_deletions +1, model.max_insertions +1),

  dim_size2[0]=_model.max_insertions+1;
  dim_size2[1]=_model.max_insertions+1;
  nPinsVDinsDJ.initialize(2, dim_size2,0.0);//zeros(model.max_insertions +1, model.max_insertions +1),

  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_V_deletions+1;
  nPVdelV.initialize(2, dim_size2,0.0);//zeros(size(model.PV,1), model.max_V_deletions+1 ),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_D_deletions+1;
  nPDdelDl.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,1), model.max_D_deletions+1 ),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_D_deletions+1;
  nPDdelDr.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,1), model.max_D_deletions+1 ),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_J_deletions+1;
  nPJdelJ.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,2), model.max_J_deletions+1 ),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_insertions+1;
  nPVinsVD.initialize(2, dim_size2,0.0);//zeros(size(model.PV,1), model.max_insertions+1 ),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_insertions+1;
  nPDinsVD.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,1), model.max_insertions+1 ),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_insertions+1;
  nPDinsDJ.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,1), model.max_insertions+1 ),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_insertions+1;
  nPJinsDJ.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,2), model.max_insertions+1 ),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_D_deletions+1;
  nPVdelDl.initialize(2, dim_size2,0.0);//zeros(size(model.PV,1), model.max_D_deletions+1 ),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_D_deletions+1;
  nPVdelDr.initialize(2, dim_size2,0.0);//zeros(size(model.PV,1), model.max_D_deletions+1 ),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_J_deletions+1;
  nPVdelJ.initialize(2, dim_size2,0.0);//zeros(size(model.PV,1), model.max_J_deletions+1 ),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_V_deletions+1;
  nPJdelV.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,2), model.max_V_deletions+1 ),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_D_deletions+1;
  nPJdelDl.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,2), model.max_D_deletions+1 ),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_D_deletions+1;
  nPJdelDr.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,2), model.max_D_deletions+1 ),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_V_deletions+1;
  nPDdelV.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,1), model.max_V_deletions+1 ),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_J_deletions+1;
  nPDdelJ.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,1), model.max_J_deletions+1 ),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_insertions+1;
  nPVinsDJ.initialize(2, dim_size2,0.0);//zeros(size(model.PV,1), model.max_insertions+1 ),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_insertions+1;
  nPJinsVD.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,2), model.max_insertions+1 ),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_insertions+1;
  nPpVinsVD.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_insertions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_insertions+1;
  nPpVinsDJ.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_insertions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPpVdelDl.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_D_deletions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPpVdelDr.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_D_deletions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_J_deletions+1;
  nPpVdelJ.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_J_deletions +1),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPVpV.initialize(2, dim_size2,0.0);//zeros(size(model.PV,1), model.max_palindrome+1 ),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPJpV.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,2), model.max_palindrome+1 ),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPDpV.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,1), model.max_palindrome+1 ),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_insertions+1;
  nPpJinsVD.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_insertions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_insertions+1;
  nPpJinsDJ.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_insertions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_V_deletions+1;
  nPpJdelV.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_V_deletions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPpJdelDl.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_D_deletions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPpJdelDr.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_D_deletions +1),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPVpJ.initialize(2, dim_size2,0.0);//zeros(size(model.PV,1), model.max_palindrome+1 ),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPJpJ.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,2), model.max_palindrome+1 ),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPDpJ.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,1), model.max_palindrome+1 ),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_insertions+1;
  nPpDlinsVD.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_insertions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_insertions+1;
  nPpDlinsDJ.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_insertions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_V_deletions+1;
  nPpDldelV.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_V_deletions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_J_deletions+1;
  nPpDldelJ.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_J_deletions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPpDldelDr.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_D_deletions +1),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPVpDl.initialize(2, dim_size2,0.0);//zeros(size(model.PV,1), model.max_palindrome+1 ),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPJpDl.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,2), model.max_palindrome+1 ),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPDpDl.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,1), model.max_palindrome+1 ),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_insertions+1;
  nPpDrinsVD.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_insertions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_insertions+1;
  nPpDrinsDJ.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_insertions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_V_deletions+1;
  nPpDrdelV.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_V_deletions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_J_deletions+1;
  nPpDrdelJ.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_J_deletions +1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_D_deletions+1;
  nPpDrdelDl.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome +1, model.max_D_deletions +1),
  dim_size2[0]=_model.number_V_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPVpDr.initialize(2, dim_size2,0.0);//zeros(size(model.PV,1), model.max_palindrome+1 ),
  dim_size2[0]=_model.number_J_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPJpDr.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,2), model.max_palindrome+1 ),
  dim_size2[0]=_model.number_D_genes;
  dim_size2[1]=_model.max_palindrome+1;
  nPDpDr.initialize(2, dim_size2,0.0);//zeros(size(model.PDJ,1), model.max_palindrome+1 ),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_palindrome+1;
  nPpVpDl.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_palindrome+1;
  nPpVpDr.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_palindrome+1;
  nPpVpJ.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_palindrome+1;
  nPpDlpDr.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_palindrome+1;
  nPpDlpJ.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  dim_size2[0]=_model.max_palindrome+1;
  dim_size2[1]=_model.max_palindrome+1;
  nPpDrpJ.initialize(2, dim_size2,0.0);//zeros(model.max_palindrome + 1, model.max_palindrome + 1),
  dim_size[0]=_model.read_length;
  nMerror_vs_position.initialize(1, dim_size,0.0);//1,_model.read_length),//zeros(model.read_length,1),
  
  nMcoverage.initialize(1, dim_size,0.0);//1, _model.read_length)//zeros(model.read_length,1);

  //finally done!!

}//end of counter constructor

VDJ_cuts_insertion_dinuc_ntbias_counter::~VDJ_cuts_insertion_dinuc_ntbias_counter()
{
  //empty so far
}
