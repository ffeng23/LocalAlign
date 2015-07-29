#include "VDJ_cuts_insertion_dinuc_ntbias_assigns.hpp"


VDJ_cuts_insertion_dinuc_ntbias_assigns::VDJ_cuts_insertion_dinuc_ntbias_assigns(const VDJ_cuts_insertion_dinuc_ntbias_model& _model, const VDJ_cuts_insertion_dinuc_ntbias_counter& _counter):
  proba(),//zeros(model.max_assignments, 1);
  event_probability(),//zeros(model.max_assignments,1);
  n_assignments(0),//0, number of valid assignments
  skips(0),//0, skipped sequences
  likelihood(0),//0, total probability of generating the sequence (sum over ssignments
  generation_probability(0),//0, total probability of generating the sequence (assuming zero error rate);
  max_proba_index(0),//0, index of the best one
  
  //-------nP****
  //% If counter variable begins with 'nP' i.e. the model variable is a probability distribution P(x),
  //% then the assigns variable should store the index of the value of x, for each assignment for a sequence.
  //% eg. counter.nPpVdelV -> assigns.pVdelV which stores the indices of palindrome half-length and deletions for V as [1 + npV, 1 + ndV]

  //-------nM****
  //% If counter variable begins with 'nM' i.e. the model variable is a mean value of x,
  //% then the assigns variable should store the value of x, for each assignment for a sequence.
  //% eg. counter.nMerror_vs_position -> assigns.error_vs_position which stores number of errors at each position for each assignment for a sequence.
  
  insVD(),
  insDJ(),

  cutV_given_V(),
  cutJ_given_J(),
  cutDlcutDr_given_D(),

  V(),
  DJ(),//Joint P(V, D, J gene choices)
  Vallele_given_gene(),//Probabilities of alleles given gene for each gene
  Dallele_given_gene(),
  Jallele_given_gene(),

  VD_left_edge_dinucleotide(),// = zeros(4,4),
  VD_right_edge_dinucleotide(),// = zeros(4,4);
  DJ_left_edge_dinucleotide(),// = zeros(4,4);
  DJ_right_edge_dinucleotide(),// = zeros(4,4);

  pVmax_delV_V(),// =  zeros(model.max_palindrome + 1, model.max_V_deletions + 1, size(model.PV,1));
  pJmax_delJ_J(),// = zeros(model.max_palindrome + 1, model.max_J_deletions + 1, size(model.PDJ,2));
  
  pDlmax_delDl_D(),// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));
  pDrmax_delDr_D(),// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));

  VDJ(),// = zeros(size(model.PV,1), size(model.PDJ,1), size(model.PDJ,2));
  pVdelV(),// = zeros( model.max_palindrome + 1, model.max_V_deletions + 1);
  pDldelDl(),// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
  pDrdelDr(),// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
  pJdelJ(),// = zeros( model.max_palindrome + 1, model.max_J_deletions + 1);


  VcutV(),// = zeros(size(model.PV,1),model.max_palindrome + model.max_V_deletions + 1);
  DcutDl(),// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
  DcutDr(),// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
  JcutJ(),// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_J_deletions + 1);

  DcutV(),// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_V_deletions + 1);
  DcutJ(),// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_J_deletions + 1);

  VcutJ(),// = zeros(size(model.PV,1),model.max_palindrome + model.max_J_deletions + 1);
  VcutDl(),// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);
  VcutDr(),// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);

  JcutDl(),// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
  JcutDr(),// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
  JcutV(),// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_V_deletions + 1);

  insVDcutV(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);
  insDJcutV(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);

  insVDcutDl(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  insDJcutDl(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

  insVDcutDr(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  insDJcutDr(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  
  insVDcutJ(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);
  insDJcutJ(),// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);

  /*
  unsigned max_V_length;// = model.read_length;
  unsigned max_D_length;// = 16;
  unsigned max_J_length;// = 18;
  */
  V_align_length(),// = zeros(max_V_length + 1,1);
  D_align_length(),// = zeros(max_D_length + 1,1);
  J_align_length(),// = zeros(max_J_length + 1,1);
  
  JJ_align_length(),// = zeros(size(model.PDJ,2), 1 + max_J_length);
  VV_align_length(),// = zeros(size(model.PV,1), 1 + max_V_length);
  

  insDJ_D_align_length(),// = zeros(model.max_insertions +1, max_D_length + 1);
  insVD_D_align_length(),// = zeros(model.max_insertions +1, max_D_length + 1);
  
  insDJ_J_align_length(),// = zeros(model.max_insertions +1, max_J_length + 1);
  insVD_J_align_length(),// = zeros(model.max_insertions +1, max_J_length + 1);

  insDJ_V_align_length(),// = zeros(model.max_insertions +1, max_V_length + 1);
  insVD_V_align_length(),// = zeros(model.max_insertions +1, max_V_length + 1);

  Dallele_D_align_length(),// = zeros(3, max_D_length + 1);


  delVinsVD(),// = zeros(model.max_V_deletions +1, model.max_insertions +1);
  delVinsDJ(),// = zeros(model.max_V_deletions +1, model.max_insertions +1);
  delVdelDl(),// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
  delVdelDr(),// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
  delVdelJ(),// = zeros(model.max_V_deletions +1, model.max_J_deletions +1);

   delJinsVD(),// = zeros(model.max_J_deletions +1, model.max_insertions +1);
   delJinsDJ(),//zeros(model.max_J_deletions +1, model.max_insertions +1);
   delJdelDl(),//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
   delJdelDr(),//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
  
   delDlinsVD(),//zeros(model.max_D_deletions +1, model.max_insertions +1);
   delDlinsDJ(),//zeros(model.max_D_deletions +1, model.max_insertions +1);
   delDldelDr(),//zeros(model.max_D_deletions +1, model.max_D_deletions +1);
  
   delDrinsVD(),//zeros(model.max_D_deletions +1, model.max_insertions +1);
   delDrinsDJ(),//zeros(model.max_D_deletions +1, model.max_insertions +1);
  
   insVDinsDJ(),//zeros(model.max_insertions +1, model.max_insertions +1);
  
   VdelV(),//zeros(size(model.PV,1), model.max_V_deletions+1 );
   DdelDl(),//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
   DdelDr(),//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
   JdelJ(),//zeros(size(model.PDJ,2), model.max_J_deletions+1 );
  
   VinsVD(),//zeros(size(model.PV,1), model.max_insertions+1 );
   DinsVD(),//zeros(size(model.PDJ,1), model.max_insertions+1 );
   DinsDJ(),//zeros(size(model.PDJ,1), model.max_insertions+1 );
   JinsDJ(),//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
  VdelDl(),//zeros(size(model.PV,1), model.max_D_deletions+1 );
   VdelDr(),//zeros(size(model.PV,1), model.max_D_deletions+1 );
   VdelJ(),//zeros(size(model.PV,1), model.max_J_deletions+1 );
  
  JdelV(),//zeros(size(model.PDJ,2), model.max_V_deletions+1 );
   JdelDl(),//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
   JdelDr(),//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
  
   DdelV(),//zeros(size(model.PDJ,1), model.max_V_deletions+1 );
   DdelJ(),//zeros(size(model.PDJ,1), model.max_J_deletions+1 );
  
   VinsDJ(),//zeros(size(model.PV,1), model.max_insertions+1 );
   JinsVD(),//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
   pVinsVD(),//zeros(model.max_palindrome +1, model.max_insertions +1);
   pVinsDJ(),//zeros(model.max_palindrome +1, model.max_insertions +1);
   pVdelDl(),//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   pVdelDr(),//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   pVdelJ(),//zeros(model.max_palindrome +1, model.max_J_deletions +1);
   VpV(),//zeros(size(model.PV,1), model.max_palindrome+1 );
   JpV(),//zeros(size(model.PDJ,2), model.max_palindrome+1 );
   DpV(),//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
   pJinsVD(),//zeros(model.max_palindrome +1, model.max_insertions +1);
   pJinsDJ(),//zeros(model.max_palindrome +1, model.max_insertions +1);
   pJdelV(),//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  pJdelDl(),//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   pJdelDr(),//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   VpJ(),//zeros(size(model.PV,1), model.max_palindrome+1 );
   JpJ(),//zeros(size(model.PDJ,2), model.max_palindrome+1 );
   DpJ(),//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
   pDlinsVD(),//zeros(model.max_palindrome +1, model.max_insertions +1);
   pDlinsDJ(),//zeros(model.max_palindrome +1, model.max_insertions +1);
   pDldelV(),//zeros(model.max_palindrome +1, model.max_V_deletions +1);
   pDldelJ(),//zeros(model.max_palindrome +1, model.max_J_deletions +1);
   pDldelDr(),//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   VpDl(),//zeros(size(model.PV,1), model.max_palindrome+1 );
   JpDl(),//zeros(size(model.PDJ,2), model.max_palindrome+1 );
   DpDl(),//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
   pDrinsVD(),//zeros(model.max_palindrome +1, model.max_insertions +1);
   pDrinsDJ(),//zeros(model.max_palindrome +1, model.max_insertions +1);
   pDrdelV(),//zeros(model.max_palindrome +1, model.max_V_deletions +1);
   pDrdelJ(),//zeros(model.max_palindrome +1, model.max_J_deletions +1);
   pDrdelDl(),//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   VpDr(),//zeros(size(model.PV,1), model.max_palindrome+1 );
   JpDr(),//zeros(size(model.PDJ,2), model.max_palindrome+1 );
   DpDr(),//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  pVpDl(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
   pVpDr(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
   pVpJ(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
   pDlpDr(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
   pDlpJ(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
   pDrpJ(),//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  
  //%------nM*** fields
      
  mononucleotideVD(),// zeros(4,1);
  mononucleotideDJ(),// = zeros(4,1);
  insertionVD(),// = zeros(4,1);
  insertionDJ(),// = zeros(4,1);

  
  trinucleotideVD(),// = zeros(4,4,4);
  trinucleotideDJ(),// = zeros(4,4,4);

  VV_err_pos(),// = zeros(size(model.PV,1), max_V_length + model.max_V_deletions);
  JJ_err_pos(),// = zeros(size(model.PDJ,2), max_J_length);

  zeroD(), 
  nucleotideVD(),  nucleotideVD_5prime(),//nucleotide distr's
  nucleotideDJ(), nucleotideDJ_3prime(),

  error(), sequenced_nucleotide (),//error rate, will be 

  error_vs_position(),//!!this is nM zeros(model.read_length,1);
  coverage()//zeros(model.read_length,1);
{
  //do initialization
  //==========================>member declaration from counter
  unsigned dim_size[]={_model.max_assignments};
  proba.initialize(1, dim_size, 0.0); //zeros(model.max_assignments, 1);
  event_probability.initialize(1, dim_size, 0.0); //zeros(model.max_assignments,1);
    
  //-------nP****
  //% If counter variable begins with 'nP' i.e. the model variable is a probability distribution P(x),
  //% then the assigns variable should store the index of the value of x, for each assignment for a sequence.
  //% eg. counter.nPpVdelV -> assigns.pVdelV which stores the indices of palindrome half-length and deletions for V as [1 + npV, 1 + ndV]

  //-------nM****
  //% If counter variable begins with 'nM' i.e. the model variable is a mean value of x,
  //% then the assigns variable should store the value of x, for each assignment for a sequence.
  //% eg. counter.nMerror_vs_position -> assigns.error_vs_position which stores number of errors at each position for each assignment for a sequence.
  
   insVD.initialize(1, dim_size, -1);
   insDJ.initialize(1, dim_size, -1);
   unsigned dim_size2[2]={_model.max_assignments, 2};
   cutV_given_V.initialize(2, dim_size2, -1);
   cutJ_given_J.initialize(2, dim_size2, -1);
   dim_size2[1]=3;
   cutDlcutDr_given_D.initialize(2, dim_size2, -1);;

   V.initialize(1, dim_size, -1);
   dim_size2[1]=2;
   DJ.initialize(2, dim_size2, -1); //Joint P(V, D, J gene choices)
   Vallele_given_gene.initialize(2, dim_size2, -1); //Probabilities of alleles given gene for each gene
   Dallele_given_gene.initialize(2, dim_size2, -1);
   Jallele_given_gene.initialize(2, dim_size2, -1);

   VD_left_edge_dinucleotide.initialize(2, dim_size2, 0.0);// = zeros(4,4);
   VD_right_edge_dinucleotide.initialize(2, dim_size2, 0.0);// = zeros(4,4);
   DJ_left_edge_dinucleotide.initialize(2, dim_size2, 0.0);// = zeros(4,4);
   DJ_right_edge_dinucleotide.initialize(2, dim_size2, 0.0);// = zeros(4,4);
  

   dim_size2[1]=3;
   pVmax_delV_V.initialize(2, dim_size2, -1);// =  zeros(model.max_palindrome + 1, model.max_V_deletions + 1, size(model.PV,1));
   pJmax_delJ_J.initialize(2, dim_size2, -1);// = zeros(model.max_palindrome + 1, model.max_J_deletions + 1, size(model.PDJ,2));

   pDlmax_delDl_D.initialize(2, dim_size2, -1);// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));
   pDrmax_delDr_D.initialize(2, dim_size2, -1);// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));

   VDJ.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1), size(model.PDJ,1), size(model.PDJ,2));
   dim_size2[1]=2;
   pVdelV.initialize(2, dim_size2, -1);// = zeros( model.max_palindrome + 1, model.max_V_deletions + 1);
   pDldelDl.initialize(2, dim_size2, -1);// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
   pDrdelDr.initialize(2, dim_size2, -1);// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
   pJdelJ.initialize(2, dim_size2, -1);// = zeros( model.max_palindrome + 1, model.max_J_deletions + 1);


   VcutV.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1),model.max_palindrome + model.max_V_deletions + 1);
   DcutDl.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
   DcutDr.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
   JcutJ.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_J_deletions + 1);

   DcutV.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_V_deletions + 1);
   DcutJ.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_J_deletions + 1);

   VcutJ.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1),model.max_palindrome + model.max_J_deletions + 1);
   VcutDl.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);
   VcutDr.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);

   JcutDl.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
   JcutDr.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
   JcutV.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_V_deletions + 1);

   insVDcutV.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);
   insDJcutV.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);

   insVDcutDl.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
   insDJcutDl.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

   insVDcutDr.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
   insDJcutDr.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

   insVDcutJ.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);
   insDJcutJ.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);

  /*
  unsigned max_V_length;// = model.read_length;
  unsigned max_D_length;// = 16;
  unsigned max_J_length;// = 18;
  */
   V_align_length.initialize(1, dim_size, -1);// = zeros(max_V_length + 1,1);
   D_align_length.initialize(1, dim_size, -1);// = zeros(max_D_length + 1,1);
   J_align_length.initialize(1, dim_size, -1);// = zeros(max_J_length + 1,1);

   JJ_align_length.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,2), 1 + max_J_length);
   VV_align_length.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1), 1 + max_V_length);


   insDJ_D_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_D_length + 1);
   insVD_D_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_D_length + 1);

   insDJ_J_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_J_length + 1);
   insVD_J_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_J_length + 1);

   insDJ_V_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_V_length + 1);
   insVD_V_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_V_length + 1);

   Dallele_D_align_length.initialize(2, dim_size2, -1);// = zeros(3, max_D_length + 1);


   delVinsVD.initialize(2, dim_size2, -1);// = zeros(model.max_V_deletions +1, model.max_insertions +1);
   delVinsDJ.initialize(2, dim_size2, -1);// = zeros(model.max_V_deletions +1, model.max_insertions +1);
   delVdelDl.initialize(2, dim_size2, -1);// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
   delVdelDr.initialize(2, dim_size2, -1);// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
   delVdelJ.initialize(2, dim_size2, -1);// = zeros(model.max_V_deletions +1, model.max_J_deletions +1);

   delJinsVD.initialize(2, dim_size2, -1);// = zeros(model.max_J_deletions +1, model.max_insertions +1);
   delJinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_J_deletions +1, model.max_insertions +1);
   delJdelDl.initialize(2, dim_size2, -1);//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
   delJdelDr.initialize(2, dim_size2, -1);//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
  
   delDlinsVD.initialize(2, dim_size2, -1);//zeros(model.max_D_deletions +1, model.max_insertions +1);
   delDlinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_D_deletions +1, model.max_insertions +1);
   delDldelDr.initialize(2, dim_size2, -1);//zeros(model.max_D_deletions +1, model.max_D_deletions +1);
  
   delDrinsVD.initialize(2, dim_size2, -1);//zeros(model.max_D_deletions +1, model.max_insertions +1);
   delDrinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_D_deletions +1, model.max_insertions +1);
  
   insVDinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_insertions +1, model.max_insertions +1);
  
   VdelV.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_V_deletions+1 );
   DdelDl.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
   DdelDr.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
   JdelJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_J_deletions+1 );
  
   VinsVD.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_insertions+1 );
   DinsVD.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_insertions+1 );
   DinsDJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_insertions+1 );
   JinsDJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
   VdelDl.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_D_deletions+1 );
   VdelDr.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_D_deletions+1 );
   VdelJ.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_J_deletions+1 );
  
   JdelV.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_V_deletions+1 );
   JdelDl.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
   JdelDr.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
  
   DdelV.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_V_deletions+1 );
   DdelJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_J_deletions+1 );
  
   VinsDJ.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_insertions+1 );
   JinsVD.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
   pVinsVD.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
   pVinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
   pVdelDl.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   pVdelDr.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   pVdelJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_J_deletions +1);
   VpV.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_palindrome+1 );
   JpV.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_palindrome+1 );
   DpV.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
   pJinsVD.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
   pJinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
   pJdelV.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_V_deletions +1);
   pJdelDl.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   pJdelDr.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   VpJ.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_palindrome+1 );
   JpJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_palindrome+1 );
   DpJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
   pDlinsVD.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
   pDlinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
   pDldelV.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_V_deletions +1);
   pDldelJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_J_deletions +1);
   pDldelDr.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   VpDl.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_palindrome+1 );
   JpDl.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_palindrome+1 );
   DpDl.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
   pDrinsVD.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
   pDrinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
   pDrdelV.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_V_deletions +1);
   pDrdelJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_J_deletions +1);
   pDrdelDl.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
   VpDr.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_palindrome+1 );
   JpDr.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_palindrome+1 );
   DpDr.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
   pVpDl.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
   pVpDr.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
   pVpJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
   pDlpDr.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
   pDlpJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
   pDrpJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  
  //%------nM*** fields
   dim_size2[1]=4;
   mononucleotideVD.initialize(2, dim_size2, 0.0);// zeros(4,1);
   mononucleotideDJ.initialize(2, dim_size2, 0.0);// = zeros(4,1);
   insertionVD.initialize(2, dim_size2, 0.0);// = zeros(4,1);what this is?? is this insertion or nucleotide dist'n among insertion??
  insertionDJ.initialize(2, dim_size2, 0.0);// = zeros(4,1);


  unsigned dim_size4[4]={_model.max_assignments, 4, 4, 4};

  trinucleotideVD.initialize(4, dim_size4, 0.0);// = zeros(4,4,4);
  trinucleotideDJ.initialize(4, dim_size4, 0.0);// = zeros(4,4,4);
  unsigned dim_size3[3]={_model.max_assignments, _model.number_V_genes, _model.model_params.max_V_depth};
  VV_err_pos.initialize(3, dim_size3, 0.0);// = zeros(size(model.PV,1), max_V_length + model.max_V_deletions);
  dim_size3[2]=_model.model_params.max_J_depth;
  dim_size3[1]=_model.number_J_genes;
  JJ_err_pos.initialize(3, dim_size3, 0.0);// = zeros(size(model.PDJ,2), max_J_length);

  
  zeroD.initialize(1, dim_size, 0.0) ;
  dim_size3[1]=4;
  dim_size3[2]=4;
  nucleotideVD.initialize(3, dim_size3, 0.0); nucleotideVD_5prime.initialize(3, dim_size3, 0.0);//nucleotide distr's
  nucleotideDJ.initialize(3, dim_size3, 0.0); nucleotideDJ_3prime.initialize(3, dim_size3, 0.0);

  error.initialize(1, dim_size,0.0); sequenced_nucleotide.initialize(1,dim_size, 0.0);//error rate, will be 
  
  dim_size2[1]=_model.model_params.maximum_read_length;
  error_vs_position.initialize(2,dim_size2,0.0);//!!this is nM zeros(model.read_length,1);
  coverage.initialize(2,dim_size2,0.0);//zeros(model.read_length,1);
}


VDJ_cuts_insertion_dinuc_ntbias_assigns::~VDJ_cuts_insertion_dinuc_ntbias_assigns()
{
  //empty
}


