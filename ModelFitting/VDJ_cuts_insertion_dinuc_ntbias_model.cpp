#include <cmath>

#include "VDJ_cuts_insertion_dinuc_ntbias_model.hpp"

#include "VDJ_cuts_insertion_dinuc_ntbias_counter.hpp"
#include "VDJ_cuts_insertion_dinuc_ntbias_assigns.hpp"

//NOTE: to do: do we need to add the parameter such max_assignments, max_insertions, etc
//in the alignmentSettings (?) or some other way to specify these values!!!

VDJ_cuts_insertion_dinuc_ntbias_model::VDJ_cuts_insertion_dinuc_ntbias_model
(const GenomicV* _genV, const unsigned& _numV, 
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ)
  :
  /*initilization list*/
  max_assignments(6000), max_insertions(30),
  max_V_deletions(16), max_D_deletions(16), max_J_deletions(18),
  number_V_genes(0), number_D_genes(0), number_J_genes(0),
  max_V_n_alleles(0), max_D_n_alleles(0), max_J_n_alleles(0),
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
  
  PcutV_given_V(), PcutJ_given_J(), PcutDlcutDr_given_D(),
  PV(), PDJ(), PVallele_given_gene(), PDallele_given_gene(),
  
  Rerror_per_sequenced_nucleotide(0)
{

  //need to do things to initialize the parameters and arrays
  min_V_cut=-1*max_palindrome;
  min_D_cut=-1*max_palindrome;
  min_J_cut=-1*max_palindrome;
  
  max_V_cut=max_V_deletions;
  max_D_cut=max_D_deletions;
  max_J_cut=max_J_deletions;

  //start intialize the various matrix
  unsigned maxGeneIndex_V=max_gene_index(_genV,_numV);
  unsigned maxGeneIndex_D=max_gene_index(_genD,_numD);
  unsigned maxGeneIndex_J=max_gene_index(_genJ,_numJ);
  number_V_genes=maxGeneIndex_V;
  number_D_genes=maxGeneIndex_D;
  number_J_genes=maxGeneIndex_J;

  unsigned dim_size[1]={max_insertions+1};
  PinsVD.initialize(1, dim_size, 1);
  PinsDJ.initialize(1, dim_size, 1);
  dim_size[0]=maxGeneIndex_D;
  PV.initialize(1, dim_size, 1);
  unsigned dim_size2[2]={maxGeneIndex_D, maxGeneIndex_J};
  PDJ.initialize(2, dim_size2, 1);

  dim_size2[1]=(unsigned)(max_V_cut-min_V_cut+1);
  dim_size2[0]=maxGeneIndex_V;
  PcutV_given_V.initialize(2, dim_size2,1);

  dim_size2[1]=(unsigned)(max_J_cut-min_J_cut+1);
  dim_size2[0]=maxGeneIndex_J;
  PcutJ_given_J.initialize(2,dim_size2, 1);
  
  unsigned dim_size3[3]={maxGeneIndex_D, (unsigned)(max_D_cut-min_D_cut+1), (unsigned)(max_D_cut-min_D_cut+1) };
  PcutDlcutDr_given_D.initialize(3, dim_size3, 1);

  //also go through to set some impossible cases impossible
  //case where cut on both side together is longer than sequence
unsigned lseq;
  for(unsigned d=0;d<maxGeneIndex_D;d++)
    {
      lseq = _genD[d].Get_Seq().GetLength();
      for(int  cutDl=0;cutDl<=max_D_cut;cutDl++)
	{
	  for(int  cutDr=0;cutDr<=max_D_cut;cutDr++)
	    {
	      if (((unsigned)(cutDl + cutDr)) > lseq)
		PcutDlcutDr_given_D(d, (unsigned)(cutDl-  min_D_cut) ,(unsigned)( cutDr - min_D_cut)) = 0;
	      //end
	    }// end
	}//       end
    }//       end

  //set up rate for what? still don't know at this point
  dim_size2[0]=4;
  dim_size2[1]=4;
  RnucleotideVD_per_nucleotideVD_5prime.initialize(2, dim_size2,1.0/4.0);
  //RnucleotideVD_per_nucleotideVD_5prime=RnucleotideVD_per_nucleotideVD_5prime/4;
  RnucleotideDJ_per_nucleotideDJ_3prime.initialize(2, dim_size2,1.0/4.0);

  //% Initialize probabilities of alleles given genes, to 1/n_alleles for each gene.
  //V_genes = max([genV.gene_index]);
  max_V_n_alleles = max_n_alleles(_genV, _numV);
  dim_size2[0]=maxGeneIndex_V;
  dim_size2[1]=max_V_n_alleles;
PVallele_given_gene.initialize(2, dim_size2, 0.0);
  for(unsigned a=0;a<_numV;a++)
    {
PVallele_given_gene(_genV[a].Get_GeneIndex(), _genV[a].Get_Allele()) = 
	1.0/_genV[a].Get_n_alleles();
    }       //end
  
  //D_genes = max([genD.gene_index]);
  max_D_n_alleles = max_n_alleles(_genD, _numD);//...n_alleles]);
  max_J_n_alleles =max_n_alleles(_genJ, _numJ);
  dim_size2[0]=maxGeneIndex_D;
  dim_size2[1]=max_D_n_alleles;
  PDallele_given_gene.initialize(2, dim_size2,0.0);// zeros(max_D_alleles , D_genes);
  for(unsigned a=0;a<_numD;a++)
    {
      PDallele_given_gene(_genD[a].Get_GeneIndex(), _genD[a].Get_Allele()) = 1.0/_genD[a].Get_n_alleles();
    }
  
Rerror_per_sequenced_nucleotide = 1.0E-7; //% error rate
//% Normalizes all distributions
//doing first round of functions
  Normalize();
  CalculateAssignmentEntropies();
}

VDJ_cuts_insertion_dinuc_ntbias_model::~VDJ_cuts_insertion_dinuc_ntbias_model
()
{
//empty so far
}

bool VDJ_cuts_insertion_dinuc_ntbias_model::ValidateModel()
{
//here inside this one, we simple need to make sure all the
  //fields has been initialized and not empty matrix and 
  //no unwanted default values
  
  //at this point, we don't do much, since we rely on the constructor to 
  //initialize/populate everything

return true;
}

bool VDJ_cuts_insertion_dinuc_ntbias_model::Normalize()//normalize the model
{
  
//Normalizes all probability distributions in model.
//NOT, the R** fields (rate fields)
  //norm_model = model;
  //          ps = fieldnames(norm_model);
  //          for p=1:length(ps)
  //              if ps{p}(1)=='P'
  //                  if isempty(strfind(ps{p},'_given_'))
  
  //=======================
  //for all the P*** field (not _given_ ones) in the model, we simply add all items
  //and then normlize it
  PinsVD=PinsVD/sum_all(PinsVD);
  PinsDJ=PinsDJ/sum_all(PinsDJ);
  PV=PV/sum_all(PV);
  PDJ=PDJ/sum_all(PDJ);

  //                      norm_model.(ps{p}) = norm_model.(ps{p})/sum( reshape(norm_model.(ps{p}),[],1) );

  //================
  // If it's a conditional probability, then normalize appropriately
  //nd = ndims(norm_model.(ps{p}));
  //                      norm_model.(ps{p}) = conditional_distribution(norm_model.(ps{p}), nd);

  //
  Matrix<double> GeneTotal=sum(PcutV_given_V, 1);//sum the dim 1 (second dimension, since it starts at 0) that pointing to cutV, so to get geneTotal
  //the return matrix is of dimension _numV x 1 (vector).
  PcutV_given_V.divide_by_dimension(GeneTotal, 0);//divide along dimension 0 (first) to normalize
  GeneTotal=sum(PcutJ_given_J, 1);
  PcutJ_given_J.divide_by_dimension(GeneTotal, 0);

  Matrix<double> GeneTotal_temp=sum(PcutDlcutDr_given_D,2);
  GeneTotal=sum(GeneTotal_temp,1);
  PcutDlcutDr_given_D.divide_by_dimension(GeneTotal, 0);
  
  GeneTotal=sum(PVallele_given_gene,1);
  PVallele_given_gene.divide_by_dimension(GeneTotal, 0);

  GeneTotal=sum(PDallele_given_gene,1);
  PDallele_given_gene.divide_by_dimension(GeneTotal, 0);
  return true;
}

void  VDJ_cuts_insertion_dinuc_ntbias_model::InitializeCounter(Counter& _c) const
{
//we don't do anything
  //the actual job was done by constructor of the Counter;
  //return ;
}

void VDJ_cuts_insertion_dinuc_ntbias_model::GetModelFromCounter(const Counter& _c) 
{
}

void VDJ_cuts_insertion_dinuc_ntbias_model::InitializeAssign(Assigns& _a) const
{
  //we don't do anything
  //the actual job was done by constructor of the assign.
  
}

void VDJ_cuts_insertion_dinuc_ntbias_model::UpdateCounter(Assigns& _a, Counter& _c) const
{

  VDJ_cuts_insertion_dinuc_ntbias_assigns& vdj_assigns=
  dynamic_cast<VDJ_cuts_insertion_dinuc_ntbias_assigns&>(_a);
  unsigned n_as=vdj_assigns.n_assignments;

  if(n_as>0 && vdj_assigns.likelihood>0)
    {
      //first we simply create a temporary counter
      VDJ_cuts_insertion_dinuc_ntbias_counter delta_counter(*this);
      
      delta_counter.logLikelihood=log(vdj_assigns.likelihood);

      delta_counter.N_processed=1;

      //normalize the probilitties of assignments
      vdj_assigns.proba=vdj_assigns.proba/vdj_assigns.likelihood;

      //now we want to update every relevant fields
      //first, nP**** fields.
      //--nPinsVD
      double tempContribution=0;
 
      for(unsigned i=0;i<n_as;i++)//for each valid assignment
	{
	  if(((signed)insVD(i))!=-1)
	    delta.nPinsVD(insVD(i))+=vdj_assigns.prob((insVD(i));
	}
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


   VD_left_edge_dinucleotide.initialize(2, dim_size2, -1);// = zeros(4,4);
   VD_right_edge_dinucleotide.initialize(2, dim_size2, -1);// = zeros(4,4);
   DJ_left_edge_dinucleotide.initialize(2, dim_size2, -1);// = zeros(4,4);
   DJ_right_edge_dinucleotide.initialize(2, dim_size2, -1);// = zeros(4,4);

   dim_size2[1]=3;
   pVmax_delV_V.initialize(2, dim_size2, -1);// =  zeros(model.max_palindrome + 1, model.max_V_deletions + 1, size(model.PV,1));
   pJmax_delJ_J.initialize(2, dim_size2, -1);// = zeros(model.max_palindrome + 1, model.max_J_deletions + 1, size(model.PDJ,2));

   pDlmax_delDl_D.initialize(2, dim_size2, -1);// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));
   pDrmax_delDr_D.initialize(2, dim_size2, -1);// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));

   VDJ.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1), size(model.PDJ,1), size(model.PDJ,2));
   dim_size[1]=2;
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
  
    }
}

  //!!!here this one needs testing............===???????????????
  //the function used to update the nP* field in the counter
  //index_fields, fields in the assigns, contains the information about
  //    the assigns that will be write to counter _fields_c
  //
  //      index_field and prob field are from assign, and _fields_c is from counter
  //return false if the input matrix is not what they should , such as the dimensions
  //    or the size are not appropriate
  bool VDJ_cuts_insertion_dinuc_ntbias_model::Update_nP_field
    (const Matrix<unsigned>& _index_field_a, const Matrix<double>& _prob_field_a, 
     Matrix<double>& _fields_c, const unsigned& _num_valid_assignments) const
  {
      bool flag=true;
      //first check for the correct input
      //_index_field_a is max_assignments x index1 x index2...x index_n
      //_prob_field_a is max_assignments x 1
      //_fields_c =index1 x index2 ... x index_n
      
      if(_index_field.size(0)<_num_valid_assignments)
	{
	  cout<<"Error: pass in not number of valid assignments or incorrect sized index fields, please check"<<endl;
	  return false;
	}
      
      if(_prob_fields_a.dim()!=1)
	{
	  cout<<"Error: pass in not correct prob_field, please check"<<endl;
	  return false;
	}
      if(_index_field_a.dim()<1||_index_field_a.dim()>2) //so far only support 1 or 2
	{
	  cout<<"Error: pass in not correct index_field (dimension not 1 or 2), please check"<<endl;
	  return false;
	}

      if(_prob_field_a.size(0)!=_index_field_a.size(0))
	{
	  cout<<"Error: the size of two fields from assing don't match, please check"<<endl;
	  return false;
	}
      if(_index_field_a.dim()==2 &&
	 (_index_field_a.size(1)!=_fields_c.dim())
	)
	{
	  cout<<"Error: the dimension of field from counter is not right, please check"<<endl;
	  return false;
	}
      
      /*
      for(unsigned i=1;i<_index_field_a.dim();i++)
	{
	  //check for the size
	  if(_index_field_a.size(i)!=_field_c.size(i-1))
	    {
	      cout<<"Error: the sizes for index and field don't match, please check"<<endl;
	      return false;
	    }
	    }*/
      
      //phew, finally, we can do the job
      //first, let's figure out num_all_good_assigns in this index of
      //valid assignment. A good assigns is the valid assign plus all the
      //indices have been assigned. no not assigned index.
      //so we can also do the correct normalization of proba
      unsigned blockSize=_field_c.nTotal();//total number of elements in the counter field
      unsigned* k_all_good_assigns=new unsigned[_num_valid_assignments];
      unsigned num_all_good_assigns=0;
      double totalProb=0;
      
      for(unsigned i=0;i<_num_valid_assignments;i++)
	{
	  
	  if(_index_field_a.dim()==1)
	    {
	      //check for good assignment
	      if(((signed)_index_field_a(i))==-1)
		    {
		      
		      break;
		    }
	      else //not a bad one, so
		{
		  //set the values, if we are here, this one is good
		      k_all_good_assigns[num_all_good_assigns]=i;
		      num_all_good_assigns++;
		      totalProb+=_prob_field_a(i);
		}
	      
	    }
	  else //2, can only be two in this case
	    {
	  
	      for(unsigned j=0;j<_index_field_a.size(1);j++)
		{
		  if(((signed)_index_field_a(i,j))==-1)
		    {
		      //badAssign=true;
		      break;
		    }
		  if(j==_index_field_a.size(1)-1)
		    {
		      //set the values, if we are here, this one is good
		      k_all_good_assigns[num_all_good_assigns]=i;
		      num_all_good_assigns++;
		      totalProb+=_prob_field_a(i);
		    }
		}
	    }
	}
      //now we do the updating
      unsigned* position_index_field_a=new unsigned[_index_field_a.dim()];
      //_index_field_a dimension has to be larger than 0. 1 or 2 only, can not be bigger than 2
      //
      unsigned* position_field_c=new unsigned[_field_c.dim()];
      for(unsigned i=0;i<num_all_good_assigns;i++)
	{
	  if(_index_field_a.dim()==1)
	    {
	      _field_c(_index_field_a(k_all_good_assigns[i]))+=_prob_field_a(k_all_good_assigns[i]);
	    }
	  else //case of dimension =2, but at this dimension 2 it could 2 or 3 or bigger, which is the dimension of the field.dim()
	    {
	      position_index_field_a[0]=k_all_good_assigns[i];
	      position_field_c[0]=0;//this is in case of scalar matrix, it is possible
	      //for(unsigned j=1;j<_index_field_c.dim();j++)
	      //	{
	      for(unsigned k=0;k<_index_field_c.size(1);k++)
		{
		  position_index_field_a[1]=k;
		  //figure out the indices of the specific entry
		  position_field_c[k]=_index_field_a(position_index_field_a, _index_field_a.dim()); 
		  
		}
	      //now we got position, let's up 
	      _field_c(position)+=_prob_field_a(k_all_good_assigns[i])/totalProb;
	      //}
	    }//end of else
	}//end of for updating

      //clear mem
      delete [] k_all_good_assigns;
      delete [] position_field_c;
      delete [] position_index_field_a;

      
      return flag;

  }
  //the function used to update the nM* field in the counter
  //count_fields, fields in the assigns, contains the information about
  //    the assigns that will be write to counter _fields_c
  //
  //      count_field and prob field are from assign, and _fields_c is from counter
  //
  //the calculation is that we add expectation value from this set
  //of assignments for this sequence to counter variable
  //
  //return false if the input matrix is not what they should , such as the dimensions
  //    or the size are not appropriate
  bool VDJ_cuts_insertion_dinuc_ntbias_model::Update_nM_field
    (Matrix<unsigned> _count_field_a, Matrix<double> _prob_field_a, 
     Matrix<double> _fields_c, const unsigned& _num_valid_assignments)
  {
    bool flag=true;
    //first check for the validity of the input
    //_count_field_a is max_assignments x index1 x index2...x index_n
      //_prob_field_a is max_assignments x 1
      //_fields_c =index1 x index2 ... x index_n
      
      if(_count_field.size(0)<_num_valid_assignments)
	{
	  cout<<"Error: pass in not number of valid assignments or incorrect sized index fields, please check"<<endl;
	  return false;
	}
      
      if(_prob_fields_a.dim()!=1)
	{
	  cout<<"Error: pass in not correct prob_field, please check"<<endl;
	  return false;
	}
      if(_count_field_a.dim()!=_field_c.dim()+1) //_count_field is one more dimension than _field_a
	{
	  cout<<"Error: pass in not correct index_field (dimension not 1 or 2), please check"<<endl;
	  return false;
	}

      if(_prob_field_a.size(0)!=_count_field_a.size(0))
	{
	  cout<<"Error: the size of two fields from assing don't match, please check"<<endl;
	  return false;
	}
      /*if(_index_field_a.dim()==2 &&
	 (_index_field_a.size(1)!=_fields_c.dim())
	)
	{
	  cout<<"Error: the dimension of field from counter is not right, please check"<<endl;
	  return false;
	  }*/
      
      
      for(unsigned i=1;i<_count_field_a.dim();i++)
	{
	  //check for the size
	  if(_count_field_a.size(i)!=_field_c.size(i-1))
	    {
	      cout<<"Error: the sizes for index and field don't match, please check"<<endl;
	      return false;
	    }
	}
      //now we are good, start doing the updating
      //now go through each element of field of each assignment/here it is a "block"
      unsigned blockSize=_field_c.nTotal();//total number of elements in the counter field
      for(unsigned i=0;i<_num_valid_assignments;i++)
	{
	  for(unsigned j=0;j<_block_size;j++)
	    {
	      double temp=_field_c.Set1DArrayElement(j)+_count_field_a.Get1DArrayElement(i*blockSize+j)*_prob_field_a(i);
	      _field_c.Set1DArrayElement(j, temp);
	    }
	}

      //done
      return flag;
  }//end of function


void VDJ_cuts_insertion_dinuc_ntbias_model::SumCounters(const Counter& _c1, const Counter& _c2, Counter& _retC) const
{
  //need to dynamic casting everything
  
}

void VDJ_cuts_insertion_dinuc_ntbias_model::CalculateAssignmentEntropies()
{
  //=========>to be implemented, for now leave it blank;
}
