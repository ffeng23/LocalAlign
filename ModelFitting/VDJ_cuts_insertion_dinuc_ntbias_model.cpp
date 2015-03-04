

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
  unsigned max_V_alleles = max_n_alleles(_genV, _numV);
  dim_size2[0]=maxGeneIndex_V;
  dim_size2[1]=max_V_alleles;
  PVallele_given_gene.initialize(2, dim_size2, 0.0);
  for(unsigned a=0;a<_numV;a++)
    {
      PVallele_given_gene(_genV[a].Get_GeneIndex(), _genV[a].Get_Allele()) = 
	1.0/_genV[a].Get_n_alleles();
    }       //end
  
  //D_genes = max([genD.gene_index]);
  unsigned max_D_alleles = max_n_alleles(_genD, _numD);//...n_alleles]);
  dim_size2[0]=maxGeneIndex_D;
  dim_size2[1]=max_D_alleles;
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

void VDJ_cuts_insertion_dinuc_ntbias_model::CalculateAssignmentEntropies()
{
  //=========>to be implemented, for now leave it blank;
}
