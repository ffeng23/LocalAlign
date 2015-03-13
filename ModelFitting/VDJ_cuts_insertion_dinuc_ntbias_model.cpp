#include <cmath>
#include <iostream>
#include <stdlib.h>
#include "VDJ_cuts_insertion_dinuc_ntbias_model.hpp"

#include "VDJ_cuts_insertion_dinuc_ntbias_counter.hpp"
#include "VDJ_cuts_insertion_dinuc_ntbias_assigns.hpp"

//NOTE: to do: do we need to add the parameter such max_assignments, max_insertions, etc
//in the alignmentSettings (?) or some other way to specify these values!!!

VDJ_cuts_insertion_dinuc_ntbias_model::VDJ_cuts_insertion_dinuc_ntbias_model
(const GenomicV* _genV, const unsigned& _numV, 
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 const VDJ_cuts_insertion_dinuc_ntbias_model_params& vdj_mps
)
  :
  /*initilization list*/
  model_params(vdj_mps),
  max_assignments(vdj_mps.max_assignments), max_insertions(vdj_mps.max_insertions),
  max_V_deletions(vdj_mps.max_V_deletions), max_D_deletions(vdj_mps.max_D_deletions), max_J_deletions(vdj_mps.max_J_deletions),
  number_V_genes(vdj_mps.number_V_genes), number_D_genes(vdj_mps.number_D_genes), number_J_genes(vdj_mps.number_J_genes),
  max_V_n_alleles(vdj_mps.max_V_n_alleles), max_D_n_alleles(vdj_mps.max_D_n_alleles), max_J_n_alleles(vdj_mps.max_J_n_alleles),
  max_excess_V_deletions(vdj_mps.max_excess_V_deletions),  max_excess_D_deletions(vdj_mps.max_excess_D_deletions), max_excess_J_deletions(vdj_mps.max_excess_J_deletions),
  max_palindrome(vdj_mps.max_palindrome), 
  min_V_cut(vdj_mps.min_V_cut) /*cut will be set later*/, min_D_cut(vdj_mps.min_D_cut)/*later*/, min_J_cut(vdj_mps.min_J_cut)/*later*/,
  max_V_cut(vdj_mps.max_V_cut), max_D_cut(vdj_mps.max_D_cut),  max_J_cut(vdj_mps.max_J_cut), /*later*/

  /*entropies*/
  S_total(0), S_gene(0), S_V(0), S_DJ(0), S_D(0), S_J(0), S_insVD(0), S_insDJ(0),
  S_insVD_length(0), S_insDJ_length(0), S_insVD_nt(0), S_insDJ_nt(0), S_delV(0),
  S_delD(0), S_delJ(0), S_delDJ(0),
  
  negative_excess_deletions_max (vdj_mps.negative_excess_deletions_max), /*not sure why we set it up as 0 in here*/
  min_J_align_length(vdj_mps.min_J_align_length), min_J_assign_length(vdj_mps.min_J_assign_length), 
  min_V_length(vdj_mps.min_V_length)/*originally in matlab code is 15*/,
  high_error_region(vdj_mps.high_error_region) /*we PROBABLY will NOT use this one*/,
  use_no_D_match_seqs(vdj_mps.use_no_D_match_seqs),
  read_length(vdj_mps.read_length)/*the minmum read length has to be 101nts, is this good*/,

  /*model parameters*/
  PinsVD(), PinsDJ(),
  RnucleotideVD_per_nucleotideVD_5prime(), RnucleotideDJ_per_nucleotideDJ_3prime(),
  
  PcutV_given_V(), PcutJ_given_J(), PcutDlcutDr_given_D(),
  PV(), PDJ(), PVallele_given_gene(), PDallele_given_gene(), PJallele_given_gene(),
  
  Rerror_per_sequenced_nucleotide(0)
{
  //change 3/7/2015 by using vdj_..._model_params
  //need to do things to initialize the parameters and arrays
  /*min_V_cut=-1*max_palindrome;
  min_D_cut=-1*max_palindrome;
  min_J_cut=-1*max_palindrome;
  
  max_V_cut=max_V_deletions;
  max_D_cut=max_D_deletions;
  max_J_cut=max_J_deletions;
  
  //start intialize the various matrix
  max_gene_index(_genV,_numV);
  =max_gene_index(_genD,_numD);
  =max_gene_index(_genJ,_numJ);*/
  unsigned maxGeneIndex_V=vdj_mps.number_V_genes;
  unsigned maxGeneIndex_D=vdj_mps.number_D_genes;
  unsigned maxGeneIndex_J=vdj_mps.number_J_genes;

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
  //max_V_n_alleles = max_n_alleles(_genV, _numV);
  dim_size2[0]=maxGeneIndex_V;
  dim_size2[1]=max_V_n_alleles;
  PVallele_given_gene.initialize(2, dim_size2, 0.0);
  for(unsigned a=0;a<_numV;a++)
    {
      PVallele_given_gene(_genV[a].Get_GeneIndex(), _genV[a].Get_Allele()) = 
	1.0/_genV[a].Get_n_alleles();
    }       //end
  
  //D_genes = max([genD.gene_index]);
  //max_D_n_alleles = max_n_alleles(_genD, _numD);//...n_alleles]);
  //max_J_n_alleles =max_n_alleles(_genJ, _numJ);
  dim_size2[0]=maxGeneIndex_D;
  dim_size2[1]=max_D_n_alleles;
  PDallele_given_gene.initialize(2, dim_size2,0.0);// zeros(max_D_alleles , D_genes);
  for(unsigned a=0;a<_numD;a++)
    {
      PDallele_given_gene(_genD[a].Get_GeneIndex(), _genD[a].Get_Allele()) = 1.0/_genD[a].Get_n_alleles();
    }
  //J gene alleles
  dim_size2[0]=maxGeneIndex_J;
  dim_size2[1]=max_J_n_alleles;
  PJallele_given_gene.initialize(2, dim_size2, 0.0);
  for(unsigned a=0;a<_numV;a++)
    {
      PJallele_given_gene(_genJ[a].Get_GeneIndex(), _genJ[a].Get_Allele()) = 
	1.0/_genJ[a].Get_n_alleles();
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

  GeneTotal=sum(PJallele_given_gene,1);
  PJallele_given_gene.divide_by_dimension(GeneTotal,0);

  return true;
}

void  VDJ_cuts_insertion_dinuc_ntbias_model::InitializeCounter(Counter& _c) const
{
//we don't do anything
  //the actual job was done by constructor of the Counter;
  //return ;
}
//doing in_situ update
//basically in here we simply based on the counter update the necessary
//field. it is kind of doing statistics
void VDJ_cuts_insertion_dinuc_ntbias_model::GetModelFromCounter(Counter& _c) 
{
  //dynamic cast first
  VDJ_cuts_insertion_dinuc_ntbias_counter& vdj_model_counter =
    dynamic_cast<VDJ_cuts_insertion_dinuc_ntbias_counter&> (_c);
  
  //now for different fields, doing different things
  //nP*** only, no _given_
  PinsVD=vdj_model_counter.nPinsVD/sum_all(vdj_model_counter.nPinsVD);//normalization
  PinsDJ=vdj_model_counter.nPinsDJ/sum_all(vdj_model_counter.nPinsDJ);//normalization
  PV=vdj_model_counter.nPV/sum_all(vdj_model_counter.nPV);//normalization
  PDJ=vdj_model_counter.nPinsDJ/sum_all(vdj_model_counter.nPDJ);//normalization

  //nP*_given_*
  //doing conditional distribution
  Matrix<double> GeneTotal=sum(vdj_model_counter.nPcutV_given_V, 1);//sum the dim 1 (second dimension, since it starts at 0) that pointing to cutV, so to get geneTotal
  //the return matrix is of dimension _numV x 1 (vector).
  vdj_model_counter.nPcutV_given_V.divide_by_dimension(GeneTotal, 0);//divide along dimension 0 (first) to normalize
  this->PcutV_given_V.CopyValidElements(vdj_model_counter.nPcutV_given_V);

  GeneTotal=sum(vdj_model_counter.nPcutJ_given_J, 1);
  vdj_model_counter.nPcutJ_given_J.divide_by_dimension(GeneTotal, 0);
  this->PcutJ_given_J.CopyValidElements(vdj_model_counter.nPcutJ_given_J);

  Matrix<double> GeneTotal_temp=sum(vdj_model_counter.nPcutDlcutDr_given_D,2);
  GeneTotal=sum(GeneTotal_temp,1);
  vdj_model_counter.nPcutDlcutDr_given_D.divide_by_dimension(GeneTotal, 0);
  this->PcutDlcutDr_given_D.CopyValidElements(vdj_model_counter.nPcutDlcutDr_given_D);

  GeneTotal=sum(vdj_model_counter.nPVallele_given_gene,1);
  vdj_model_counter.nPVallele_given_gene.divide_by_dimension(GeneTotal, 0);
  this->PVallele_given_gene.CopyValidElements(vdj_model_counter.nPVallele_given_gene);

  GeneTotal=sum(vdj_model_counter.nPDallele_given_gene,1);
  vdj_model_counter.nPDallele_given_gene.divide_by_dimension(GeneTotal, 0);
  this->PDallele_given_gene.CopyValidElements(vdj_model_counter.nPDallele_given_gene);

  GeneTotal=sum(vdj_model_counter.nPJallele_given_gene, 1);
  vdj_model_counter.nPJallele_given_gene.divide_by_dimension(GeneTotal, 0);
  this->PJallele_given_gene.CopyValidElements(vdj_model_counter.nPJallele_given_gene);

  //now doing the R****** fields
  RnucleotideVD_per_nucleotideVD_5prime=
    vdj_model_counter.nMnucleotideVD/vdj_model_counter.nMnucleotideVD_5prime;
  RnucleotideDJ_per_nucleotideDJ_3prime=
    vdj_model_counter.nMnucleotideDJ/vdj_model_counter.nMnucleotideDJ_3prime;
  
  Rerror_per_sequenced_nucleotide=
    vdj_model_counter.nMerror/vdj_model_counter.nMsequenced_nucleotide;
  //M******** field
  //there is none for now!!???
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
  VDJ_cuts_insertion_dinuc_ntbias_counter& vdj_counter=
    dynamic_cast<VDJ_cuts_insertion_dinuc_ntbias_counter&>(_c);
  unsigned n_as=vdj_assigns.n_assignments;
  cout<<vdj_counter.nPinsVD.dim();
  if(n_as>0 && vdj_assigns.likelihood>0)
    {
      //first we simply create a temporary counter
      VDJ_cuts_insertion_dinuc_ntbias_counter delta_counter(this->model_params);
      
      delta_counter.logLikelihood=log(vdj_assigns.likelihood);

      delta_counter.N_processed=1;

      //normalize the probilitties of assignments
      vdj_assigns.proba=vdj_assigns.proba/vdj_assigns.likelihood;
      
      //now we want to update every relevant fields
      //first, nP**** fields.
      //--nPinsVD 
      if(!Update_nP_field(vdj_assigns.insVD,vdj_assigns.proba, delta_counter.nPinsVD, n_as ))
	{
	  cout<<"error on updating counter field insVD"<<endl;
	  exit(-1);
	}
      if(!Update_nP_field(vdj_assigns.insDJ,vdj_assigns.proba, delta_counter.nPinsDJ, n_as ))
	{
	  cout<<"error on updating counter field insDJ"<<endl;
	  exit(-1);
	}
      if(!Update_nP_field(vdj_assigns.cutV_given_V,vdj_assigns.proba, delta_counter.nPcutV_given_V, n_as ))
	{
	  cout<<"error on updating counter field cutV_given_V"<<endl;
	  exit(-1);
	}
      
      //cutJ_given_J.initialize(2, dim_size2, -1);
      if(!Update_nP_field(vdj_assigns.cutJ_given_J,vdj_assigns.proba, delta_counter.nPcutJ_given_J, n_as ))
	{
	  cout<<"error on updating counter field cutJ_given_J"<<endl;
	  exit(-1);
	}
      
      //cutDlcutDr_given_D.initialize(2, dim_size2, -1);;
      if(!Update_nP_field(vdj_assigns.cutDlcutDr_given_D,vdj_assigns.proba, delta_counter.nPcutDlcutDr_given_D, n_as ))
	{
	  cout<<"error on updating counter field cutDlcutDr_given_D"<<endl;
	  exit(-1);
	}
      
      //V.initialize(1, dim_size, -1);
      if(!Update_nP_field(vdj_assigns.V,vdj_assigns.proba, delta_counter.nPV, n_as ))
	{
	  cout<<"error on updating counter field V"<<endl;
	  exit(-1);
	}

      //DJ.initialize(2, dim_size2, -1); //Joint P(V, D, J gene choices)
      if(!Update_nP_field(vdj_assigns.DJ,vdj_assigns.proba, delta_counter.nPDJ, n_as ))
	{
	  cout<<"error on updating counter field DJ"<<endl;
	  exit(-1);
	}
      
      //Vallele_given_gene.initialize(2, dim_size2, -1); //Probabilities of alleles given gene for each gene
      if(!Update_nP_field(vdj_assigns.Vallele_given_gene, vdj_assigns.proba, delta_counter.nPVallele_given_gene, n_as ))
	{
	  cout<<"error on updating counter field Vallele_given_gene"<<endl;
	  exit(-1);
	}
      
      //Dallele_given_gene.initialize(2, dim_size2, -1);
      if(!Update_nP_field(vdj_assigns.Dallele_given_gene, vdj_assigns.proba, delta_counter.nPDallele_given_gene, n_as ))
	{
	  cout<<"error on updating counter field Dallele_given_gene"<<endl;
	  exit(-1);
	}

      //Jallele_given_gene.initialize(2, dim_size2, -1);
      if(!Update_nP_field(vdj_assigns.Jallele_given_gene, vdj_assigns.proba, delta_counter.nPJallele_given_gene, n_as ))
	{
	  cout<<"error on updating counter field Jallele_given_gene"<<endl;
	  exit(-1);
	}

      //VD_left_edge_dinucleotide.initialize(2, dim_size2, -1);// = zeros(4,4);
      if(!Update_nP_field(vdj_assigns.VD_left_edge_dinucleotide, vdj_assigns.proba, delta_counter.nPVD_left_edge_dinucleotide, n_as ))
	{
	  cout<<"error on updating counter field VD_left_edge_dinucleotide"<<endl;
	  exit(-1);
	}
      
      //VD_right_edge_dinucleotide.initialize(2, dim_size2, -1);// = zeros(4,4);
      if(!Update_nP_field(vdj_assigns.VD_right_edge_dinucleotide, vdj_assigns.proba, delta_counter.nPVD_right_edge_dinucleotide, n_as ))
	{
	  cout<<"error on updating counter field VD_right_edge_dinucleotide"<<endl;
	  exit(-1);
	}
      //DJ_left_edge_dinucleotide.initialize(2, dim_size2, -1);// = zeros(4,4);
      if(!Update_nP_field(vdj_assigns.DJ_left_edge_dinucleotide, vdj_assigns.proba, delta_counter.nPDJ_left_edge_dinucleotide, n_as ))
	{
	  cout<<"error on updating counter field DJ_left_edge_dinucleotide"<<endl;
	  exit(-1);
	}

      //DJ_right_edge_dinucleotide.initialize(2, dim_size2, -1);// = zeros(4,4);
      if(!Update_nP_field(vdj_assigns.DJ_right_edge_dinucleotide, vdj_assigns.proba, delta_counter.nPDJ_right_edge_dinucleotide, n_as ))
	{
	  cout<<"error on updating counter field DJ_right_edge_dinucleotide"<<endl;
	  exit(-1);
	}
   
      //pVmax_delV_V.initialize(2, dim_size2, -1);// =  zeros(model.max_palindrome + 1, model.max_V_deletions + 1, size(model.PV,1));
      if(!Update_nP_field(vdj_assigns.pVmax_delV_V, vdj_assigns.proba, delta_counter.nPpVmax_delV_V, n_as ))
	{
	  cout<<"error on updating counter field pVmax_delV_V"<<endl;
	  exit(-1);
	}
      //pJmax_delJ_J.initialize(2, dim_size2, -1);// = zeros(model.max_palindrome + 1, model.max_J_deletions + 1, size(model.PDJ,2));
      if(!Update_nP_field(vdj_assigns.pJmax_delJ_J, vdj_assigns.proba, delta_counter.nPpJmax_delJ_J, n_as ))
	{
	  cout<<"error on updating counter field pJmax_delJ_J"<<endl;
	  exit(-1);
	}

      //pDlmax_delDl_D.initialize(2, dim_size2, -1);// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));
      if(!Update_nP_field(vdj_assigns.pDlmax_delDl_D, vdj_assigns.proba, delta_counter.nPpDlmax_delDl_D, n_as ))
	{
	  cout<<"error on updating counter field pDlmax_delDl_D"<<endl;
	  exit(-1);
	}
      //pDrmax_delDr_D.initialize(2, dim_size2, -1);// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));
      if(!Update_nP_field(vdj_assigns.pDrmax_delDr_D, vdj_assigns.proba, delta_counter.nPpDrmax_delDr_D, n_as ))
	{
	  cout<<"error on updating counter field pDrmax_delDr_D"<<endl;
	  exit(-1);
	}

      //VDJ.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1), size(model.PDJ,1), size(model.PDJ,2));
      if(!Update_nP_field(vdj_assigns.VDJ, vdj_assigns.proba, delta_counter.nPVDJ, n_as ))
	{
	  cout<<"error on updating counter field VDJ"<<endl;
	  exit(-1);
	}
      //pVdelV.initialize(2, dim_size2, -1);// = zeros( model.max_palindrome + 1, model.max_V_deletions + 1);
      if(!Update_nP_field(vdj_assigns.pVdelV, vdj_assigns.proba, delta_counter.nPpVdelV, n_as ))
	{
	  cout<<"error on updating counter field pVdelV"<<endl;
	  exit(-1);
	}
      //pDldelDl.initialize(2, dim_size2, -1);// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
      if(!Update_nP_field(vdj_assigns.pDldelDl, vdj_assigns.proba, delta_counter.nPpDldelDl, n_as ))
	{
	  cout<<"error on updating counter field pDldelDl"<<endl;
	  exit(-1);
	}

      //pDrdelDr.initialize(2, dim_size2, -1);// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
      if(!Update_nP_field(vdj_assigns.pDrdelDr, vdj_assigns.proba, delta_counter.nPpDrdelDr, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}

      //pJdelJ.initialize(2, dim_size2, -1);// = zeros( model.max_palindrome + 1, model.max_J_deletions + 1);
      if(!Update_nP_field(vdj_assigns.pJdelJ, vdj_assigns.proba, delta_counter.nPpJdelJ, n_as ))
	{
	  cout<<"error on updating counter field pJdelJ"<<endl;
	  exit(-1);
	}

      //VcutV.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1),model.max_palindrome + model.max_V_deletions + 1);
      if(!Update_nP_field(vdj_assigns.VcutV, vdj_assigns.proba, delta_counter.nPVcutV, n_as ))
	{
	  cout<<"error on updating counter field VcutV"<<endl;
	  exit(-1);
	}

      //DcutDl.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
      if(!Update_nP_field(vdj_assigns.DcutDl, vdj_assigns.proba, delta_counter.nPDcutDl, n_as ))
	{
	  cout<<"error on updating counter field DcutDl"<<endl;
	  exit(-1);
	}
      //DcutDr.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
if(!Update_nP_field(vdj_assigns.DcutDr, vdj_assigns.proba, delta_counter.nPDcutDr, n_as ))
	{
	  cout<<"error on updating counter field DcutDr"<<endl;
	  exit(-1);
	}
//JcutJ.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_J_deletions + 1);
if(!Update_nP_field(vdj_assigns.JcutJ, vdj_assigns.proba, delta_counter.nPJcutJ, n_as ))
	{
	  cout<<"error on updating counter field JcutJ"<<endl;
	  exit(-1);
	}
//DcutV.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_V_deletions + 1);
if(!Update_nP_field(vdj_assigns.DcutV, vdj_assigns.proba, delta_counter.nPDcutV, n_as ))
	{
	  cout<<"error on updating counter field DcutV"<<endl;
	  exit(-1);
	}
//DcutJ.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_J_deletions + 1);
if(!Update_nP_field(vdj_assigns.DcutJ, vdj_assigns.proba, delta_counter.nPDcutJ, n_as ))
	{
	  cout<<"error on updating counter field DcutJ"<<endl;
	  exit(-1);
	}
//VcutJ.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1),model.max_palindrome + model.max_J_deletions + 1);
if(!Update_nP_field(vdj_assigns.VcutJ, vdj_assigns.proba, delta_counter.nPVcutJ, n_as ))
	{
	  cout<<"error on updating counter field VcutJ"<<endl;
	  exit(-1);
	}
//VcutDl.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);
if(!Update_nP_field(vdj_assigns.VcutDl, vdj_assigns.proba, delta_counter.nPVcutDl, n_as ))
	{
	  cout<<"error on updating counter field VcutDl"<<endl;
	  exit(-1);
	}
//VcutDr.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);
if(!Update_nP_field(vdj_assigns.VcutDr, vdj_assigns.proba, delta_counter.nPVcutDr, n_as ))
	{
	  cout<<"error on updating counter field VcutDr"<<endl;
	  exit(-1);
	}

//JcutDl.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
if(!Update_nP_field(vdj_assigns.JcutDl, vdj_assigns.proba, delta_counter.nPJcutDl, n_as ))
	{
	  cout<<"error on updating counter field JcutDl"<<endl;
	  exit(-1);
	}
//JcutDr.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
if(!Update_nP_field(vdj_assigns.JcutDr, vdj_assigns.proba, delta_counter.nPJcutDr, n_as ))
	{
	  cout<<"error on updating counter field JcutDr"<<endl;
	  exit(-1);
	}
//JcutV.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_V_deletions + 1);
if(!Update_nP_field(vdj_assigns.JcutV, vdj_assigns.proba, delta_counter.nPJcutV, n_as ))
	{
	  cout<<"error on updating counter field JcutV"<<endl;
	  exit(-1);
	}
//insVDcutV.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);
if(!Update_nP_field(vdj_assigns.insVDcutV, vdj_assigns.proba, delta_counter.nPinsVDcutV, n_as ))
	{
	  cout<<"error on updating counter field insVDcutV"<<endl;
	  exit(-1);
	}
//insDJcutV.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);
if(!Update_nP_field(vdj_assigns.insDJcutV, vdj_assigns.proba, delta_counter.nPinsDJcutV, n_as ))
	{
	  cout<<"error on updating counter field insDJcutV"<<endl;
	  exit(-1);
	}
//insVDcutDl.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
if(!Update_nP_field(vdj_assigns.insVDcutDl, vdj_assigns.proba, delta_counter.nPinsVDcutDl, n_as ))
	{
	  cout<<"error on updating counter field insVDcutDl"<<endl;
	  exit(-1);
	}
//insDJcutDl.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
if(!Update_nP_field(vdj_assigns.insDJcutDl, vdj_assigns.proba, delta_counter.nPinsDJcutDl, n_as ))
	{
	  cout<<"error on updating counter field insDJcutDl"<<endl;
	  exit(-1);
	}

//insVDcutDr.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
if(!Update_nP_field(vdj_assigns.insVDcutDr, vdj_assigns.proba, delta_counter.nPinsVDcutDr, n_as ))
	{
	  cout<<"error on updating counter field insVDcutDr"<<endl;
	  exit(-1);
	}
//insDJcutDr.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
if(!Update_nP_field(vdj_assigns.insDJcutDr, vdj_assigns.proba, delta_counter.nPinsDJcutDr, n_as ))
	{
	  cout<<"error on updating counter field insDJcutDr"<<endl;
	  exit(-1);
	}

//insVDcutJ.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);
if(!Update_nP_field(vdj_assigns.insVDcutJ, vdj_assigns.proba, delta_counter.nPinsVDcutJ, n_as ))
	{
	  cout<<"error on updating counter field insVDcutJ"<<endl;
	  exit(-1);
	}
//insDJcutJ.initialize(2, dim_size2, -1);// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);
if(!Update_nP_field(vdj_assigns.insDJcutJ, vdj_assigns.proba, delta_counter.nPinsDJcutJ, n_as ))
	{
	  cout<<"error on updating counter field insDJcutJ"<<endl;
	  exit(-1);
	}
  
//V_align_length.initialize(1, dim_size, -1);// = zeros(max_V_length + 1,1);
if(!Update_nP_field(vdj_assigns.V_align_length, vdj_assigns.proba, delta_counter.nPV_align_length, n_as ))
	{
	  cout<<"error on updating counter field V_align_length"<<endl;
	  exit(-1);
	}
//D_align_length.initialize(1, dim_size, -1);// = zeros(max_D_length + 1,1);
if(!Update_nP_field(vdj_assigns.D_align_length, vdj_assigns.proba, delta_counter.nPD_align_length, n_as ))
	{
	  cout<<"error on updating counter field D_align_length"<<endl;
	  exit(-1);
	}
//J_align_length.initialize(1, dim_size, -1);// = zeros(max_J_length + 1,1);
if(!Update_nP_field(vdj_assigns.J_align_length, vdj_assigns.proba, delta_counter.nPJ_align_length, n_as ))
	{
	  cout<<"error on updating counter field J_align_length"<<endl;
	  exit(-1);
	}

//JJ_align_length.initialize(2, dim_size2, -1);// = zeros(size(model.PDJ,2), 1 + max_J_length);
if(!Update_nP_field(vdj_assigns.JJ_align_length, vdj_assigns.proba, delta_counter.nPJJ_align_length, n_as ))
	{
	  cout<<"error on updating counter field JJ_align_length"<<endl;
	  exit(-1);
	}
//   VV_align_length.initialize(2, dim_size2, -1);// = zeros(size(model.PV,1), 1 + max_V_length);
if(!Update_nP_field(vdj_assigns.VV_align_length, vdj_assigns.proba, delta_counter.nPVV_align_length, n_as ))
	{
	  cout<<"error on updating counter field VV_align_length"<<endl;
	  exit(-1);
	}


//insDJ_D_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_D_length + 1);
if(!Update_nP_field(vdj_assigns.insDJ_D_align_length, vdj_assigns.proba, delta_counter.nPinsDJ_D_align_length, n_as ))
	{
	  cout<<"error on updating counter field insDJ_D_align_length"<<endl;
	  exit(-1);
	}
//insVD_D_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_D_length + 1);
if(!Update_nP_field(vdj_assigns.insVD_D_align_length, vdj_assigns.proba, delta_counter.nPinsVD_D_align_length, n_as ))
	{
	  cout<<"error on updating counter field insVD_D_align_length"<<endl;
	  exit(-1);
	}

//   insDJ_J_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_J_length + 1);
if(!Update_nP_field(vdj_assigns.insDJ_J_align_length, vdj_assigns.proba, delta_counter.nPinsDJ_J_align_length, n_as ))
	{
	  cout<<"error on updating counter field insDJ_J_align_length"<<endl;
	  exit(-1);
	}

//insVD_J_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_J_length + 1);
if(!Update_nP_field(vdj_assigns.insVD_J_align_length, vdj_assigns.proba, delta_counter.nPinsVD_J_align_length, n_as ))
	{
	  cout<<"error on updating counter field insVD_J_align_length"<<endl;
	  exit(-1);
	}

//  insDJ_V_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_V_length + 1);
if(!Update_nP_field(vdj_assigns.insDJ_V_align_length, vdj_assigns.proba, delta_counter.nPinsDJ_V_align_length, n_as ))
	{
	  cout<<"error on updating counter field insDJ_V_align_length"<<endl;
	  exit(-1);
	}
//insVD_V_align_length.initialize(2, dim_size2, -1);// = zeros(model.max_insertions +1, max_V_length + 1);
if(!Update_nP_field(vdj_assigns.insVD_V_align_length, vdj_assigns.proba, delta_counter.nPinsVD_V_align_length, n_as ))
	{
	  cout<<"error on updating counter field insVD_V_align_length"<<endl;
	  exit(-1);
	}

//Dallele_D_align_length.initialize(2, dim_size2, -1);// = zeros(3, max_D_length + 1);
if(!Update_nP_field(vdj_assigns.Dallele_D_align_length, vdj_assigns.proba, delta_counter.nPDallele_D_align_length, n_as ))
	{
	  cout<<"error on updating counter field Dallele_D_align_length"<<endl;
	  exit(-1);
	}


//delVinsVD.initialize(2, dim_size2, -1);// = zeros(model.max_V_deletions +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.delVinsVD, vdj_assigns.proba, delta_counter.nPdelVinsVD, n_as ))
	{
	  cout<<"error on updating counter field delVinsVD"<<endl;
	  exit(-1);
	}
//delVinsDJ.initialize(2, dim_size2, -1);// = zeros(model.max_V_deletions +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.delVinsDJ, vdj_assigns.proba, delta_counter.nPdelVinsDJ, n_as ))
	{
	  cout<<"error on updating counter field delVinsDJ"<<endl;
	  exit(-1);
	}
//delVdelDl.initialize(2, dim_size2, -1);// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.delVdelDl, vdj_assigns.proba, delta_counter.nPdelVdelDl, n_as ))
	{
	  cout<<"error on updating counter field delVdelDl"<<endl;
	  exit(-1);
	}

//delVdelDr.initialize(2, dim_size2, -1);// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.delVdelDr, vdj_assigns.proba, delta_counter.nPdelVdelDr, n_as ))
	{
	  cout<<"error on updating counter field delVdelDr"<<endl;
	  exit(-1);
	}
//delVdelJ.initialize(2, dim_size2, -1);// = zeros(model.max_V_deletions +1, model.max_J_deletions +1);
if(!Update_nP_field(vdj_assigns.delVdelJ, vdj_assigns.proba, delta_counter.nPdelVdelJ, n_as ))
	{
	  cout<<"error on updating counter field delVdelJ"<<endl;
	  exit(-1);
	}

//delJinsVD.initialize(2, dim_size2, -1);// = zeros(model.max_J_deletions +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.delJinsVD, vdj_assigns.proba, delta_counter.nPdelJinsVD, n_as ))
	{
	  cout<<"error on updating counter field delJinsVD"<<endl;
	  exit(-1);
	}
//delJinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_J_deletions +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.delJinsDJ, vdj_assigns.proba, delta_counter.nPdelJinsDJ, n_as ))
	{
	  cout<<"error on updating counter field delJinsDJ"<<endl;
	  exit(-1);
	}
//delJdelDl.initialize(2, dim_size2, -1);//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.delJdelDl, vdj_assigns.proba, delta_counter.nPdelJdelDl, n_as ))
	{
	  cout<<"error on updating counter field delJdelDl"<<endl;
	  exit(-1);
	}
//   delJdelDr.initialize(2, dim_size2, -1);//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.delJdelDr, vdj_assigns.proba, delta_counter.nPdelJdelDr, n_as ))
	{
	  cout<<"error on updating counter field delJdelDr"<<endl;
	  exit(-1);
	}
  
//delDlinsVD.initialize(2, dim_size2, -1);//zeros(model.max_D_deletions +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.delDlinsVD, vdj_assigns.proba, delta_counter.nPdelDlinsVD, n_as ))
	{
	  cout<<"error on updating counter field delDlinsVD"<<endl;
	  exit(-1);
	}
//delDlinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_D_deletions +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.delDlinsDJ, vdj_assigns.proba, delta_counter.nPdelDlinsDJ, n_as ))
	{
	  cout<<"error on updating counter field delDlinsDJ"<<endl;
	  exit(-1);
	}
//delDldelDr.initialize(2, dim_size2, -1);//zeros(model.max_D_deletions +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.delDldelDr, vdj_assigns.proba, delta_counter.nPdelDldelDr, n_as ))
	{
	  cout<<"error on updating counter field delDldelDr"<<endl;
	  exit(-1);
	}
  
//delDrinsVD.initialize(2, dim_size2, -1);//zeros(model.max_D_deletions +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.delDrinsVD, vdj_assigns.proba, delta_counter.nPdelDrinsVD, n_as ))
	{
	  cout<<"error on updating counter field delDrinsVD"<<endl;
	  exit(-1);
	}
//delDrinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_D_deletions +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.delDrinsDJ, vdj_assigns.proba, delta_counter.nPdelDrinsDJ, n_as ))
	{
	  cout<<"error on updating counter field delDrinsDJ"<<endl;
	  exit(-1);
	}
  
//insVDinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_insertions +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.insVDinsDJ, vdj_assigns.proba, delta_counter.nPinsVDinsDJ, n_as ))
	{
	  cout<<"error on updating counter field insVDinsDJ"<<endl;
	  exit(-1);
	}
  
//VdelV.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_V_deletions+1 );
if(!Update_nP_field(vdj_assigns.VdelV, vdj_assigns.proba, delta_counter.nPVdelV, n_as ))
	{
	  cout<<"error on updating counter field VdelV"<<endl;
	  exit(-1);
	}
//DdelDl.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
if(!Update_nP_field(vdj_assigns.DdelDl, vdj_assigns.proba, delta_counter.nPDdelDl, n_as ))
	{
	  cout<<"error on updating counter field DdelDl"<<endl;
	  exit(-1);
	}

//DdelDr.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
if(!Update_nP_field(vdj_assigns.DdelDr, vdj_assigns.proba, delta_counter.nPDdelDr, n_as ))
	{
	  cout<<"error on updating counter field DdelDr"<<endl;
	  exit(-1);
	}

//JdelJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_J_deletions+1 );
if(!Update_nP_field(vdj_assigns.JdelJ, vdj_assigns.proba, delta_counter.nPJdelJ, n_as ))
	{
	  cout<<"error on updating counter field JdelJ"<<endl;
	  exit(-1);
	}
  
//VinsVD.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_insertions+1 );
if(!Update_nP_field(vdj_assigns.VinsVD, vdj_assigns.proba, delta_counter.nPVinsVD, n_as ))
	{
	  cout<<"error on updating counter field VinsVD"<<endl;
	  exit(-1);
	}
//DinsVD.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_insertions+1 );
if(!Update_nP_field(vdj_assigns.DinsVD, vdj_assigns.proba, delta_counter.nPDinsVD, n_as ))
	{
	  cout<<"error on updating counter field DinsVD"<<endl;
	  exit(-1);
	}
//DinsDJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_insertions+1 );
if(!Update_nP_field(vdj_assigns.DinsDJ, vdj_assigns.proba, delta_counter.nPDinsDJ, n_as ))
	{
	  cout<<"error on updating counter field DinsDJ"<<endl;
	  exit(-1);
	}
//JinsDJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_insertions+1 );
if(!Update_nP_field(vdj_assigns.JinsDJ, vdj_assigns.proba, delta_counter.nPJinsDJ, n_as ))
	{
	  cout<<"error on updating counter field JinsDJ"<<endl;
	  exit(-1);
	}
  
//VdelDl.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_D_deletions+1 );
if(!Update_nP_field(vdj_assigns.VdelDl, vdj_assigns.proba, delta_counter.nPVdelDl, n_as ))
	{
	  cout<<"error on updating counter field VdelDl"<<endl;
	  exit(-1);
	}
//VdelDr.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_D_deletions+1 );
if(!Update_nP_field(vdj_assigns.VdelDr, vdj_assigns.proba, delta_counter.nPVdelDr, n_as ))
	{
	  cout<<"error on updating counter field VdelDr"<<endl;
	  exit(-1);
	}
//VdelJ.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_J_deletions+1 );
if(!Update_nP_field(vdj_assigns.VdelJ, vdj_assigns.proba, delta_counter.nPVdelJ, n_as ))
	{
	  cout<<"error on updating counter field VdelJ"<<endl;
	  exit(-1);
	}
  
//JdelV.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_V_deletions+1 );
if(!Update_nP_field(vdj_assigns.JdelV, vdj_assigns.proba, delta_counter.nPJdelV, n_as ))
	{
	  cout<<"error on updating counter field JdelV"<<endl;
	  exit(-1);
	}
//JdelDl.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
if(!Update_nP_field(vdj_assigns.JdelDl, vdj_assigns.proba, delta_counter.nPJdelDl, n_as ))
	{
	  cout<<"error on updating counter field JdelDl"<<endl;
	  exit(-1);
	}
//JdelDr.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
if(!Update_nP_field(vdj_assigns.JdelDr, vdj_assigns.proba, delta_counter.nPJdelDr, n_as ))
	{
	  cout<<"error on updating counter field JdelDr"<<endl;
	  exit(-1);
	}
  
//DdelV.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_V_deletions+1 );
if(!Update_nP_field(vdj_assigns.DdelV, vdj_assigns.proba, delta_counter.nPDdelV, n_as ))
	{
	  cout<<"error on updating counter field DdelV"<<endl;
	  exit(-1);
	}
//DdelJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_J_deletions+1 );
if(!Update_nP_field(vdj_assigns.DdelJ, vdj_assigns.proba, delta_counter.nPDdelJ, n_as ))
	{
	  cout<<"error on updating counter field DdelJ"<<endl;
	  exit(-1);
	}
  
//VinsDJ.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_insertions+1 );
if(!Update_nP_field(vdj_assigns.VinsDJ, vdj_assigns.proba, delta_counter.nPVinsDJ, n_as ))
	{
	  cout<<"error on updating counter field VinsDJ"<<endl;
	  exit(-1);
	}
//JinsVD.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_insertions+1 );
if(!Update_nP_field(vdj_assigns.JinsVD, vdj_assigns.proba, delta_counter.nPJinsVD, n_as ))
	{
	  cout<<"error on updating counter field JinsVD"<<endl;
	  exit(-1);
	}
  
//pVinsVD.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.pVinsVD, vdj_assigns.proba, delta_counter.nPpVinsVD, n_as ))
	{
	  cout<<"error on updating counter field pVinsVD"<<endl;
	  exit(-1);
	}
//pVinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.pVinsDJ, vdj_assigns.proba, delta_counter.nPpVinsDJ, n_as ))
	{
	  cout<<"error on updating counter field pVinsDJ"<<endl;
	  exit(-1);
	}
//pVdelDl.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.pVdelDl, vdj_assigns.proba, delta_counter.nPpVdelDl, n_as ))
	{
	  cout<<"error on updating counter field pVdelDl"<<endl;
	  exit(-1);
	}
//pVdelDr.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.pVdelDr, vdj_assigns.proba, delta_counter.nPpVdelDr, n_as ))
	{
	  cout<<"error on updating counter field pVdelDr"<<endl;
	  exit(-1);
	}
//pVdelJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_J_deletions +1);
if(!Update_nP_field(vdj_assigns.pVdelJ, vdj_assigns.proba, delta_counter.nPpVdelJ, n_as ))
	{
	  cout<<"error on updating counter field pVdelJ"<<endl;
	  exit(-1);
	}
//VpV.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.VpV, vdj_assigns.proba, delta_counter.nPVpV, n_as ))
	{
	  cout<<"error on updating counter field VpV"<<endl;
	  exit(-1);
	}
//JpV.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.JpV, vdj_assigns.proba, delta_counter.nPJpV, n_as ))
	{
	  cout<<"error on updating counter field JpV"<<endl;
	  exit(-1);
	}
//DpV.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.DpV, vdj_assigns.proba, delta_counter.nPDpV, n_as ))
	{
	  cout<<"error on updating counter field DpV"<<endl;
	  exit(-1);
	}
  
//pJinsVD.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.pJinsVD, vdj_assigns.proba, delta_counter.nPpJinsVD, n_as ))
	{
	  cout<<"error on updating counter field pJinsVD"<<endl;
	  exit(-1);
	}

//pJinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.pJinsDJ, vdj_assigns.proba, delta_counter.nPpJinsDJ, n_as ))
	{
	  cout<<"error on updating counter field pJinsDJ"<<endl;
	  exit(-1);
	}
//pJdelV.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_V_deletions +1);
if(!Update_nP_field(vdj_assigns.pJdelV, vdj_assigns.proba, delta_counter.nPpJdelV, n_as ))
	{
	  cout<<"error on updating counter field pJdelV"<<endl;
	  exit(-1);
	}
//pJdelDl.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.pJdelDl, vdj_assigns.proba, delta_counter.nPpJdelDl, n_as ))
	{
	  cout<<"error on updating counter field pJdelDl"<<endl;
	  exit(-1);
	}
//pJdelDr.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.pJdelDr, vdj_assigns.proba, delta_counter.nPpJdelDr, n_as ))
	{
	  cout<<"error on updating counter field pJdelDr"<<endl;
	  exit(-1);
	}
//VpJ.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.VpJ, vdj_assigns.proba, delta_counter.nPVpJ, n_as ))
	{
	  cout<<"error on updating counter field VpJ"<<endl;
	  exit(-1);
	}
//JpJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.JpJ, vdj_assigns.proba, delta_counter.nPJpJ, n_as ))
	{
	  cout<<"error on updating counter field JpJ"<<endl;
	  exit(-1);
	}
//DpJ.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.DpJ, vdj_assigns.proba, delta_counter.nPDpJ, n_as ))
	{
	  cout<<"error on updating counter field DpJ"<<endl;
	  exit(-1);
	}
  
//pDlinsVD.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.pDlinsVD, vdj_assigns.proba, delta_counter.nPpDlinsVD, n_as ))
	{
	  cout<<"error on updating counter field pDlinsVD"<<endl;
	  exit(-1);
	}
//pDlinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.pDlinsDJ, vdj_assigns.proba, delta_counter.nPpDlinsDJ, n_as ))
	{
	  cout<<"error on updating counter field pDlinsDJ"<<endl;
	  exit(-1);
	}
//pDldelV.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_V_deletions +1);
if(!Update_nP_field(vdj_assigns.pDldelV, vdj_assigns.proba, delta_counter.nPpDldelV, n_as ))
	{
	  cout<<"error on updating counter field pDldelV"<<endl;
	  exit(-1);
	}

//pDldelJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_J_deletions +1);
if(!Update_nP_field(vdj_assigns.pDldelJ, vdj_assigns.proba, delta_counter.nPpDldelJ, n_as ))
	{
	  cout<<"error on updating counter field pDldelJ"<<endl;
	  exit(-1);
	}
//pDldelDr.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.pDldelDr, vdj_assigns.proba, delta_counter.nPpDldelDr, n_as ))
	{
	  cout<<"error on updating counter field pDldelDr"<<endl;
	  exit(-1);
	}
//VpDl.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.VpDl, vdj_assigns.proba, delta_counter.nPVpDl, n_as ))
	{
	  cout<<"error on updating counter field VpDl"<<endl;
	  exit(-1);
	}
//JpDl.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.JpDl, vdj_assigns.proba, delta_counter.nPJpDl, n_as ))
	{
	  cout<<"error on updating counter field JpDl"<<endl;
	  exit(-1);
	}
//DpDl.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.DpDl, vdj_assigns.proba, delta_counter.nPDpDl, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}
  
//pDrinsVD.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.pDrinsVD, vdj_assigns.proba, delta_counter.nPpDrinsVD, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}
//pDrinsDJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_insertions +1);
if(!Update_nP_field(vdj_assigns.pDrinsDJ, vdj_assigns.proba, delta_counter.nPpDrinsDJ, n_as ))
	{
	  cout<<"error on updating counter field pDrinsDJ"<<endl;
	  exit(-1);
	}

//pDrdelV.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_V_deletions +1);
if(!Update_nP_field(vdj_assigns.pDrdelV, vdj_assigns.proba, delta_counter.nPpDrdelV, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}
//pDrdelJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_J_deletions +1);
if(!Update_nP_field(vdj_assigns.pDrdelJ, vdj_assigns.proba, delta_counter.nPpDrdelJ, n_as ))
	{
	  cout<<"error on updating counter field pDrdelJ"<<endl;
	  exit(-1);
	}
//pDrdelDl.initialize(2, dim_size2, -1);//zeros(model.max_palindrome +1, model.max_D_deletions +1);
if(!Update_nP_field(vdj_assigns.pDrdelDl, vdj_assigns.proba, delta_counter.nPpDrdelDl, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}
//VpDr.initialize(2, dim_size2, -1);//zeros(size(model.PV,1), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.VpDr, vdj_assigns.proba, delta_counter.nPVpDr, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}
//JpDr.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,2), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.JpDr, vdj_assigns.proba, delta_counter.nPJpDr, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}
//DpDr.initialize(2, dim_size2, -1);//zeros(size(model.PDJ,1), model.max_palindrome+1 );
if(!Update_nP_field(vdj_assigns.DpDr, vdj_assigns.proba, delta_counter.nPDpDr, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}
  
//pVpDl.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
if(!Update_nP_field(vdj_assigns.pVpDl, vdj_assigns.proba, delta_counter.nPpVpDl, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}
//pVpDr.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
if(!Update_nP_field(vdj_assigns.pVpDr, vdj_assigns.proba, delta_counter.nPpVpDr, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}
//pVpJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
if(!Update_nP_field(vdj_assigns.pVpJ, vdj_assigns.proba, delta_counter.nPpVpJ, n_as ))
	{
	  cout<<"error on updating counter field pVpJ"<<endl;
	  exit(-1);
	}
//pDlpDr.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
if(!Update_nP_field(vdj_assigns.pDlpDr, vdj_assigns.proba, delta_counter.nPpDlpDr, n_as ))
	{
	  cout<<"error on updating counter field pDrdelDr"<<endl;
	  exit(-1);
	}
//pDlpJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
if(!Update_nP_field(vdj_assigns.pDlpJ, vdj_assigns.proba, delta_counter.nPpDlpJ, n_as ))
	{
	  cout<<"error on updating counter field pDlpJ"<<endl;
	  exit(-1);
	}
//pDrpJ.initialize(2, dim_size2, -1);//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
if(!Update_nP_field(vdj_assigns.pDrpJ, vdj_assigns.proba, delta_counter.nPpDrpJ, n_as ))
	{
	  cout<<"error on updating counter field pDrpJ"<<endl;
	  exit(-1);
	}

//now we are doing nM*******fields
//mononucleotideVD.initialize(2, dim_size2, 0.0);// zeros(4,1);
 if(!Update_nM_field(vdj_assigns.mononucleotideVD, vdj_assigns.proba, delta_counter.nMmononucleotideVD, n_as))
   {
     cout<<"error on updating counter field mononucleotideVD"<<endl;
	  exit(-1);
   }
 //mononucleotideDJ.initialize(2, dim_size2, 0.0);// = zeros(4,1);
if(!Update_nM_field(vdj_assigns.mononucleotideDJ, vdj_assigns.proba, delta_counter.nMmononucleotideDJ, n_as))
   {
     cout<<"error on updating counter field mononucleotideDJ"<<endl;
	  exit(-1);
   }

//insertionVD.initialize(2, dim_size2, 0.0);// = zeros(4,1);what this is?? is this insertion or nucleotide dist'n among insertion??
if(!Update_nM_field(vdj_assigns.insertionVD, vdj_assigns.proba, delta_counter.nMinsertionVD, n_as))
   {
     cout<<"error on updating counter field insertionVD"<<endl;
	  exit(-1);
   }
//insertionDJ.initialize(2, dim_size2, 0.0);// = zeros(4,1);
if(!Update_nM_field(vdj_assigns.insertionDJ, vdj_assigns.proba, delta_counter.nMinsertionDJ, n_as))
   {
     cout<<"error on updating counter field insertionDJ"<<endl;
	  exit(-1);
   }

//trinucleotideVD.initialize(4, dim_size4, 0.0);// = zeros(4,4,4);
if(!Update_nM_field(vdj_assigns.trinucleotideVD, vdj_assigns.proba, delta_counter.nMtrinucleotideVD, n_as))
   {
     cout<<"error on updating counter field trinucleotideVD"<<endl;
	  exit(-1);
   }
//trinucleotideDJ.initialize(4, dim_size4, 0.0);// = zeros(4,4,4);
if(!Update_nM_field(vdj_assigns.trinucleotideDJ, vdj_assigns.proba, delta_counter.nMtrinucleotideDJ, n_as))
   {
     cout<<"error on updating counter field trinucleotideDJ"<<endl;
	  exit(-1);
   }

//VV_err_pos.initialize(3, dim_size3, 0.0);// = zeros(size(model.PV,1), max_V_length + model.max_V_deletions);
if(!Update_nM_field(vdj_assigns.VV_err_pos, vdj_assigns.proba, delta_counter.nMVV_err_pos, n_as))
   {
     cout<<"error on updating counter field VV_err_pos"<<endl;
	  exit(-1);
   }

//JJ_err_pos.initialize(3, dim_size3, 0.0);// = zeros(size(model.PDJ,2), max_J_length);
if(!Update_nM_field(vdj_assigns.JJ_err_pos, vdj_assigns.proba, delta_counter.nMJJ_err_pos, n_as))
   {
     cout<<"error on updating counter field JJ_err_pos"<<endl;
	  exit(-1);
   }
  
//zeroD.initialize(1, dim_size, 0.0) ;
if(!Update_nM_field(vdj_assigns.zeroD, vdj_assigns.proba, delta_counter.nMzeroD, n_as))
   {
     cout<<"error on updating counter field zeroD"<<endl;
	  exit(-1);
   }

// nucleotideVD.initialize(3, dim_size3, 0.0); nucleotideVD_5prime.initialize(3, dim_size3, 0.0);//nucleotide distr's
if(!Update_nM_field(vdj_assigns.nucleotideVD, vdj_assigns.proba, delta_counter.nMnucleotideVD, n_as))
   {
     cout<<"error on updating counter field nucleotideVD"<<endl;
	  exit(-1);
   }
//nucleotideDJ.initialize(3, dim_size3, 0.0); nucleotideDJ_3prime.initialize(3, dim_size3, 0.0);
if(!Update_nM_field(vdj_assigns.nucleotideDJ, vdj_assigns.proba, delta_counter.nMnucleotideDJ, n_as))
   {
     cout<<"error on updating counter field nucleotideDJ"<<endl;
	  exit(-1);
   }
//error.initialize(1, dim_size,0.0); sequenced_nucleotide.initialize(1,dim_size, 0.0);//error rate, will be 
  if(!Update_nM_field(vdj_assigns.error, vdj_assigns.proba, delta_counter.nMerror, n_as))
   {
     cout<<"error on updating counter field error"<<endl;
	  exit(-1);
   }
  
  //error_vs_position.initialize(2,dim_size2,0.0);//!!this is nM zeros(model.read_length,1);
if(!Update_nM_field(vdj_assigns.error_vs_position, vdj_assigns.proba, delta_counter.nMerror_vs_position, n_as))
   {
     cout<<"error on updating counter field error_vs_position"<<endl;
	  exit(-1);
   }
//coverage.initialize(2,dim_size2,0.0);//zeros(model.read_length,1); 
 if(!Update_nM_field(vdj_assigns.coverage, vdj_assigns.proba, delta_counter.nMcoverage, n_as))
   {
     cout<<"error on updating counter field mononucleotideVD"<<endl;
     exit(-1);
   }
    }//end of outer loop for check the valid assignmnet
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
    (Matrix<unsigned>& _index_field_a, Matrix<double>& _prob_field_a, 
     Matrix<double>& _field_c, const unsigned& _num_valid_assignments) const
  {
      bool flag=true;
      //first check for the correct input
      //_index_field_a is max_assignments x index1 x index2...x index_n
      //_prob_field_a is max_assignments x 1
      //_fields_c =index1 x index2 ... x index_n
      
      if(_index_field_a.size(0)<_num_valid_assignments)
	{
	  cout<<"Error: pass in not number of valid assignments or incorrect sized index fields, please check"<<endl;
	  return false;
	}
      
      if(_prob_field_a.dim()!=1)
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
	 (_index_field_a.size(1)!=_field_c.dim())
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
      //unsigned blockSize=_field_c.nTotal();//total number of elements in the counter field
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
	      for(unsigned k=0;k<_index_field_a.size(1);k++)
		{
		  position_index_field_a[1]=k;
		  //figure out the indices of the specific entry
		  position_field_c[k]=_index_field_a(position_index_field_a, _index_field_a.dim()); 
		  
		}
	      //now we got position, let's up 
	      _field_c(position_field_c,_field_c.dim() )+=_prob_field_a(k_all_good_assigns[i])/totalProb;
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
    (Matrix<double> _count_field_a, Matrix<double> _prob_field_a, 
     Matrix<double> _field_c, const unsigned& _num_valid_assignments) const
  {
    bool flag=true;
    //first check for the validity of the input
    //_count_field_a is max_assignments x index1 x index2...x index_n
      //_prob_field_a is max_assignments x 1
      //_fields_c =index1 x index2 ... x index_n
      
      if(_count_field_a.size(0)<_num_valid_assignments)
	{
	  std::cout<<"Error: pass in not number of valid assignments or incorrect sized index fields, please check"<<endl;
	  return false;
	}
      
      if(_prob_field_a.dim()!=1)
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
      //unsigned blockSize=_field_c.nTotal();//total number of elements in the counter field
      int* dim_array=new int[_count_field_a.dim()];
      for(unsigned j=1;j<_count_field_a.dim();j++)
	{
	  dim_array[j]=-1;
	}
      for(unsigned i=0;i<_num_valid_assignments;i++)
	{
	  //determine the dim and construct the dim_array for
	  // the sub matrix. basicall, it is the first dimenion i
	  //plus all -1 (take all ) for the rest
	  dim_array[0]=i;
	  
	  Matrix<double> subM=_count_field_a.SubMatrix(_count_field_a.dim(), dim_array);
	  
	  //we got the subMatrix holding the counts
	  //then get the weighted counter by probability of this assignment
	  subM=subM*_prob_field_a(i);
	  //add to the counter
	  _field_c=_field_c+subM;
	}

      //done
      return flag;
  }//end of function

//another version of above one. for the case where nM field is a vector and 
  //    counter field is a scalar
  bool VDJ_cuts_insertion_dinuc_ntbias_model::Update_nM_field(Matrix<double> _count_field_a, Matrix<double> _prob_field_a, 
		       double _fields_c,
		       const unsigned& _num_valid_assignments) const
  {
    bool flag=true;
    //first check for the validity of the input
    //_count_field_a is max_assignments x 1
      //_prob_field_a is max_assignments x 1
      //_fields_c =scalar value
      
      if(_count_field_a.size(0)<_num_valid_assignments)
	{
	  std::cout<<"Error: pass in not number of valid assignments or incorrect sized index fields, please check"<<endl;
	  return false;
	}
      
      if(_prob_field_a.dim()!=1)
	{
	  cout<<"Error: pass in not correct prob_field, please check"<<endl;
	  return false;
	}
      if(_count_field_a.dim()!=1) //_count_field is one more dimension than _field_a
	{
	  cout<<"Error: pass in not correct index_field (dimension not 1), please check"<<endl;
	  return false;
	}

      if(_prob_field_a.size(0)!=_count_field_a.size(0))
	{
	  cout<<"Error: the size of two fields from assing don't match, please check"<<endl;
	  return false;
	}
      
      //now we are good, start doing the updating
      for(unsigned i=0;i<_num_valid_assignments;i++)
	{
	  _fields_c+=_count_field_a(i)*_prob_field_a(i);
	}

      //done
      return flag;
  }

//TODO::::::
//=========??????? will later change the signature, removing _model, but fill in
//with other necessary parameter for building a vdj counter
//_model is only used to build a new Vdj counter

VDJ_cuts_insertion_dinuc_ntbias_counter VDJ_cuts_insertion_dinuc_ntbias_model::SumCounter
(const VDJ_cuts_insertion_dinuc_ntbias_counter& _c1, const VDJ_cuts_insertion_dinuc_ntbias_counter& _c2)
{
  //
  VDJ_cuts_insertion_dinuc_ntbias_counter c_ret(this->model_params);
  
  //copy over all the fields
//==========================>member declaration, from the model
  c_ret.logLikelihood=_c1.logLikelihood+_c2.logLikelihood;//This adds up the log likelihood of generating a sequence given the model
  c_ret.N_processed=_c1.N_processed+_c2.N_processed;

  //% For all parameters of model that begin with a 'P' or an 'M' AND don't have the word 'specific' in them,
  //% make a corresponding variable in counter, with an 'n' prefixed to it.
  c_ret.nPinsVD=_c1.nPinsVD+_c2.nPinsVD;
  c_ret.nPinsDJ=_c1.nPinsDJ+_c2.nPinsDJ;

  c_ret.nPcutV_given_V=_c1.nPcutV_given_V+_c2.nPcutV_given_V;
  c_ret.nPcutJ_given_J=_c1.nPcutJ_given_J+_c2.nPcutJ_given_J;
  c_ret.nPcutDlcutDr_given_D=_c1.nPcutDlcutDr_given_D+_c2.nPcutDlcutDr_given_D;

  c_ret.nPV=_c1.nPV+_c2.nPV;
  c_ret.nPDJ=_c1.nPDJ+_c2.nPDJ; //Joint P(V, D, J gene choices)
  c_ret.nPVallele_given_gene=_c1.nPVallele_given_gene+_c2.nPVallele_given_gene; //Probabilities of alleles given gene for each gene
  c_ret.nPDallele_given_gene=_c1.nPDallele_given_gene+_c2.nPDallele_given_gene;
  c_ret.nPJallele_given_gene=_c1.nPJallele_given_gene+_c2.nPJallele_given_gene;
  
  //_per_
  //% Rate parameter. So make two 'nM' variables in counter, one for numerator and one for denominator.
  //% Names of these two should be separated by '_per_'.
  //% The ratio of these will be estimate of the model parameter.
  c_ret.nMnucleotideVD=_c1.nMnucleotideVD+_c2.nMnucleotideVD; 
  c_ret.nMnucleotideVD_5prime=_c1.nMnucleotideVD_5prime+_c2.nMnucleotideVD_5prime;//nucleotide distr's
  c_ret.nMnucleotideDJ=_c1.nMnucleotideDJ+_c2.nMnucleotideDJ; 
  c_ret.nMnucleotideDJ_3prime=_c1.nMnucleotideDJ_3prime+_c2.nMnucleotideDJ_3prime;

  c_ret.nMerror=_c1.nMerror+_c2.nMerror; 
  c_ret.nMsequenced_nucleotide=_c1.nMsequenced_nucleotide+_c2.nMsequenced_nucleotide ;//error rate, will be 
  //============the above members are from model

  //==========by the counter
  //add any other fields you might want to keep track of, but are not in model,
  //below, begin the name with an 'nP' or 'nM',
  //if it is a distribution or mean value, respectively;
  c_ret.nMmononucleotideVD=_c1.nMmononucleotideVD+_c2.nMmononucleotideVD;// zeros(4,1);
  c_ret.nMmononucleotideDJ=_c1.nMmononucleotideDJ+_c2.nMmononucleotideDJ;// = zeros(4,1);
  c_ret.nMinsertionVD=_c1.nMinsertionVD+_c2.nMinsertionVD;// = zeros(4,1);
  c_ret.nMinsertionDJ=_c1.nMinsertionDJ+_c2.nMinsertionDJ;// = zeros(4,1);

  c_ret.nPVD_left_edge_dinucleotide=_c1.nPVD_left_edge_dinucleotide+_c2.nPVD_left_edge_dinucleotide;// = zeros(4,4);
  c_ret.nPVD_right_edge_dinucleotide=_c1.nPVD_right_edge_dinucleotide+_c2.nPVD_right_edge_dinucleotide;// = zeros(4,4);
  c_ret.nPDJ_left_edge_dinucleotide=_c1.nPDJ_left_edge_dinucleotide+_c2.nPDJ_left_edge_dinucleotide;// = zeros(4,4);
  c_ret.nPDJ_right_edge_dinucleotide=_c1.nPDJ_right_edge_dinucleotide+_c2.nPDJ_right_edge_dinucleotide;// = zeros(4,4);

  c_ret.nMtrinucleotideVD=_c1.nMtrinucleotideVD+_c2.nMtrinucleotideVD;// = zeros(4,4,4);
  c_ret.nMtrinucleotideDJ=_c1.nMtrinucleotideDJ+_c2.nMtrinucleotideDJ;// = zeros(4,4,4);
  
  c_ret.nPpVmax_delV_V=_c1.nPpVmax_delV_V+_c2.nPpVmax_delV_V;// =  zeros(model.max_palindrome + 1, model.max_V_deletions + 1, size(model.PV,1));
  c_ret.nPpJmax_delJ_J=_c1.nPpJmax_delJ_J+_c2.nPpJmax_delJ_J;// = zeros(model.max_palindrome + 1, model.max_J_deletions + 1, size(model.PDJ,2));

  c_ret.nPpDlmax_delDl_D=_c1.nPpDlmax_delDl_D+_c2.nPpDlmax_delDl_D;// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));
  c_ret.nPpDrmax_delDr_D=_c1.nPpDrmax_delDr_D+_c2.nPpDrmax_delDr_D;// = zeros(model.max_palindrome + 1, model.max_D_deletions + 1, size(model.PDJ,1));

  c_ret.nMzeroD=_c1.nMzeroD+_c2.nMzeroD ;
  c_ret.nPVDJ=_c1.nPVDJ+_c2.nPVDJ;// = zeros(size(model.PV,1), size(model.PDJ,1), size(model.PDJ,2));
  c_ret.nPpVdelV=_c1.nPpVdelV+_c2.nPpVdelV;// = zeros( model.max_palindrome + 1, model.max_V_deletions + 1);
  c_ret.nPpDldelDl=_c1.nPpDldelDl+_c2.nPpDldelDl;// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
  c_ret.nPpDrdelDr=_c1.nPpDrdelDr+_c2.nPpDrdelDr;// = zeros( model.max_palindrome + 1, model.max_D_deletions + 1);
  c_ret.nPpJdelJ=_c1.nPpJdelJ+_c2.nPpJdelJ;// = zeros( model.max_palindrome + 1, model.max_J_deletions + 1);


  c_ret.nPVcutV=_c1.nPVcutV+_c2.nPVcutV;// = zeros(size(model.PV,1),model.max_palindrome + model.max_V_deletions + 1);
  c_ret.nPDcutDl=_c1.nPDcutDl+_c2.nPDcutDl;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
  c_ret.nPDcutDr=_c1.nPDcutDr+_c2.nPDcutDr;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_D_deletions + 1);
  c_ret.nPJcutJ=_c1.nPJcutJ+_c2.nPJcutJ;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_J_deletions + 1);

  c_ret.nPDcutV=_c1.nPDcutV+_c2.nPDcutV;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_V_deletions + 1);
  c_ret.nPDcutJ=_c1.nPDcutJ+_c2.nPDcutJ;// = zeros(size(model.PDJ,1),model.max_palindrome + model.max_J_deletions + 1);

  c_ret.nPVcutJ=_c1.nPVcutJ+_c2.nPVcutJ;// = zeros(size(model.PV,1),model.max_palindrome + model.max_J_deletions + 1);
  c_ret.nPVcutDl=_c1.nPVcutDl+_c2.nPVcutDl;// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);
  c_ret.nPVcutDr=_c1.nPVcutDr+_c2.nPVcutDr;// = zeros(size(model.PV,1),model.max_palindrome + model.max_D_deletions + 1);

  c_ret.nPJcutDl=_c1.nPJcutDl+_c2.nPJcutDl;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
  c_ret.nPJcutDr=_c1.nPJcutDr+_c2.nPJcutDr;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_D_deletions + 1);
  c_ret.nPJcutV=_c1.nPJcutV+_c2.nPJcutV;// = zeros(size(model.PDJ,2),model.max_palindrome + model.max_V_deletions + 1);

  c_ret.nPinsVDcutV=_c1.nPinsVDcutV+_c2.nPinsVDcutV;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);
  c_ret.nPinsDJcutV=_c1.nPinsDJcutV+_c2.nPinsDJcutV;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_V_deletions + 1);

  c_ret.nPinsVDcutDl=_c1.nPinsVDcutDl+_c2.nPinsVDcutDl;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  c_ret.nPinsDJcutDl=_c1.nPinsDJcutDl+_c2.nPinsDJcutDl;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

  c_ret.nPinsVDcutDr=_c1.nPinsVDcutDr+_c2.nPinsVDcutDr;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);
  c_ret.nPinsDJcutDr=_c1.nPinsDJcutDr+_c2.nPinsDJcutDr;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_D_deletions + 1);

  c_ret.nPinsVDcutJ=_c1.nPinsVDcutJ+_c2.nPinsVDcutJ;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);
  c_ret.nPinsDJcutJ=_c1.nPinsDJcutJ+_c2.nPinsDJcutJ;// = zeros(model.max_insertions + 1,model.max_palindrome + model.max_J_deletions + 1);

  c_ret.nPV_align_length=_c1.nPV_align_length+_c2.nPV_align_length;// = zeros(max_V_length + 1,1);
  c_ret.nPD_align_length=_c1.nPD_align_length+_c2.nPD_align_length;// = zeros(max_D_length + 1,1);
  c_ret.nPJ_align_length=_c1.nPJ_align_length+_c2.nPJ_align_length;// = zeros(max_J_length + 1,1);

  c_ret.nMVV_err_pos=_c1.nMVV_err_pos+_c2.nMVV_err_pos;// = zeros(size(model.PV,1), max_V_length + model.max_V_deletions);
  c_ret.nMJJ_err_pos=_c1.nMJJ_err_pos+_c2.nMJJ_err_pos;// = zeros(size(model.PDJ,2), max_J_length);

  c_ret.nPJJ_align_length=_c1.nPJJ_align_length+_c2.nPJJ_align_length;// = zeros(size(model.PDJ,2), 1 + max_J_length);
  c_ret.nPVV_align_length=_c1.nPVV_align_length+_c2.nPVV_align_length;// = zeros(size(model.PV,1), 1 + max_V_length);


  c_ret.nPinsDJ_D_align_length=_c1.nPinsDJ_D_align_length+_c2.nPinsDJ_D_align_length;// = zeros(model.max_insertions +1, max_D_length + 1);
  c_ret.nPinsVD_D_align_length=_c1.nPinsVD_D_align_length+_c2.nPinsVD_D_align_length;// = zeros(model.max_insertions +1, max_D_length + 1);

  c_ret.nPinsDJ_J_align_length=_c1.nPinsDJ_J_align_length+_c2.nPinsDJ_J_align_length;// = zeros(model.max_insertions +1, max_J_length + 1);
  c_ret.nPinsVD_J_align_length=_c1.nPinsVD_J_align_length+_c2.nPinsVD_J_align_length;// = zeros(model.max_insertions +1, max_J_length + 1);

  c_ret.nPinsDJ_V_align_length=_c1.nPinsDJ_V_align_length+_c2.nPinsDJ_V_align_length;// = zeros(model.max_insertions +1, max_V_length + 1);
  c_ret.nPinsVD_V_align_length=_c1.nPinsVD_V_align_length+_c2.nPinsVD_V_align_length;// = zeros(model.max_insertions +1, max_V_length + 1);

  c_ret.nPDallele_D_align_length=_c1.nPDallele_D_align_length+_c2.nPDallele_D_align_length;// = zeros(3, max_D_length + 1);


  c_ret.nPdelVinsVD=_c1.nPdelVinsVD+_c2.nPdelVinsVD;// = zeros(model.max_V_deletions +1, model.max_insertions +1);
  c_ret.nPdelVinsDJ=_c1.nPdelVinsDJ+_c2.nPdelVinsDJ;// = zeros(model.max_V_deletions +1, model.max_insertions +1);
  c_ret.nPdelVdelDl=_c1.nPdelVdelDl+_c2.nPdelVdelDl;// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
  c_ret.nPdelVdelDr=_c1.nPdelVdelDr+_c2.nPdelVdelDr;// = zeros(model.max_V_deletions +1, model.max_D_deletions +1);
  c_ret.nPdelVdelJ=_c1.nPdelVdelJ+_c2.nPdelVdelJ;// = zeros(model.max_V_deletions +1, model.max_J_deletions +1);

  c_ret.nPdelJinsVD=_c1.nPdelJinsVD+_c2.nPdelJinsVD;// = zeros(model.max_J_deletions +1, model.max_insertions +1);
  c_ret.nPdelJinsDJ=_c1.nPdelJinsDJ+_c2.nPdelJinsDJ;//zeros(model.max_J_deletions +1, model.max_insertions +1);
  c_ret.nPdelJdelDl=_c1.nPdelJdelDl+_c2.nPdelJdelDl;//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
  c_ret.nPdelJdelDr=_c1.nPdelJdelDr+_c2.nPdelJdelDr;//zeros(model.max_J_deletions +1, model.max_D_deletions +1);
  
  c_ret.nPdelDlinsVD=_c1.nPdelDlinsVD+_c2.nPdelDlinsVD;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  c_ret.nPdelDlinsDJ=_c1.nPdelDlinsDJ+_c2.nPdelDlinsDJ;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  c_ret.nPdelDldelDr=_c1.nPdelDldelDr+_c2.nPdelDldelDr;//zeros(model.max_D_deletions +1, model.max_D_deletions +1);
  
  c_ret.nPdelDrinsVD=_c1.nPdelDrinsVD+_c2.nPdelDrinsVD;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  c_ret.nPdelDrinsDJ=_c1.nPdelDrinsDJ+_c2.nPdelDrinsDJ;//zeros(model.max_D_deletions +1, model.max_insertions +1);
  
  c_ret.nPinsVDinsDJ=_c1.nPinsVDinsDJ+_c2.nPinsVDinsDJ;//zeros(model.max_insertions +1, model.max_insertions +1);
  
  c_ret.nPVdelV=_c1.nPVdelV+_c2.nPVdelV;//zeros(size(model.PV,1), model.max_V_deletions+1 );
  c_ret.nPDdelDl=_c1.nPDdelDl+_c2.nPDdelDl;//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
  c_ret.nPDdelDr=_c1.nPDdelDr+_c2.nPDdelDr;//zeros(size(model.PDJ,1), model.max_D_deletions+1 );
  c_ret.nPJdelJ=_c1.nPJdelJ+_c2.nPJdelJ;//zeros(size(model.PDJ,2), model.max_J_deletions+1 );
  
  c_ret.nPVinsVD=_c1.nPVinsVD+_c2.nPVinsVD;//zeros(size(model.PV,1), model.max_insertions+1 );
  c_ret.nPDinsVD=_c1.nPDinsVD+_c2.nPDinsVD;//zeros(size(model.PDJ,1), model.max_insertions+1 );
  c_ret.nPDinsDJ=_c1.nPDinsDJ+_c2.nPDinsDJ;//zeros(size(model.PDJ,1), model.max_insertions+1 );
  c_ret.nPJinsDJ=_c1.nPJinsDJ+_c2.nPJinsDJ;//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
  c_ret.nPVdelDl=_c1.nPVdelDl+_c2.nPVdelDl;//zeros(size(model.PV,1), model.max_D_deletions+1 );
  c_ret.nPVdelDr=_c1.nPVdelDr+_c2.nPVdelDr;//zeros(size(model.PV,1), model.max_D_deletions+1 );
  c_ret.nPVdelJ=_c1.nPVdelJ+_c2.nPVdelJ;//zeros(size(model.PV,1), model.max_J_deletions+1 );
  
  c_ret.nPJdelV=_c1.nPJdelV+_c2.nPJdelV;//zeros(size(model.PDJ,2), model.max_V_deletions+1 );
  c_ret.nPJdelDl=_c1.nPJdelDl+_c2.nPJdelDl;//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
  c_ret.nPJdelDr=_c1.nPJdelDr+_c2.nPJdelDr;//zeros(size(model.PDJ,2), model.max_D_deletions+1 );
  
  c_ret.nPDdelV=_c1.nPDdelV+_c2.nPDdelV;//zeros(size(model.PDJ,1), model.max_V_deletions+1 );
  c_ret.nPDdelJ=_c1.nPDdelJ+_c2.nPDdelJ;//zeros(size(model.PDJ,1), model.max_J_deletions+1 );
  
  c_ret.nPVinsDJ=_c1.nPVinsDJ+_c2.nPVinsDJ;//zeros(size(model.PV,1), model.max_insertions+1 );
  c_ret.nPJinsVD=_c1.nPJinsVD+_c2.nPJinsVD;//zeros(size(model.PDJ,2), model.max_insertions+1 );
  
  c_ret.nPpVinsVD=_c1.nPpVinsVD+_c2.nPpVinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  c_ret.nPpVinsDJ=_c1.nPpVinsDJ+_c2.nPpVinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  c_ret.nPpVdelDl=_c1.nPpVdelDl+_c2.nPpVdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  c_ret.nPpVdelDr=_c1.nPpVdelDr+_c2.nPpVdelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  c_ret.nPpVdelJ=_c1.nPpVdelJ+_c2.nPpVdelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  c_ret.nPVpV=_c1.nPVpV+_c2.nPVpV;//zeros(size(model.PV,1), model.max_palindrome+1 );
  c_ret.nPJpV=_c1.nPJpV+_c2.nPJpV;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  c_ret.nPDpV=_c1.nPDpV+_c2.nPDpV;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  c_ret.nPpJinsVD=_c1.nPpJinsVD+_c2.nPpJinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  c_ret.nPpJinsDJ=_c1.nPpJinsDJ+_c2.nPpJinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  c_ret.nPpJdelV=_c1.nPpJdelV+_c2.nPpJdelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  c_ret.nPpJdelDl=_c1.nPpJdelDl+_c2.nPpJdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  c_ret.nPpJdelDr=_c1.nPpJdelDr+_c2.nPpJdelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  c_ret.nPVpJ=_c1.nPVpJ+_c2.nPVpJ;//zeros(size(model.PV,1), model.max_palindrome+1 );
  c_ret.nPJpJ=_c1.nPJpJ+_c2.nPJpJ;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  c_ret.nPDpJ=_c1.nPDpJ+_c2.nPDpJ;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  c_ret.nPpDlinsVD=_c1.nPpDlinsVD+_c2.nPpDlinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  c_ret.nPpDlinsDJ=_c1.nPpDlinsDJ+_c2.nPpDlinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  c_ret.nPpDldelV=_c1.nPpDldelV+_c2.nPpDldelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  c_ret.nPpDldelJ=_c1.nPpDldelJ+_c2.nPpDldelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  c_ret.nPpDldelDr=_c1.nPpDldelDr+_c2.nPpDldelDr;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  c_ret.nPVpDl=_c1.nPVpDl+_c2.nPVpDl;//zeros(size(model.PV,1), model.max_palindrome+1 );
  c_ret.nPJpDl=_c1.nPJpDl+_c2.nPJpDl;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  c_ret.nPDpDl=_c1.nPDpDl+_c2.nPDpDl;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  c_ret.nPpDrinsVD=_c1.nPpDrinsVD+_c2.nPpDrinsVD;//zeros(model.max_palindrome +1, model.max_insertions +1);
  c_ret.nPpDrinsDJ=_c1.nPpDrinsDJ+_c2.nPpDrinsDJ;//zeros(model.max_palindrome +1, model.max_insertions +1);
  c_ret.nPpDrdelV=_c1.nPpDrdelV+_c2.nPpDrdelV;//zeros(model.max_palindrome +1, model.max_V_deletions +1);
  c_ret.nPpDrdelJ=_c1.nPpDrdelJ+_c2.nPpDrdelJ;//zeros(model.max_palindrome +1, model.max_J_deletions +1);
  c_ret.nPpDrdelDl=_c1.nPpDrdelDl+_c2.nPpDrdelDl;//zeros(model.max_palindrome +1, model.max_D_deletions +1);
  c_ret.nPVpDr=_c1.nPVpDr+_c2.nPVpDr;//zeros(size(model.PV,1), model.max_palindrome+1 );
  c_ret.nPJpDr=_c1.nPJpDr+_c2.nPJpDr;//zeros(size(model.PDJ,2), model.max_palindrome+1 );
  c_ret.nPDpDr=_c1.nPDpDr+_c2.nPDpDr;//zeros(size(model.PDJ,1), model.max_palindrome+1 );
  
  c_ret.nPpVpDl=_c1.nPpVpDl+_c2.nPpVpDl;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  c_ret.nPpVpDr=_c1.nPpVpDr+_c2.nPpVpDr;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  c_ret.nPpVpJ=_c1.nPpVpJ+_c2.nPpVpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  c_ret.nPpDlpDr=_c1.nPpDlpDr+_c2.nPpDlpDr;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  c_ret.nPpDlpJ=_c1.nPpDlpJ+_c2.nPpDlpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  c_ret.nPpDrpJ=_c1.nPpDrpJ+_c2.nPpDrpJ;//zeros(model.max_palindrome + 1, model.max_palindrome + 1);
  
  //%
  c_ret.nMerror_vs_position=_c1.nMerror_vs_position+_c2.nMerror_vs_position;//zeros(model.read_length,1);
  c_ret.nMcoverage=_c1.nMcoverage+_c2.nMcoverage;//zeros(model.read_length,1);

  
  return c_ret;
}

//=========??????? will later change the signature, removing _model, but fill in
//with other necessary parameter for building a vdj counter
//_model is only used to build a new Vdj counter
VDJ_cuts_insertion_dinuc_ntbias_counter VDJ_cuts_insertion_dinuc_ntbias_model::SumCounter
(const VDJ_cuts_insertion_dinuc_ntbias_counter* _c, const unsigned _size) 
{
  
  VDJ_cuts_insertion_dinuc_ntbias_counter c_ret(this->model_params);

  if(_size==0)
    {
      std::cout<<"WARNING:: passing in the zero sized counter array, return a null one too"<<endl;
      std::cerr<<"WARNING:: passing in the zero sized counter array, return a null one too"<<endl;
      return c_ret;
    }
  //if(_size==1)
  c_ret=_c[0];
  
  //now do the addition
  for(unsigned i=1;i<_size;i++)
    {
      c_ret=SumCounter(c_ret, _c[i]);
    }
  
  return c_ret;
}


void VDJ_cuts_insertion_dinuc_ntbias_model::CalculateAssignmentEntropies()
{
  //=========>to be implemented, for now leave it blank;
}
