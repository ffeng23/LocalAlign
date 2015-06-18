#include "VDJ_cuts_insertion_dinuc_ntbias_model_params.hpp"
#include "../SIGPIG/genomicSegments.hpp"
#include <iostream>
using namespace std;

//VDJ_cuts_insertion_dinuc_ntbias_model_params vdj_mps;

//vdj_mps.max_assignments=5000;
VDJ_cuts_insertion_dinuc_ntbias_model_params::VDJ_cuts_insertion_dinuc_ntbias_model_params
(
const GenomicV* _genV, const unsigned& _numV, 
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ
)
  :
/*initilization list*/
  max_assignments(6000), max_insertions(50),
  max_V_deletions(16), max_D_deletions(16), max_J_deletions(18),
  number_V_genes(0), number_D_genes(0), number_J_genes(0),
  max_V_n_alleles(0), max_D_n_alleles(0), max_J_n_alleles(0),
  max_excess_V_deletions(12),  max_excess_D_deletions(15), max_excess_J_deletions(10),
  max_palindrome(6), 
  min_V_cut(0) /*cut will be set later*/, min_D_cut(0)/*later*/, min_J_cut(0)/*later*/,
  max_V_cut(0), max_D_cut(0),  max_J_cut(0) /*later*/,
  negative_excess_deletions_max (0), /*not sure why we set it up as 0 in here*/
  min_J_align_length(2), min_J_assign_length(1), 
  min_V_length(20)/*originally in matlab code is 15*/,
  high_error_region(15) /*we PROBABLY will NOT use this one*/,
  use_no_D_match_seqs(true),
  read_length(101),//???????????????is this good??? was 101 previous
  maximum_read_length(500)
{
  //empty
  min_V_cut=-1*max_palindrome;
  min_D_cut=-1*max_palindrome;
  min_J_cut=-1*max_palindrome;
  
  max_V_cut=max_V_deletions;
  max_D_cut=max_D_deletions;
  max_J_cut=max_J_deletions;
  //start intialize the various matrix
  //cout<<"---inside  vdj parameter:"<<number_V_genes<<endl;
  number_V_genes=max_gene_index(_genV,_numV)+1;
  //cout<<"---after"<<number_V_genes<<endl;
  //cout<<"---before"<<number_D_genes<<endl;
  number_D_genes=max_gene_index(_genD,_numD)+1;
  //cout<<"---aafter"<<number_D_genes<<endl;
  //cout<<"---before"<<number_J_genes<<endl;
  /*for(unsigned i=0;i<_numJ;i++)
    {
      cout<<_genJ[i].toString()<<endl;
      cout<<_genJ[i].Get_GeneIndex()<<endl;
    }
    cout<<"*****pointer address:"<<_genJ<<endl;*/
  number_J_genes=max_gene_index(_genJ,_numJ)+1;
  //cout<<"---after"<<number_J_genes<<endl;
  max_V_n_alleles = max_n_alleles(_genV, _numV);
  max_D_n_alleles = max_n_alleles(_genD, _numD);//...n_alleles]);
  max_J_n_alleles =max_n_alleles(_genJ, _numJ);
  /*for(unsigned i=0;i<_numJ;i++)
    {
      cout<<_genJ[i].toString()<<endl;
      }*/
}
