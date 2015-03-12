
#include <cmath>
#include "VDJ_cuts_insertion_dinuc_ntbias_model_assignments.hpp"

#include "VDJ_model_assignments_settings.hpp"
#include "../SIGPIG/AlignmentSettings.hpp"

bool VDJ_model_assignments
(
 const VDJ_cuts_insertion_dinuc_ntbias_model& _model,
 const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _geneV, const unsigned& _numV,
 const GenomicD* _geneD, const unsigned& _numD,
 const GenomicJ* _geneJ, const unsigned& _numJ,
 const double& _probability_threshold_factor, const bool& _no_err,
 const bool& _ignore_deep_error, const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
  /*output*/ VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns
 )
{
  //start preparing for parameters
  VDJ_model_assignments_settings assignment_params;

  assignment_params.log_probability_threshold_factor=log(_probability_threshold_factor);
  assignment_params.log_probability_hopeless_threshold_factor=assignment_params.log_probability_threshold_factor+log(1E-6);
  assignment_params.in=0;
  assignment_params.skips=0;
 
  assignment_params.max_J_depth=AlignmentSettings::max_J_length;
  assignment_params.max_V_depth=_seq.GetLength()+ AlignmentSettings::max_V_length;
  
  assignment_params.J_max_error=0;
  assignment_params.D_max_error=0;

  assignment_params.assume_palindrome_negative_deletions=true;

  assignment_params.log_max_model_p_nt_DJ=matrix_log(max(model.RnucleotideDJ_per_nucleotideDJ_3_prime,2));
  assignment_params.log_max_model_p_nt_VD=matrix_log(max(model.RnucleotideVD_per_nucleotideVD_5_prime,2));
  assignment_params.log_max_model_p_nt=max(assignment_params.log_max_model_p_nt_DJ);
  double temp=max(assignment_params.log_max_model_p_nt_VD);
  if(assignment_params.log_max_model_p_nt<temp)
    {
      assignment_params.log_max_model_p_nt=temp;
    }

  assignment_params.log_RnucleotideDJ_per_nucleotideDJ_3prime=matrix_log(_model.RnucleotideDJ_per_nucleotideDJ_3prime.m2vec());
  assignment_params.log_RnucleotideVD_per_nucleotideVD_5prime=matrix_log(_model.RnucleotideVD_per_nucleotideVD_5prime.m2vec());

  if(_model.Rerror_per_sequenced_nucleotide>0)
    {
      assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3=
	log(_model.Rerror_per_sequenced_nucleotide/3);
    }
  else
    {
      assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3=-1000.0
    }
  
  assignment_params.L_err=_seq.GetLength();
  assignment_params.log_proba_Rerror_normalization=-1.0*L_err*log(1+_model.Rerror_per_sequenced_nucleotide);
  assignment_params.log_max_model_pcutV_given_gene=matrix_log(max(_model.PcutV_given_V,1));
  assignment_params.log_max_model_pcutJ_given_gene=matrix_log(max(_model.PcutJ_given_J,1));
  Matrix<double> temp_m=max(_model.PcutDlcutDr_given_D,2);
  assignment_params.log_max_model_pcutD_given_gene=matrix_log(max(temp_m,1));
  
  assignment_params.log_max_model_pinsVD=log(max(model.PinsVD));
  assignment_params.log_max_model_pinsDJ=log(max(model.PinsDJ));
  assignment_params.log_max_model_pins=assignment_params.log_max_model_pinsVD+assignment_params.log_max_model_pinsDJ;
  
  if(do_smoothing&&!no_error)
    {
      assignment_params.nd_start=-1* _model.negative_excess_deletions_max;
      assignment_params.np_start_from_max=false;
    }
  else
    {
      assignment_params.nd_start=0;
      assignment_params.np_start_from_max=true;
    }
  assignment_params.log_highest_probability=-1000.0;
  assignment_params.best_D_align_length=0;
  //assignment_params.max_nerrorsv_1=0;
  

  return true;

}
			   

bool assign_VDJ_alleles
(const VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const bool ignore_deep_error,
 /*output, input*/VDJ_model_assignment_settings& assignment_params,
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns )
{
  unsigned dim_size[4]={0,0,0,0};
  unsigned deletable_errorss=0;
  unsigned max_nerrorsv_1=0;
  unsigned min_nerrorsv, max_nerrorsv;

  double log_max_pcutV_loop_v;
  for(unsigned v=0;v<_V.numOfAligned;i++)
    {
      assignment_params.v=v;
      //check for validity
      if(assignment_params.skips>assignment_params.max_skips)
	{
	  _assign.n_assignments=assignment_params.in;
	  _assign.skips=assignment_params.skips;
	  return true;
	}
      assignment_params.v_a=_V.alleles_all(v);
      assignment_params.v_g=_genV[v_a].Get_GeneIndex();

      //for high_error_region, we don't use in my code,
      //but keep it for now anyway
      if(ignore_deep_error)
	assignment_params.high_error_region=_model.high_error_region;
      else
	assignment_params.high_error_region=_model.high_error_region;

      if(v==0)//for first round, set some starting point
	{
	  dim_size[0]=_V.n_errors[v];
	  v_err_pos.initialize(1, dim_size, _V.error_positions[v]);
	  deletable_errors=sum_all_bool(v_err_pos > (_V.align_length[v]-_model.max_excess_V_deletions));
	  min_nerrorsv=_V.n_errors[v]-deletable_errors;
	  //here it is different from Matlab code,
	  //since we don't have high_error region anyway.
	  //We simply don't account for it
	  max_nerrorsv_1=_V.n_errors[v];//again, this is different from Matlab code.
	  //since we don't have high error region
	  if(no_error && min_nerrorsv>0)
	    {
	      //we don't allow errors , but do have errors in here, so do next
	      continue;
	    }
	}
      else  //for other cases not v==1
	{
	  //check if this v has too many more eerors than the best V
	  v_error_pos.clear();
	  dim_size[0]=_V.n_errors[v];
	  v_err_pos.initialize(1, dim_size, _V.error_positions[v]);
	  deletable_errors=sum_all_bool(v_err_pos > (_V.align_length[v]-_model.max_excess_V_deletions));
	  min_nerrorsv=_V.n_errors[v]-deletable_errors;
	  //here it is different from Matlab code,
	  //since we don't have high_error region anyway.
	  //We simply don't account for it
	  max_nerrorsv=_V.n_errors[v];//again, this is different from Matlab code.
	  //since we don't have high error region
	  if(no_error && min_nerrorsv>0)
	    {
	      //we don't allow errors , but do have errors in here, so do next
	      continue;
	    }
	  //check for too many more errors compared with best V
	  double log_min_diff_perrv=(max_nerrorsv-max_nerrorsv_1) * assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
	  if(log_min_diff_perrv<assignment.params.log_probability_factor)
	    {
	      assignment_params.skips+=1;
	      continue;//why did not count skips in the above no_error position
	    }
	  
	}//end of else
      log_max_pcutV_loop_V=assignment_params.log_max_model_pcutV_given_gene(v_g);

      assignment_params.log_highest_probability_GIVEN_V_allele=-1000.0;
      assignment_params.v_break_out=false;
      assignment_params.n_assignments_v_gene=0;

      //========loop over J alleles
      for(unsigned j=0;j<_J.numOfAligned;j++)
	{
	  assignment_params.j=j;
	  if(assignment_params.skips>assignment_params.max_skips)
	    {
	      _assigns.n_assignments=assignment_params.in;
	      _assigns.skips=assignment_param.skips;
	      return;
	    }
	  if(_J.align_length(j)<_model.min_J_align_length)
	    continue;
	  assignment_params.j_a=_J.alleles_all[j];
	  assignment_params.j_g=_genJ[assignment_params.j_a].Get_GeneIndex();
	  
	  assignment_params.niVD_DJ0=_J.align_position[j][0] -(_V.align_position[v][0]+ _V.align_length[v]-1)-1;
	  //NOTE:::here the code is different, since the alignment of V not starting from zero. By adding _V.align_position[v][0]

	  assignment_params.log_max_pcutJ_loop_j=assignment_params.log_max_odel_pcutJ_given_gene(j_g);

	  assignment_params.log_highest_probability_GIVEN_current_J_allele=-1000.0;
	  assignment_params.n_assignment_j_gene=0;
	  assignment_params.j_break_out=false;

	  //=======loop over D alleles
	  for(unsigned d_i=0;d_i<_numD;d_i++)
	    {
	      if(assignment_params.skips>assignment_params.max_skips)
		{
		  _assigns.n_assignments=assignment_params.in;
		  _assigns.skips=assignment_param.skips;
		  return;
		}
	      unsigned d=D.allele_order[d_i];
	      assignment_params.d=d;
	      assignment_params.d_a=d;
	      assignment_params.d_g=_genD(assignment_params.d_a).Get_GeneIndex();

	      assignment_params.n_assignment_d_gene=0;
	      assignment_params.d_break_out=false;

	      //base probability of gene choices, multiplied by conditional probability for allele choices, given genes
	      //J have no distinugishable alles for original code, but not for project
	      double probabase=_model.PV((v_g)*_model.PDJ(d_g, j_g)*_model.PVallele
	    }//end of d loop ====
	  
	}//end of j loop =====
      
    }//end of for outer for v loop  ========
  
}





