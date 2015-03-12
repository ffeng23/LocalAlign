
#include <cmath>
#include "VDJ_cuts_insertion_dinuc_ntbias_model_assignments.hpp"

#include "VDJ_model_assignments_settings.hpp"
#include "../SIGPIG/AlignmentSettings.hpp"

bool VDJ_model_assignments
(
 const VDJ_cuts_insertion_dinuc_ntbias_model& _model,
 const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
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
  assignment_params.max_nerrorsv_1=0;
  

  return true;

}
			   

bool assign_V_alleles
(const VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 /*output, input*/VDJ_model_assignment_settings& _params,
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns )
{
  for(unsigned i=0;i<_V.numOfAligned;i++)
    {
      //check for validity
      if(_params.skips>_params.max_skips)
	{
	  _assign.n_assignments=_params.in;
	  _assign.skips=_params.skips;
	  return true;
	}
    }
  
}





