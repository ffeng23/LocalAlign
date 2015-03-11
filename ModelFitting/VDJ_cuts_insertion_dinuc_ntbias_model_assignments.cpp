
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

  //

  return true;

}
			   




