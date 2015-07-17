#include "VDJ_model_assignments_settings.hpp"


VDJ_model_assignments_settings::VDJ_model_assignments_settings()
  :
  log_probability_threshold_factor(0),
  log_probability_hopeless_threshold_factor(0),
  max_skips(0), skips(0),deep_error_limit(0),
  J_max_error(0), D_max_error(0), assume_palindrome_negative_deletions(true),
  
  log_max_model_p_nt_DJ(),
  log_max_model_p_nt_VD(),
  log_max_model_p_nt(0),
  log_RnucleotideDJ_per_nucleotideDJ_3prime(),
  log_RnucleotideVD_per_nucleotideVD_5prime(),
  log_max_model_pcutV_given_gene(),
  log_max_model_pcutJ_given_gene(),
  log_max_model_pcutD_given_gene(),
  np_start_from_max(true),
  v_err_pos(), 
  j_ex_errs_i(),
  d_ex_errs_left_i(),
  d_errs_i(),
  d_ex_errs_right_i()
  //d_errs_i()
{
  //initialize everything?


}
