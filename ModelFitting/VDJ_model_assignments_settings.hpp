#ifndef VDJ_MODEL_ASSIGNMENTS_SETTINGS_HPP
#define VDJ_MODEL_ASSIGNMENTS_SETTINGS_HPP

#include "../matrix/Matrix.hpp"

//holding the params to pass around doing the assignment of model
struct VDJ_model_assignments_settings
{
  //constructor
  VDJ_model_assignments_settings();

  //define member
  double log_probability_threshold_factor;
  double log_probability_hopeless_threshold_factor;

  unsigned max_skips; //
  unsigned skips; //running counter of how many skips has happened
  unsigned in; //assignment number
  unsigned deep_error_limit; //50, for now we don't use this, but keep it for the future
  unsigned READ_LENGTH_CORRECTION;//for now, it is zero.

  unsigned max_J_depth;
  unsigned max_V_depth;

  unsigned J_max_error;
  unsigned D_max_error;

  bool assume_palindrome_negative_deletions;

  Matrix<double> log_max_model_p_nt_DJ;
  Matrix<double> log_max_model_p_nt_VD;
  Matrix<double> log_max_model_p_nt;

  Matrix<double> log_RnucleotideDJ_per_nucleotideDJ_3prime;
  Matrix<double> log_RnucleotideVD_per_nucleotideVD_5prime;

  double log_Rerror_per_sequenced_nucleotide_div_by_3;

  unsigned L_err;
  double log_proba_Rerror_normalization;
  Matrix<double> log_max_model_pcutV_given_gene;
  Matrix<double> log_max_model_pcutJ_given_gene;
  Matrix<double> log_max_model_pcutD_given_gene;

  double log_max_model_pinsVD;
  double log_max_model_pinsDJ;
  double log_max_model_pins;

  int nd_start;
  bool np_start_from_max;

  double log_highest_probability;
  unsigned best_D_align_length;
  unsigned max_nerrorsv_1;
};



#endif
