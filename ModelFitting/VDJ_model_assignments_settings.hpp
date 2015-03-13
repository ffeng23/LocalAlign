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
  //unsigned max_nerrorsv_1;

  unsigned high_error_region;//don't actually use it in my code,
  //but keep it for now for compatibility purpose

  //for v gene
  unsigned v;//this is the v allele index in the alignment
  unsigned v_a; //v allele number
  unsigned v_g; //v gene index
  Matrix<unsigned> v_err_pos;
  double log_highest_probability_GIVEN_current_V_allele;
  bool v_break_out;
  unsigned n_assignments_v_gene;

  //for J gene
  unsigned j;
  unsigned j_a; //j allele
  unsigned j_g; //J gene index
  double log_max_pcutJ_loop_J;
  double log_highest_probability_GIVEN_current_J_allele;
  bool j_break_out;
  unsigned n_assignment_j_gene;
  
  unsigned niVD_DJ0; //first insertion between VD and DJ
  unsigned niVD_DJ_min;

  //for D gene
  unsigned d;
  unsigned d_a;
  unsigned d_g;
  unsigned l_d_seq;
  unsigned n_assignments_d_gene;
  bool d_break_out=false;
  double log_max_pcutD_loop_d;
  //double log_max_pcutVDJ_loop_d; <=== here for this one, we did not define, since it is the sum of pcutV/D/J together.
  double log_highest_probability_GIVEN_D_allele;

  unsigned p_max_Dl;
  unsigned p_max_Dr;
  unsigned  p_max_J

  double log_probabase;//is the probability of vdj choices with allele
  
};

//NOTE: p is palindrome
//     P is probability

#endif
