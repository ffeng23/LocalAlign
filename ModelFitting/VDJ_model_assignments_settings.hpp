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
  //unsigned READ_LENGTH_CORRECTION;//for now, it is zero.

  unsigned max_J_depth;
  unsigned max_V_depth;

  unsigned J_max_error;
  unsigned D_max_error;

  bool assume_palindrome_negative_deletions;

  Matrix<double> log_max_model_p_nt_DJ;
  Matrix<double> log_max_model_p_nt_VD;
  double log_max_model_p_nt;

  Matrix<double> log_RnucleotideDJ_per_nucleotideDJ_3prime;
  Matrix<double> log_RnucleotideVD_per_nucleotideVD_5prime;

  double log_Rerror_per_sequenced_nucleotide_divided_by_3;

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
  double log_max_pcutV_loop_v;

  unsigned V_align_length;
  unsigned J_align_length;
  unsigned D_align_length;
  //for J gene
  unsigned j;
  unsigned j_a; //j allele
  unsigned j_g; //J gene index
  double log_max_pcutJ_loop_j;
  double log_highest_probability_GIVEN_current_J_allele;
  bool j_break_out;
  unsigned n_assignments_j_gene;
  
  unsigned niVD_DJ0; //first insertion between VD and DJ
  //unsigned niVD_DJ_min;
  unsigned nerrorsv;
  unsigned v_ex_errs;
  Matrix<bool> v_ex_errs_i;

  //for D gene
  unsigned d;
  unsigned d_a;
  unsigned d_g;
  unsigned l_d_seq;
  unsigned n_assignments_d_gene;
  bool d_break_out;
  double log_max_pcutD_loop_d;
  //double log_max_pcutVDJ_loop_d; <=== here for this one, we did not define, since it is the sum of pcutV/D/J together.
  double log_highest_probability_GIVEN_current_D_allele;

  unsigned p_max_Dl;
  unsigned p_max_Dr;
  unsigned p_max_J;

  double log_probabase;//is the probability of vdj choices with allele

  //v deletions
  int ndV; int ndV1;
  double log_perrv;
  double log_highest_probability_GIVEN_current_V_deletions;

  //palidrome
  unsigned npV_max;
  unsigned npV_potential_max;
  unsigned npJ_max;
  unsigned npJ_potential_max;

  //cut variable
  double log_PcutV;
  double log_PcutJ;
  //==> double log_max_pcutVDJ_loop_pJ; we don't remember this since it is
  // the sum of  assignment_params.log_PcutV+
  // assignment_params.log_PcutJ+assignment_params.log_max_pcutD_loop_d;
  //and log_max_pcutD_loop_d is the expected model.log_max_model_pcutD_given_gene(d);

  //j deletions
  int ndJ1;
  int ndJ;
  
  Matrix<bool> j_ex_errs_i;
  unsigned j_ex_errs;
  unsigned nerrorsj;
  double log_perrj;
  double log_highest_probability_GIVEN_current_J_deletions;

  unsigned npV;
  unsigned npDl;
  unsigned npDr;
  unsigned npJ;
  
  unsigned nerrorsd;
  bool zeroD;
  unsigned genic_length;

  //===>no nerrors define here, since it is simply the sum of all v d j errors
  //same thing as log_perr;
  int V_end, J_start;  

  //for D alignment 
  unsigned n_D_aligns;
  unsigned start_n_D_aligns;
  unsigned end_n_D_aligns;

  int ndDl, ndDl1;
  int ndDr, ndDr1;

  unsigned npDl_potential_max;
  unsigned npDr_potential_max;

  Matrix<bool> d_ex_errs_left_i;
  //not for d_err_excess_pos_left, we simply don't record in here, but instead do it on the fly, since it is basically pointing to the positions.
  unsigned d_ex_errs_left;

  Maxtrix<bool> d_errs_i;
  unsigned d_errs;

  Matrix<bool> d_ex_errs_right_i;
  unsigned d_ex_errs_right;

  double log_highest_probability_GIVEN_current_Dl_deletions;
  double log_highest_probability_GIVEN_current_Dr_deletions;
};

//NOTE: p is palindrome
//     P is probability

#endif

