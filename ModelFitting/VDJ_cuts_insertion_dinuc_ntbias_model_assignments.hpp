#ifndef VDJ_CUTS_INSERTION_DINUC_NTBIAS_MODEL_ASSIGNMENTS_HPP
#define VDJ_CUTS_INSERTION_DINUC_NTBIAS_MODEL_ASSIGNMENTS_HPP

#include "VDJ_cuts_insertion_dinuc_ntbias_model.hpp"
#include "../SequenceString.hpp"
#include "VDJ_model_assignments_settings.hpp"
#include "VDJ_cuts_insertion_dinuc_ntbias_assigns.hpp"
#include "../SIGPIG/Alignment.hpp"
#include "../SIGPIG/Alignment_V.hpp"
#include "../SIGPIG/Alignment_D.hpp"

//define the function to run the assignment for each alignment from the input
bool VDJ_model_assignments
(VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 const double& _probability_threshold_factor, const bool& _no_error,
 const bool& _ignore_deep_error, const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
 /*output*/ VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns); 
			   

//==============start the function do the assignment
bool assign_VDJ_alleles
(VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 /*const double& _probability_threshold_factor,*/ const bool& _no_error,
 const bool& _ignore_deep_error, const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
 /*output, input*/VDJ_model_assignments_settings& assignment_params,
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns );

bool assign_VJ_deletions
(VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 /*const double& _probability_threshold_factor,*/ const bool& _no_error,
 const bool& _ignore_deep_error, const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
 /*output, input*/VDJ_model_assignments_settings& assignment_params,
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns );
		    
bool assign_VJ_palindrome
(VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 /*const double& _probability_threshold_factor,*/ const bool& _no_error,
 /*const bool& _ignore_deep_error,*/ const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
 /*output, input*/VDJ_model_assignments_settings& assignment_params,
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns);

//start doing assignment with D seg on both sides
bool assign_D
(VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 /*const double& _probability_threshold_factor,*/ const bool& _no_error,
 /*const bool& _ignore_deep_error,*/ const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
 /*output, input*/VDJ_model_assignments_settings& assignment_params,
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns
 );

//start doing assignment with D seg on both sides
bool run_stats_for_assignment
(VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 /*const double& _probability_threshold_factor,*/ const bool& _no_error,
 /*const bool& _ignore_deep_error,*/ const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
 /*output, input*/VDJ_model_assignments_settings& assignment_params,
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns
 );

#endif
