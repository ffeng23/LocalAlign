#ifndef DO_PROBABILISTIC_MODEL_FITTING_HPP
#define DO_PROBABILISTIC_MODEL_FITTING_HPP

bool do_probabilistic_model_fitting
(
 const SequenceString* _seq,
 const Alignment_Object* _V, const Alignment_D* _D, const Alignment_Object* _J,
 const unsigned numOfAlignments,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 const bool& start_from_flat_prior
 /*const double& _probability_threshold_factor, const bool& _no_error,
   const bool& _ignore_deep_error, const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,*/
 // /*output*/ VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns
 
 );				    


#endif
