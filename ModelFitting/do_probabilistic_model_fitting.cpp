#include "do_probabilistic_model_fitting.hpp"

bool do_probabilistic_model_fitting
(
 const SequenceString* _seq,
 const Alignment_Object* _V, const Alignment_D* _D, const Alignment_Object* _J,
 const unsigned numOfAlignments,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 const bool& start_from_flat_prior, const unsigned& _max_iters 
 /*const double& _probability_threshold_factor, const bool& _no_error,
   const bool& _ignore_deep_error, const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,*/
 // /*output*/ VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns
 )
{
  VDJ_cuts_insertion_dinuc_ntbias_model model;
  bool do_smoothing=true;
  bool force_all_alleles=false;
  //start the looping
  for(unsigned iteration_no=0;iteration_no<_max_iters;iteration_no++)
    {
      if((iteration_no==1) && start_from_flat_prior)
	{
	  //place holder, start a new model
	  do_smoothing=true;
	}
      else
	{
	  //read from the disk and continue;
	  do_smoothing=true;
	}

      //set up variables
      if(!do_smoothing)
	{
	  model.max_assignments=6000;
	}
      model.high_error_region=0;
      model.min_V_length=25;
      force_all_alleles=false;
      double probability_threshold_factor;
      if(iteration_no==1&&start_from_flat_prior)
	{
	  probability_threshold_factor=0.01;
	  model.negative_excess_deletions_max=3;//0 in matlab
	}
      else
	{
	  if(iteration_no<5)
	    {
	      probability_threshold_factor=0.005;
	      model.negative_excess_deletions_max=3;//0 in matlab
	    }
	  else
	    {
	      probability_threshold_factor=0.001;
	      model.negative_excess_deletions_max=3;//0 in matlab
	    }
	}
      
      bool no_error=false;
      bool ignore_deep_error=false;
      unsigned READ_LENGTH_CORRECTION=0;
      
      VDJ_cuts_insertion_dinuc_ntbias_counter counter;
      model.InitializeCounter(counter);
      VDJ_cuts_insertion_dinuc_ntbias_assigns assigns;
      model.InitializeAssigns(assigns);

      double* assignment_entropy=new double[numOfAlignments];
      double* assignment_entropy_no_errors=new double[numOfAlignments];
      //now go through the alignments
      for(unsigned k=0;k<numOfAlignments;k++)
	{//for each alignment, do the assign
	  
	  bool assign_flag=VDJ_model_assignments
	    (model, _seq[k], _V[k], _D[k], _J[k], _genV, _numV, _genD, _numD, _genJ, _numJ,
	     probability_threshold_factor, no_error, ignore_deep_error, do_smoothing,
	     force_all_alleles, READ_LENGTH_CORRECTION, _assigns);
	  
	  if(assigns.n_assignments>0)
	    {
	      Matrix<double> p_ass=assigns.proba/assigns.likelihood;
	      assignment_entropy[k]=Entropy(p_ass, exp(1));
	      if(assigns.generation_probability>0)
		{
		  Matrix<double> p_gen=assigns.event_probability/assigns.generation_probability;
		  assignment_entropy_no_errors[k]=entropy(p_gen, exp(1));
		}
	      else
		assignment_entropy_no_errors[k]=nan("");
	      
	      if(assigns.likelihood>0)
		{
		  model.UpdateCounter(assigns, counter);
		}
	      
	    }//end of if >0 loop
	  
	  
	}//end of alignment
      model.GetModelFromCounter(counter);
      model.CalculateAssignmentEntropies();
      
      //save model!!!
    }//end of iteration

  return true;
}				  

