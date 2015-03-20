#include "do_probabilistic_model_fitting.hpp"
#include "VDJ_cuts_insertion_dinuc_ntbias_model_params.hpp"
#include "Entropy.hpp"
#include "VDJ_cuts_insertion_dinuc_ntbias_model_assignments.hpp"
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
  cout<<"start the model fitting......"<<endl;
  //define and set the parameter passing to the model
  //
  /*for(unsigned i=0;i<_numJ;i++)
    {
      cout<<_genJ[i].toString()<<endl;
      }*/
  VDJ_cuts_insertion_dinuc_ntbias_model_params vdj_mps
    ( _genV, _numV, _genD, _numD, _genJ, _numJ);
  //set the params special value instead of the default one
  //
  //empty for now?
  //==========
  /*cout<<"=======intialize model......."<<endl;
  for(unsigned i=0;i<_numJ;i++)
    {
      cout<<_genJ[i].toString()<<endl;
      }*/
  VDJ_cuts_insertion_dinuc_ntbias_model model
    (_genV, _numV, _genD, _numD, _genJ, _numJ, vdj_mps
     );
  bool do_smoothing=true;
  bool force_all_alleles=false;
  //start the looping
  cout<<"ready to do the iterations...."<<endl;
  for(unsigned iteration_no=0;iteration_no<_max_iters;iteration_no++)
    {
      cout<<"iteration loop==>"<<iteration_no<<endl;
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
      //cout<<"\t modeling loop 1"<<endl;
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
      //cout<<"\t modeling loop 1"<<endl;
      bool no_error=false;
      bool ignore_deep_error=false;
      unsigned READ_LENGTH_CORRECTION=0;
      
      VDJ_cuts_insertion_dinuc_ntbias_counter counter(vdj_mps);
      model.InitializeCounter(counter);
      //cout<<"\t modeling loop 2"<<endl;
      VDJ_cuts_insertion_dinuc_ntbias_assigns assigns(model,counter);
      model.InitializeAssign(assigns);
      //cout<<"\t modeling loop 3"<<endl;
      model.assignment_entropy=new double[numOfAlignments];
      model.assignment_entropy_no_errors=new double[numOfAlignments];
      //now go through the alignments
      for(unsigned k=0;k<numOfAlignments;k++)
	{//for each alignment, do the assign
	  //cout<<"\t\t alignment loop:"<<k<<endl;
	  bool assign_flag=VDJ_model_assignments
	    (model, _seq[k], _V[k], _D[k], _J[k], _genV, _numV, _genD, _numD, _genJ, _numJ,
	     probability_threshold_factor, no_error, ignore_deep_error, do_smoothing,
	     force_all_alleles, READ_LENGTH_CORRECTION, assigns);

	  if(assign_flag)
	    {
	      cout<<"Assignment Done"<<endl;
	    }
	  else
	    {
	      cout<<"assignement failed with a return value of false, please check, but we assume nothing bad happend so far"<<endl;
	    }
	  cout<<"n of assigned:"<<assigns.n_assignments<<endl;
	  cout<<" assignment likelihood:"<<assigns.likelihood<<endl;
	  if(assigns.n_assignments>0)
	    {
	      //cout<<"inside updateingggggg block............"<<endl;
	      Matrix<double> p_ass=assigns.proba/assigns.likelihood;
	      model.assignment_entropy[k]=Entropy(p_ass, exp(1));
	      if(assigns.generation_probability>0)
		{
		  Matrix<double> p_gen=assigns.event_probability/assigns.generation_probability;
		  model.assignment_entropy_no_errors[k]=Entropy(p_gen, exp(1));
		}
	      else
		model.assignment_entropy_no_errors[k]=nan("");
	      
	      if(assigns.likelihood>0)
		{
		  //cout<<"---------->update counter"<<endl;
		  model.UpdateCounter(assigns, counter);
		}
	      
	    }//end of if n_assignment >0 loop
	  
	  
	}//end of alignment
      //cout<<"PV:"<<model.PV.toString()<<endl;
      //done with one round of all alignment/iteration, need to update the model
      cout<<"Done with one iteration, up date model......"<<endl;
      cout<<"counter.nPV:"<<counter.nPV.toString()<<endl;
      cout<<"counter.nPinsVD:"<<counter.nPinsVD.toString()<<endl;
      model.GetModelFromCounter(counter);
      //cout<<"PV:"<<model.PV.toString()<<endl;
      model.CalculateAssignmentEntropies();
      cout<<"finish one iteration"<<endl;
      //cout<<"PV:"<<model.PV.toString()<<endl;
      //save model!!!
    }//end of iteration
  cout<<"writing the model data to disk, \"vdj_model_output.txt\"......."<<endl;
  model.output(numOfAlignments);

  return true;
}				  

