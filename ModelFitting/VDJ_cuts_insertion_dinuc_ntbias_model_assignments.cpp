#include <exception>
#include <stdexcept>
#include <cmath>
#include "VDJ_cuts_insertion_dinuc_ntbias_model_assignments.hpp"

#include "VDJ_model_assignments_settings.hpp"
#include "../SIGPIG/AlignmentSettings.hpp"

//always return true;
bool VDJ_model_assignments
(
 VDJ_cuts_insertion_dinuc_ntbias_model& _model,
 const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 const double& _probability_threshold_factor, const bool& _no_error,
 const bool& _ignore_deep_error, const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
  /*output*/ VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns
 )
{
  //cout<<"&&&inside assignments:"<<endl;
  //start preparing for parameters
  VDJ_model_assignments_settings assignment_params;

  assignment_params.log_probability_threshold_factor=log(_probability_threshold_factor);
  assignment_params.log_probability_hopeless_threshold_factor=assignment_params.log_probability_threshold_factor+log(1E-6);

  assignment_params.max_skips=300000;
  assignment_params.in=0;
  assignment_params.skips=0;
  assignment_params.deep_error_limit=50;//will not use it.???
  //cout<<"&&&inside assignments:2"<<endl;
    
  //here for setting up the max_depth, we using maximum possible length of each 
  //gene segments. this is different from Matlab code.
  assignment_params.max_J_depth=AlignmentSettings::max_J_length;
  assignment_params.max_V_depth=AlignmentSettings::max_V_length;
  
  assignment_params.J_max_error=6;
  assignment_params.D_max_error=3;

  assignment_params.assume_palindrome_negative_deletions=true;
//cout<<"&&&inside assignments:3"<<endl;
  assignment_params.log_max_model_p_nt_DJ=matrix_log(max(_model.RnucleotideDJ_per_nucleotideDJ_3prime,1));
  //cout<<"&&&inside assignments:3.4"<<endl;
  assignment_params.log_max_model_p_nt_VD=matrix_log(max(_model.RnucleotideVD_per_nucleotideVD_5prime,1));
  assignment_params.log_max_model_p_nt=max(assignment_params.log_max_model_p_nt_DJ);
  double temp=max(assignment_params.log_max_model_p_nt_VD);
  //cout<<"&&&inside assignments:4"<<endl;
  if(assignment_params.log_max_model_p_nt<temp)
    {
      assignment_params.log_max_model_p_nt=temp;
    }
//cout<<"&&&inside assignments:5"<<endl;
  assignment_params.log_RnucleotideDJ_per_nucleotideDJ_3prime=matrix_log(_model.RnucleotideDJ_per_nucleotideDJ_3prime.m2vec());
  assignment_params.log_RnucleotideVD_per_nucleotideVD_5prime=matrix_log(_model.RnucleotideVD_per_nucleotideVD_5prime.m2vec());
//cout<<"&&&inside assignments:6"<<endl;
  if(_model.Rerror_per_sequenced_nucleotide>0)
    {
      assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3=
	log(_model.Rerror_per_sequenced_nucleotide/3);
    }
  else
    {
      assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3=-1000.0;
    }
  //cout<<"&&&inside assignments:7"<<endl;  
  assignment_params.L_err=_seq.GetLength()-_model.high_error_region;
  assignment_params.log_proba_Rerror_normalization=-1.0*assignment_params.L_err*log(1+_model.Rerror_per_sequenced_nucleotide);
  assignment_params.log_max_model_pcutV_given_gene=matrix_log(max(_model.PcutV_given_V,1));
  assignment_params.log_max_model_pcutJ_given_gene=matrix_log(max(_model.PcutJ_given_J,1));
  //cout<<"&&&inside assignments:8"<<endl;
  Matrix<double> temp_m=max(_model.PcutDlcutDr_given_D,2);
  assignment_params.log_max_model_pcutD_given_gene=matrix_log(max(temp_m,1));
  
  assignment_params.log_max_model_pinsVD=log(max(_model.PinsVD));
  assignment_params.log_max_model_pinsDJ=log(max(_model.PinsDJ));
  assignment_params.log_max_model_pins=assignment_params.log_max_model_pinsVD+assignment_params.log_max_model_pinsDJ;
  //cout<<"&&&inside assignments:9"<<endl;
  if(_do_smoothing&&!_no_error)
    {
      assignment_params.nd_start=-1* _model.negative_excess_deletions_max;//usually is -3, but weird enough, in matlab code it is set to be zero?? anyway, we reset it to be 3 at 
      //do_probalistic_modeling functions before calling this here.
      assignment_params.np_start_from_max=false;
    }
  else
    {
      assignment_params.nd_start=0;
      assignment_params.np_start_from_max=true;
    }
  //cout<<"&&&inside assignments:10"<<endl;
  assignment_params.log_highest_probability=-1000.0;
  assignment_params.best_D_align_length=0;
  //assignment_params.max_nerrorsv_1=0;
  //====VDJ alleles
  cout<<"&&&inside assignments:VDJ allele"<<endl;

  //calling to do the vdj allele assignments. it is possibly to return false.
  //but we don't care here, since nothing to do in the false cases.we simply
  //ignore the return values
  assign_VDJ_alleles(_model, _seq, _V, _D, _J, _genV, _numV, _genD, _numD,
		     _genJ, _numJ, _no_error, _ignore_deep_error, _do_smoothing,
		     _force_all_alleles, _READ_LENGTH_CORRECTION,
		     assignment_params, _assigns
		     );
    
  return true;

}
			   
//it is possibly to return false.
//for cases where there are too many skips or assignments
//
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
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns )
{
  //cout<<">>>>>inside VDJ:1"<<endl;
  unsigned dim_size[4]={0,0,0,0};
  unsigned deletable_errors=0;
  unsigned max_nerrorsv_1=0;
  unsigned min_nerrorsv, max_nerrorsv;

  //double log_max_pcutV_loop_v;
  cout<<"seq :"<<_seq.toString()<<endl;
  cout<<"_V.numOfAligned:"<<_V.numOfAligned<<endl;
  for(unsigned v=0;v<_V.numOfAligned;v++)
    {
      cout<<"---->>>>loop v :"<<v<<endl;
      assignment_params.v=v;
      //check for validity
      if(assignment_params.skips>assignment_params.max_skips)
	{
	  _assigns.n_assignments=assignment_params.in;
	  _assigns.skips=assignment_params.skips;
	  //assignment_params.v_break_out=true; //we don't set it here, since it is not useful here
	  return false;
	}
      //cout<<"\tinside VDJ:1"<<endl;
      cout<<_V.toString()<<endl;
      
      assignment_params.v_a=_V.alleles_all[v];
      assignment_params.v_g=_genV[assignment_params.v_a].Get_GeneIndex();
      
      cout<<"v seg:"<<_genV[assignment_params.v_a].toString()<<endl;
      //for high_error_region, we don't use in my code,
      //but keep it for now anyway
      if(_ignore_deep_error)
	assignment_params.high_error_region=_model.high_error_region;
      else
	assignment_params.high_error_region=_model.high_error_region;
      //count<<"\tinside VDJ:2"<<endl;

      //this following check is moved to here from below. this is used to check for the correct ndV_max.
      //we moved this second part to here, since now we can stop earlier if this is not reasonable
      if(/*(_model.max_V_deletions<_V.min_deletions[assignment_params.v])||*/(_V.align_length[assignment_params.v]<_model.min_V_length))
	{
	  //something wrong, skip this one
	  //we did not count this a skip
	  cout<<"******BREAK OUT:  align length too short!!"<<endl;
	  cout<<" v align length:"<<_V.align_length[assignment_params.v]<<"mode min v length:"<<_model.min_V_length<<endl;
	  return true;//like a continue;
	}
      
      if(v==0)//for first round, set some starting point
	{
	  cout<<"in v=0"<<endl;
	  dim_size[0]=_V.n_errors[v];

	  //here, we try to set up the upper or lower limit of possible nerrors
	  //the reason that we do this is that we are not doing deletions yet at this point. so we 
	  //have to figure out the limits.
	  if(_V.n_errors[v]==0) //no errors
	    {
	      deletable_errors=0;
	    }
	  else
	    {
	      assignment_params.v_err_pos.initialize(1, dim_size, _V.error_positions[v]);
	      deletable_errors=sum_all_bool(assignment_params.v_err_pos > (_V.align_length[v]-_model.max_excess_V_deletions));
	    }


	  //cout<<"could be her"<<endl;
	  min_nerrorsv=_V.n_errors[v]-deletable_errors;
	  //here it is different from Matlab code,
	  //since we don't have high_error region anyway.
	  //We simply don't account for it
	  max_nerrorsv_1=_V.n_errors[v];//again, this is different from Matlab code.
	  //since we don't have high error region
	  if(_no_error && min_nerrorsv>0)
	    {
	      //we don't allow errors , but do have errors in here, so do next
	      cout<<"Quit here for the v==0 case;no_error"<<_no_error<<", min_nerrorsv:"<<min_nerrorsv<<endl; 
	      continue;
	    }
	  //cout<<"v=0,1"<<endl;
	}
      else  //for other cases not v==1
	{
	  //check if this v has too many more eerors than the best V
	  assignment_params.v_err_pos.clear();
	  dim_size[0]=_V.n_errors[v];
	  if(dim_size[0]==0)
	    {
	      deletable_errors=0;
	    }
	  else
	    {
	      assignment_params.v_err_pos.initialize(1, dim_size, _V.error_positions[v]);
	      deletable_errors=sum_all_bool(assignment_params.v_err_pos > (_V.align_length[v]-_model.max_excess_V_deletions));
	    }
	  min_nerrorsv=_V.n_errors[v]-deletable_errors;
	  //here it is different from Matlab code,
	  //since we don't have high_error region anyway.
	  //We simply don't account for it
	  max_nerrorsv=_V.n_errors[v];//again, this is different from Matlab code.
	  //since we don't have high error region
	  if(_no_error && min_nerrorsv>0)
	    {
	      cout<<"Quit here for the v==0 case;no_error"<<_no_error<<", min_nerrorsv:"<<min_nerrorsv<<endl; 
	      //we don't allow errors , but do have errors in here, so do next
	      continue;
	    }
	  //check for too many more errors compared with best V
	  double log_min_diff_perrv=(max_nerrorsv-max_nerrorsv_1) *assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
	  if(log_min_diff_perrv<assignment_params.log_probability_threshold_factor)
	    {
	      assignment_params.skips+=1;
	      cout<<"too many errors compared with the best V"<<endl;
	      //cout<<"max_nerrorsv:"<<max_nerrorsv<<";max_nerrorsv_1:"<<max_nerrorsv_1<<endl;
	      continue;//why did not count skips in the above no_error position
	    }
	  
	}//end of else

      //cout<<"\tinside VDJ:3"<<endl;
      assignment_params.log_max_pcutV_loop_v=assignment_params.log_max_model_pcutV_given_gene(assignment_params.v_g);

      assignment_params.log_highest_probability_GIVEN_current_V_allele=-1000.0;
      assignment_params.v_break_out=false;
      assignment_params.n_assignments_v_gene=0;
      //cout<<"\tinside VDJ:4"<<endl;
      //========loop over J alleles
      cout<<"=============================j numOf aligned:"<<_J.numOfAligned<<endl;
      for(unsigned j=0;j<_J.numOfAligned;j++)
	{
	  cout<<"\tloop j:"<<j<<endl;
	  assignment_params.j=j;
	  if(assignment_params.skips>assignment_params.max_skips)
	    {
	      _assigns.n_assignments=assignment_params.in;
	      _assigns.skips=assignment_params.skips;
	      assignment_params.v_break_out=true;
	      return false; //=====is true good here?no, has to be false
	    }
	  //cout<<"\tinside VDJ:j-1"<<endl;
	  if(_J.align_length[j]<_model.min_J_align_length)
	    continue;
	  assignment_params.j_a=_J.alleles_all[j];
	  assignment_params.j_g=_genJ[assignment_params.j_a].Get_GeneIndex();
	  cout<<"j alignment:"<<_J.toString()<<endl;
	  cout<<"j seg:"<<_genJ[assignment_params.j_a].toString()<<endl;
	  assignment_params.niVD_DJ0=_J.align_position[j][0] -(_V.align_position[v][0]+ _V.align_length[v]-1)-1;
	  //cout<<"\tinside VDJ:j-2"<<endl;
	  //NOTE:::here the code is different, since the alignment of V not starting from zero. By adding _V.align_position[v][0]

	  assignment_params.log_max_pcutJ_loop_j=assignment_params.log_max_model_pcutJ_given_gene(assignment_params.j_g);
	  //cout<<"\tinside VDJ:j-3"<<endl;
	  assignment_params.log_highest_probability_GIVEN_current_J_allele=-1000.0;
	  assignment_params.n_assignments_j_gene=0;
	  assignment_params.j_break_out=false;
	  //cout<<"\tinside VDJ:j-4"<<endl;
	  //=======loop over D alleles
	  for(unsigned d_i=0;d_i<_numD;d_i++)
	    {
	      cout<<"\t\tloop D (numD):"<<d_i<<"("<<_numD<<")"<<endl;
	      if(assignment_params.skips>assignment_params.max_skips)
		{
		  _assigns.n_assignments=assignment_params.in;
		  _assigns.skips=assignment_params.skips;
		  assignment_params.v_break_out=true;//in this case, just all out
		  return false; //====is true good here?
		}
	      unsigned d=_D.allele_order[d_i];
	      assignment_params.d=d;
	      assignment_params.d_a=d;
	      assignment_params.d_g=_genD[assignment_params.d_a].Get_GeneIndex();

	      cout<<"_D seg:"<<_genD[d].toString()<<endl;
	      //cout<<"align n_D_allele:"<<_D.n_D_alleles<<endl;
	      //cout<<"align numOfAligned:"<<_D.numOfAligned[d_i]<<endl;

	      if(d_i==0)
		cout<<"D alignment:"<<_D.toString()<<endl;

	      assignment_params.n_assignments_d_gene=0;
	      assignment_params.d_break_out=false;
	      //cout<<"\t\tinside VDJ:D-1"<<endl;
	      //base probability of gene choices, multiplied by conditional probability for allele choices, given genes
	      //J have no distinugishable alles for original code, but not for project
	      //cout<<"assignment_params.v_g:"<<assignment_params.v_g<<endl;
	      //cout<<"model.PV"<<_model.PV.toString()<<endl;;
	      
	      double probabase=_model.PV(assignment_params.v_g);
	      //cout<<"\tinside VDJ:D_1.1"<<endl;
	      probabase=probabase*_model.PDJ(assignment_params.d_g, assignment_params.j_g);
	      //cout<<"\tinside VDJ:D_1.2"<<endl;
	      probabase*=_model.PVallele_given_gene(assignment_params.v_g, _genV[assignment_params.v_a].Get_Allele() );
	      //cout<<"\tinside VDJ:D_1.3"<<endl;
	      probabase*=_model.PDallele_given_gene(assignment_params.d_g, _genD[assignment_params.d_a].Get_Allele());
	      //cout<<"\tinside VDJ:D_1.4"<<endl;
	      probabase*=_model.PJallele_given_gene( assignment_params.j_g, _genJ[assignment_params.j_a].Get_Allele());

	      //cout<<"probabase:---->"<<probabase<<endl;
	      //cout<<"\tinside VDJ:D_2"<<endl;
	      //code change*****, different from matlab, since we add PJallele_given_gene. 
	      //since we can distinguish J allele now.
	      //unsigned dim_pos[]={AlignmentSettings::D_maximum_deletion+1,AlignmentSettings::D_maximum_deletion+1,0};
	      assignment_params.p_max_Dl=max_mf2(_D.p_region_max_length_left[d],_D.numOfAligned[d],AlignmentSettings::D_maximum_deletion+1);
	      assignment_params.p_max_Dr=max_mf2(_D.p_region_max_length_left[d],_D.numOfAligned[d],AlignmentSettings::D_maximum_deletion+1);
	      //dim_pos[0]=j;
	      assignment_params.p_max_J=max_mf(_J.p_region_max_length[j],AlignmentSettings::J_maximum_deletion+1);
	      //count<<"\tinside VDJ:D_3"<<endl;
	      //upper boun on insertion nt bias factor
	      int niVD_DJ_min=0;
	      //	      dim_pos[0]=v;
	      int temp_int=assignment_params.niVD_DJ0- _D.align_length[d][0]-max_mf(_V.p_region_max_length[v], AlignmentSettings::V_maximum_deletion+1)-assignment_params.p_max_J - assignment_params.p_max_Dl- assignment_params.p_max_Dr;
	      if(temp_int>0)
		niVD_DJ_min=temp_int;
	      double log_p_max_nt_VD_DJ_d_allele=niVD_DJ_min*assignment_params.log_max_model_p_nt;
	      //cout<<"\t\tinside VDJ:D_4"<<endl;
	      //cout<<"log_probabase:"<<log(probabase)<<endl;
	      //cout<<"log_highest_probability:"<<assignment_params.log_highest_probability<<endl;
	      assignment_params.log_max_pcutD_loop_d=assignment_params.log_max_model_pcutD_given_gene(assignment_params.d_g);
	      //assignment_params.log_max_pcutVDJ_loop_d=assignment_params.log

	      if(probabase==0||
		 (log(probabase)+
		  assignment_params.log_max_pcutV_loop_v+
		  assignment_params.log_max_pcutD_loop_d+
		  assignment_params.log_max_pcutJ_loop_j+
		  assignment_params.log_max_model_pins+log_p_max_nt_VD_DJ_d_allele)<
		 (assignment_params.log_highest_probability
		  +assignment_params.log_probability_threshold_factor))
		{
		  
		  assignment_params.skips+=1;
		  continue;
		}
	      assignment_params.log_probabase=log(probabase);
		
	      assignment_params.log_highest_probability_GIVEN_current_D_allele=-1000.0;
	      //cout<<"\t\tinside VDJ:D_5"<<endl;
	      //=======now loop over V deletions
	      if(!assign_VJ_deletions(_model, _seq, _V, _D, _J,
				      _genV, _numV, _genD, _numD, _genJ, _numJ,
				      _no_error,
				      _ignore_deep_error, _do_smoothing,
				      _force_all_alleles, _READ_LENGTH_CORRECTION,
				      assignment_params,_assigns
				      )
		 )
		{
		  /*if(assignment_params.v_break_out)//||assignment_param.j_break_out||assignment_param.d_break.out)
		    {
		      return false;
		    }
		  else
		    {
		      break;
		      }*/
		  return false;
		  
		}
	      if(assignment_params.v_break_out||assignment_params.j_break_out/*||assignment_params.d_break_out*/)
		{
		  break;
		}
	      //cout<<">>>>>>>>end of d loop"<<endl;
	    }//end of d loop ====
	  if(assignment_params.v_break_out/*||assignment_params.j_break_out*/)
	    {
	      break;
	    }
	}//end of j loop =====
      /*if(assignment_params.v_break_out)
	{
	  break;
	  }*/
    }//end of for outer for v loop  ========
  cout<<"n_assign:"<<assignment_params.in<<endl;
  _assigns.n_assignments=assignment_params.in;
  _assigns.skips=assignment_params.skips;
  return true;//return true or false doesn't matter
}//end of function of assign VDJ alleles

//-----------------------------------------------------
//************************************************
//return false if something not right, like too many skips and so we
//want the outer caller finish too. just stop doing this alignment.
/*    only 3 cases, a false value will be returned. 1) skips and assigments is too large
 *    at V deletions, 2)same thing for J deletions; 3)false from inner function call.
 *    other cases are all giving out true; for the third case above, it is basically the 
 *    covers the similar cases, all for skips and assignments for the inner function calls
 *
 */
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
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns )

{
  //% Loop over V deletions.
  //% ndV1 is the deviation from the no. of deletions implied by the alignment (min_deletions).
  //% Lower bound is nd_start (usually -3).
  //% Upper bound is set by maximum deletions allowed and by minimum match length required.
  cout<<"<<<<inside DJ deletions:"<<endl;
  unsigned ndV_max=_model.max_excess_V_deletions;

  //this following checking was not in the original code. but it intends to check to make sure the condintion in setting up
  //the ndV_max is not zero. If it is zero or negative, then in the original code, the loop below will not iterate.
  //in here, we add this checking because we want to avoid (unsign) vs. (signed) integer code. 
  //also the second part of the checking is moved to the above code so that we can quit early.
  //it is making more sense!!
  if((_model.max_V_deletions<_V.min_deletions[assignment_params.v])/*||(_V.align_length[assignment_params.v]<_model.min_V_length)*/)
    {
      //something wrong, skip this one
      //we did not count this a skip
      cout<<"******BREAK OUT: due to minimum deletion too big or align length too short!!"<<endl;
      cout<<"model max v deletions:"<<_model.max_V_deletions<<";_V min deletion:"<<_V.min_deletions[assignment_params.v]
	  <<endl;
      return true;//like a continue;
    }
  //get a reasonable and accurate ndV_max value, was the min() part of the matlab code.
  if((signed)ndV_max>(signed)(_V.align_length[assignment_params.v]-_model.min_V_length)&&(_V.align_length[assignment_params.v]>=_model.min_V_length))
    {
      ndV_max=_V.align_length[assignment_params.v]-_model.min_V_length;
    }
  if((signed)ndV_max>(signed)(_model.max_V_deletions-_V.min_deletions[assignment_params.v])&&(_model.max_V_deletions>=_V.min_deletions[assignment_params.v]))
    {
      ndV_max=_model.max_V_deletions-_V.min_deletions[assignment_params.v];
    }

  
  //loop over V deletions
  for(assignment_params.ndV1=assignment_params.nd_start;assignment_params.ndV1<(signed)ndV_max;assignment_params.ndV1++)
    {
      cout<<"\t\t\t\tloop v_deletion:"<<assignment_params.ndV1<<endl;
      if (assignment_params.skips > assignment_params.max_skips)
	{
	  _assigns.n_assignments = assignment_params.in;
	  _assigns.skips =  assignment_params.skips;
	  assignment_params.v_break_out=true;//===in this case, we want to jump all out.
	  return false;
	}//       end
                
      //% Actual number of deletions
      assignment_params.ndV = _V.min_deletions[assignment_params.v] + assignment_params.ndV1;

      //cout<<"\t\tV deletion:vd_1"<<endl;
      //% This is a violation of a constraint.
      //% Happens for eg. when min_deletions is 0 and I consider 'negative' deletions.
      if( assignment_params.ndV < 0 || assignment_params.ndV > ((signed)_model.max_V_deletions))
	{
	  //                    % Go to next iteration of loop.
	  continue;
	  //        end
	}
      //% Calculate number of errors in V section, by leaving out errors that are in deleted region for this assignment.
      //% Also leaves out errors in left most 'high_error_region'.
      //cout<<"\t\tV deletion:vd_1"<<endl;          
      if( _ignore_deep_error)
	{
	  assignment_params.high_error_region = _model.high_error_region;//max(model.high_error_region, V.align_length(v) - ndV1 - deep_error_limit);
	}
      else
	assignment_params.high_error_region = _model.high_error_region;

      
      assignment_params.v_err_pos.clear();
      unsigned dim_size[]={_V.n_errors[assignment_params.v]};
      assignment_params.v_err_pos.initialize(1, dim_size, _V.error_positions[assignment_params.v]); //% error positions in V alignment
      dim_size[0]=3;
      Matrix<unsigned> v_err_excess_pos(1, dim_size,_V.excess_error_positions[assignment_params.v]); //% error positions in extended V alignment for 'negative' deletions, array format so far, not Matrix format yet
//cout<<"\t\tV deletion:vd_1"<<endl;       
      if(assignment_params.ndV1 < 0)
	{
	  //% Case of 'negative' deletions
	  assignment_params.v_ex_errs_i = (v_err_excess_pos <= (_V.align_position[assignment_params.v][0]+_V.align_length[assignment_params.v] -1 - assignment_params.ndV1) ) & (v_err_excess_pos > 0 );
	  assignment_params.v_ex_errs = sum_all_bool(assignment_params.v_ex_errs_i);
	  if (assignment_params.v_ex_errs > 1)//, % If the nts beyond the 1st -ve nt don't match, just skip.
	    {	    
	      continue;
	    }
	  assignment_params.nerrorsv = _V.n_errors[assignment_params.v] + assignment_params.v_ex_errs;// - sum(v_err_pos <= high_error_region & v_err_pos > 0);
	  //here we change the code, since we don't have to account for high_error_region
	}
      else
	{//% Case of positive deletions
	  assignment_params.v_ex_errs=0;
	  assignment_params.v_ex_errs_i.clear();//remove everything
	  if(_V.n_errors[assignment_params.v]==0)
	    assignment_params.nerrorsv=0;
	  else
	    assignment_params.nerrorsv = _V.n_errors[assignment_params.v] - sum_all_bool(assignment_params.v_err_pos > (_V.align_position[assignment_params.v][0]+_V.align_length[assignment_params.v]-1 - assignment_params.ndV1));// - sum_all_bool(assignment_params.v_err_pos > 0);
	  //            here change the very last part of the sum to not including
	  //the high error region. since we know we don't have this.
	}
                
      if( _no_error && assignment_params.nerrorsv > 0)
                    continue;
                           
      double log_perrv = assignment_params.nerrorsv*assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
                
      if (assignment_params.ndV==0 || ~assignment_params.assume_palindrome_negative_deletions)
	{
	  assignment_params.npV_max = _model.max_palindrome;
	  if(assignment_params.npV_max > _V.p_region_max_length[assignment_params.v][ assignment_params.ndV] ) //<---here we change the code it was assignment_params.ndV+1 in the second index
	    {
	      assignment_params.npV_max=_V.p_region_max_length[assignment_params.v][ assignment_params.ndV]; //<---here we change the code it was assignment_params.ndV+1 in the second index
	    }
	  assignment_params.npV_potential_max = assignment_params.npV_max;
	}
      else
	{
	  assignment_params.npV_max = 0;
	  assignment_params.npV_potential_max = _model.max_palindrome;
	  if(assignment_params.npV_potential_max> _V.p_region_max_length[assignment_params.v] [assignment_params.ndV]);//<---here we change the code it was assignment_params.ndV+1 in the second index
	  {
	    assignment_params.npV_potential_max=_V.p_region_max_length[assignment_params.v][assignment_params.ndV]; //<---here we change the code it was assignment_params.ndV+1 in the second index
	  }
	}                

      //cout<<"\t\tV deletion:vd_2"<<endl;
      //Upper bound on insertions nt bias factor
      int niVD_DJ_min= assignment_params.niVD_DJ0+assignment_params.ndV1-_D.align_length[assignment_params.d][0] - assignment_params.npV_max - assignment_params.p_max_J-assignment_params.p_max_Dl-assignment_params.p_max_Dr;
      
      if(niVD_DJ_min <0)
	{
	  niVD_DJ_min=0;
	}
      double log_p_max_nt_VD_DJ_dV = niVD_DJ_min*assignment_params.log_max_model_p_nt;
                //cout<<"\t\tV deletion:vd_3"<<endl;
      assignment_params.log_highest_probability_GIVEN_current_V_deletions = -1000;// % Highest probability so far GIVEN current V deletions AND outer loop variables {V, J, D alleles}                
      //cout<<"\t\tV deletion:vd_4"<<endl;
      double test_proba = assignment_params.log_probabase + log_perrv + (assignment_params.log_max_pcutV_loop_v+assignment_params.log_max_pcutJ_loop_j+assignment_params.log_max_pcutD_loop_d) + assignment_params.log_max_model_pins + log_p_max_nt_VD_DJ_dV;               
                
      if (test_proba < assignment_params.log_probability_threshold_factor + assignment_params.log_highest_probability)
	{
                    assignment_params.skips  ++;
                    continue;
	}
                
      assignment_params.V_align_length = _V.align_length[assignment_params.v] - assignment_params.ndV1 - _READ_LENGTH_CORRECTION;
      //this is an extra checking to make sure the alignment is not too long,
      //this is not in the original matlab code. it is based on the possible
      //V segment length
      if(assignment_params.V_align_length>=350)
	{
	  //assignment_params.v_break_out=true;
	  continue; //do the next one, this one is not right, next one here means next ndV1, not next v alignment
	  //return true; //we return true here, since we have set up the v_breakout
	  //here we just need to jump out to next V_allele and don't do this 
	}
//cout<<"\t\tV deletion:vd_6"<<endl;

      //======now start doing j deletions loop for
      //% Loop over J right deletions.
      //% ndJ1 is the deviation from the no. of deletions implied by the alignment (min_deletions).
      //% Lower bound is nd_start (usually -3).
      //% Upper bound is set by maximum deletions allowed and by minimum match length required
      signed temp_array[]={(signed)_model.max_excess_J_deletions,(signed)(_J.align_length[assignment_params.j]-_model.min_J_assign_length) ,(signed) (_model.max_J_deletions - _J.min_deletions[assignment_params.j])};
      cout<<"--->max_excess_J_deletions:"<<_model.max_excess_J_deletions<<endl;
      cout<<"--->min_J_assign_length:"<<_model.min_J_assign_length<<endl;
      cout<<"---->max_J_deletions:"<<_model.max_J_deletions<<endl;
      cout<<"_J.min_deletions"<<_J.min_deletions[assignment_params.j]<<endl;
      unsigned ndJ_max = min_mf(temp_array, 3);
      for( assignment_params.ndJ1=assignment_params.nd_start;assignment_params.ndJ1<(signed)ndJ_max;assignment_params.ndJ1++)
	{	//		%ndJ1=nd_start
	  cout<<"\t\t\t\t\tloop j deletion"<<assignment_params.ndJ1<<endl;
	  if (assignment_params.skips > assignment_params.max_skips)
	    {
	      _assigns.n_assignments = assignment_params.in;
	      _assigns.skips = assignment_params.skips;
	      assignment_params.v_break_out=true;
	      return false;//reason we report false is that we want the outside caller to jump out too. since we are doing too many skips
	    }
	  //% actual number of J deletions
	  assignment_params.ndJ = _J.min_deletions[assignment_params.j] + assignment_params.ndJ1;
	  //% This is a violation of a constraint.
	  //% Happens for eg. when min_deletions is 0 and I consider 'negative' deletions.
	  //cout<<"ndJ:"<<assignment_params.ndJ<<";model.max_J_deletions:"<<_model.max_J_deletions<<endl;
	  if (assignment_params.ndJ < 0 || assignment_params.ndJ >(signed) _model.max_J_deletions)
	    {
	      //% Go to next iteration of loop.
	      cout<<"*********CONTINUE:"<<endl;
	      
	      continue;
	    }
	  //cout<<"\t\t\tJ deletion:J_1"<<endl;
	  //niVD_DJ_total_with_ps = niVD_DJ0 + ndV1 + ndJ1;
          //          % Violation of constraint. V and J overlap in this assignment. So move on until this is not true.
	  if (assignment_params.niVD_DJ0+ assignment_params.ndV1+assignment_params.ndJ1 < 0)
	    continue;
                                      
	  //% Calculate number of errors in J section, by leaving out errors that are in deleted region for this assignment.
	  //count<<"\t\t\tJ deletion:J_2"<<endl;
	  dim_size[0]=_J.n_errors[assignment_params.j];
	  //count<<"dim_size[0]:"<<dim_size[0]<<endl;
	  
	  Matrix<unsigned> j_err_pos(1, dim_size, _J.error_positions[assignment_params.j]); //% error positions in J alignment
	  dim_size[0]=3;
	  //count<<"j_err_pos:dim:"<<j_err_pos.dim()<<endl;
	  //count<<"dims-zie[0]:"<<dim_size[0]<<endl;
	  Matrix<unsigned> j_err_excess_pos(1, dim_size, _J.excess_error_positions[assignment_params.j]);// % error positions in extended J alignment for 'negative' deletions
	  //count<<"\t\t\tJ deletion:J_3"<<endl;
	  //cout<<"assignment_params.ndJ1:"<<assignment_params.ndJ1<<endl;
	  if (assignment_params.ndJ1 < 0)
	    {
	      //% Case of 'negative' deletions
	      assignment_params.j_ex_errs_i = j_err_excess_pos >= (_J.align_position[assignment_params.j][0] + assignment_params.ndJ1);
	      assignment_params.j_ex_errs = sum_all_bool(assignment_params.j_ex_errs_i);
	      if (assignment_params.j_ex_errs > 1)
		{
		  cout<<"******GO NEXT J since this is negative deletions, but with too many errors in the \"negative\" zone. "<<endl;
		  continue;
		}
	      assignment_params.nerrorsj = _J.n_errors[assignment_params.j] + assignment_params.j_ex_errs;
	    }
	  else
	    {//                        % Case of positive deletions
	      assignment_params.j_ex_errs = 0;
	      assignment_params.j_ex_errs_i.clear();
 //count<<"\t\t\tJ deletion:J_4.1"<<endl;
	      if(_J.n_errors[assignment_params.j]==0)
		assignment_params.nerrorsj=0;
	      else
		assignment_params.nerrorsj = sum_all_bool( j_err_pos >= (_J.align_position[assignment_params.j][0] + assignment_params.ndJ1));
	    }//      end
                     //count<<"\t\t\tJ deletion:J_4.2"<<endl;
	  if (_no_error && assignment_params.nerrorsj > 0)
	    continue;
	            //count<<"\t\t\tJ deletion:J_4.3"<<endl;           
	  if( assignment_params.nerrorsj > assignment_params.J_max_error || ((double)assignment_params.nerrorsj)/(_J.align_length[assignment_params.j]-assignment_params.ndJ1) > 0.3)
	    {
	      //					disp('in the condition')
	      continue;
	    }
                    //count<<"\t\t\tJ deletion:J_4"<<endl;
                    
	  double log_perrj = assignment_params.nerrorsj*assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
          
	  //%                     % Upper bound on insertions nt bias factor
	  niVD_DJ_min =  assignment_params.niVD_DJ0 + assignment_params.ndV1 + assignment_params.ndJ1 - _D.align_length[assignment_params.d][0] - _V.p_region_max_length[assignment_params.v][1 + assignment_params.ndV] - _J.p_region_max_length[assignment_params.j][1 + assignment_params.ndJ] - assignment_params.p_max_Dl - assignment_params.p_max_Dr;
	  if(niVD_DJ_min<0)
	    {
	      niVD_DJ_min=0;
	    }
	  double log_p_max_nt_VD_DJ_dJ = niVD_DJ_min*assignment_params.log_max_model_p_nt;
          
	  assignment_params.log_highest_probability_GIVEN_current_J_deletions = -1000;                   //cout<<"\t\t\tJ deletion:J_5"<<endl;
	  
	  double test_proba = assignment_params.log_probabase + log_perrv + log_perrj + (assignment_params.log_max_pcutV_loop_v+assignment_params.log_max_pcutJ_loop_j+assignment_params.log_max_pcutD_loop_d) + assignment_params.log_max_model_pins + log_p_max_nt_VD_DJ_dJ;
	            
	  if( test_proba < assignment_params.log_probability_threshold_factor + assignment_params.log_highest_probability)
	    continue;
	       
                    //cout<<"\t\t\tJ deletion:J_6"<<endl;
	  // % Maximum possible half-lengths of palindromes, given current deletions values
	  if (assignment_params.ndJ==0 || ~assignment_params.assume_palindrome_negative_deletions)
	    {
	      assignment_params.npJ_max = _model.max_palindrome;
	      if(assignment_params.npJ_max> _J.p_region_max_length[assignment_params.j][assignment_params.ndJ])  //<---here we change the code it was assignment_params.ndj+1 in the second index
		
		assignment_params.npJ_potential_max = assignment_params.npJ_max;
	    }
	  else
	    {
              assignment_params.npJ_max = 0;
	      assignment_params.npJ_potential_max = _model.max_palindrome; 
	      if(assignment_params.npJ_potential_max>_J.p_region_max_length[assignment_params.j][ assignment_params.ndJ]) //<---here we change the code it was assignment_params.ndj+1 in the second index
		assignment_params.npJ_potential_max=_J.p_region_max_length[assignment_params.j][ assignment_params.ndJ]; //<---here we change the code it was assignment_params.ndj+1 in the second index
		
	    }//             end*/
	  //        cout<<"\t\t\tJ deletion:J_7"<<endl;
	  
	  assignment_params.J_align_length = _J.align_length[assignment_params.j] - assignment_params.ndJ1;
	  //cout<<"\t\t\tJ deletion:J_9"<<endl;
	  //% Loop over half-length of V palindrome
	  //for now it is just a place holder for the following code
	  //it simply return true and hopefully this will work for now.
	  if(!assign_VJ_palindrome(_model, _seq, _V, _D, _J,
				   _genV, _numV, _genD, _numD,
				   _genJ, _numJ,
				   _no_error, /*_ignore_deep_error,*/ _do_smoothing,
				   _force_all_alleles, _READ_LENGTH_CORRECTION,
				   assignment_params, _assigns)
	     )
	    {
	      /*if(assignment_params.v_break_out||assignment_params.j_break_out||assignment_params.d_break_out)
		{
		  return false;
		}
	      else
		{
		  break;
		  }*/
	      return false;
	    }

	  if(assignment_params.v_break_out||assignment_params.j_break_out||assignment_params.d_break_out)
	    {
	      break;
	    }

	  //here in ending the J deletions, we need to check for exit/breakingup cases
	  //  =======>  a little different, we add 2 to ndJ1>2 condition, instead of 1.
	  if(assignment_params.log_highest_probability_GIVEN_current_J_deletions>-900  &&
	     ( ( (assignment_params.ndJ1>2 && assignment_params.nerrorsj==0)&&(assignment_params.log_highest_probability_GIVEN_current_J_deletions < (assignment_params.log_highest_probability+assignment_params.log_probability_threshold_factor))
		 )
		||
	       ( (assignment_params.ndJ1>0)&&(assignment_params.log_highest_probability_GIVEN_current_J_deletions < (assignment_params.log_highest_probability+assignment_params.log_probability_hopeless_threshold_factor))
		 )
	       )
	     )
	    {
		break;
	    }
	  
	}//end of J deletion loop

      if(assignment_params.v_break_out||assignment_params.j_break_out||assignment_params.d_break_out)
	{
	  break;
	  //return false;
	}
      
      //here in ending the J deletions, we need to check for exit/breakingup cases
      //  =======>  a little different, we add 2 to ndJ1>2 condition, instead of 1.
      if(assignment_params.log_highest_probability_GIVEN_current_V_deletions>-900  &&
	 ( ( (assignment_params.ndV1>2 && assignment_params.nerrorsv==0)&&(assignment_params.log_highest_probability_GIVEN_current_V_deletions < (assignment_params.log_highest_probability+assignment_params.log_probability_threshold_factor))
	     )
	   ||
	   ( (assignment_params.ndV1>0)&&(assignment_params.log_highest_probability_GIVEN_current_V_deletions < (assignment_params.log_highest_probability+assignment_params.log_probability_hopeless_threshold_factor))
	     )
	   )
	 )
	{
	  break;
	}
      
    }//end of V deletions loop for
  
  return true;
}//end of assign_V_deletions

//=======================================
//***************************************
/*output:
 *     returning false for 3 cases, 1) npV and 2)npJ for maximum skips and assignments
 *         3) for returning false from inner function of D_alignment
 *       upon these cases, we will return false to upper.
 *     other cases, we will return true;
 */
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
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns)
{
  //set the direction for VJ_palindrome search, IS THIS REALLY NECESSARY or IMPORTANT???
  //do it anyway for now
  //assignment_params.nd_start=0;
  //assignment_params.np_start_from_max=true;<=====

  //cout<<"inside: VJ palindrome"<<endl;
  int np_step, npV_start, npV_end;
  int npJ_start, npJ_end;
  if (assignment_params.np_start_from_max)
    {
      np_step = -1;
      npV_start = assignment_params.npV_max;
      npV_end = 0;
      npJ_start = assignment_params.npJ_max;
      npJ_end = 0;
    }
  else
    {
      np_step = 1;
      npV_start = 0;
      npV_end = assignment_params.npV_max;
      npJ_start = 0;
      npJ_end = assignment_params.npJ_max;
    }//       end
  //cout<<" np_step:"<<np_step<<";npV_start:"<<npV_start<<";npV_end:"<<npV_end<<endl;
  //now looping the possible cases of npV first
  for(assignment_params.npV=npV_start;;assignment_params.npV+=np_step)
    {
      cout<<"\t\t\t\t\t\tlooping: assignment_params.npV:"<<assignment_params.npV<<endl;
      if(assignment_params.np_start_from_max&&(signed)assignment_params.npV<npV_end)
	{
	  break; //we are done
	}
      if(!assignment_params.np_start_from_max&&(signed)assignment_params.npV>npV_end)
	{
	  break;//we are done.
	}
      
      if (assignment_params.skips > assignment_params.max_skips)
	{
	  _assigns.n_assignments = assignment_params.in;
	  _assigns.skips =  assignment_params.skips;
	  //assignment_params.v_break_out=true;//===in this case, we want to jump all out.
	  return false;
	}//       end
      
      assignment_params.V_end=_V.align_position[assignment_params.v][0]+_V.align_length[assignment_params.v]-1-assignment_params.ndV1+assignment_params.npV;
      if(((unsigned)assignment_params.V_end)>=((unsigned)_J.align_position[assignment_params.j][0]+assignment_params.ndJ1))
	//V with palindrome overlaps with J, not possible
	{
	  continue;
	}
      //cout<<"888888inside vj palindrome"<<endl;
      //now start doing the cut variable, Note: cut variable is the observed deletion since P_nt often is confused with negative deletions.
      int ncutV=assignment_params.ndV-assignment_params.npV; //ncutV could be negative, that is why below, we need to add min_V_cut (which is negative value) to make sure it is positive.  the ncutV will never be smaller than -min_V_cut, which is negative model_max_palindrome (-6)
      double PcutV=_model.PcutV_given_V(assignment_params.v_g, ncutV-_model.min_V_cut);

      //cout<<"ncutV: "<<ncutV<<";PcutV="<<PcutV<<endl;
      //cout<<"v_g:"<<assignment_params.v
      /*if(assignment_params.npV==0&&assignment_params.v==0)
	{
	  cout<<"pcutv_given_v matrix:"<<_model.PcutV_given_V.toString()<<endl;
	  }*/
      if(PcutV==0)
	{
	  assignment_params.skips++;
	  continue;
	}
      assignment_params.log_PcutV=log(PcutV);
      double log_max_pcutVDJ_loop_pV=assignment_params.log_PcutV+assignment_params.log_max_pcutD_loop_d+assignment_params.log_max_pcutJ_loop_j;
      //cout<<">>>>>>>>>>inside vj palindrome"<<endl;
      //upper bound on insertions nt bias factor
      int niVD_DJ_min=(signed)(assignment_params.niVD_DJ0)+assignment_params.ndV1+assignment_params.ndJ1-_D.align_length[assignment_params.d][0]-assignment_params.npV-_J.p_region_max_length[assignment_params.j][assignment_params.ndJ]-assignment_params.p_max_Dl-assignment_params.p_max_Dr;
      if( niVD_DJ_min<0)
	{
	  niVD_DJ_min=0;
	}
      double log_p_max_nt_VD_DJ_pV=niVD_DJ_min*assignment_params.log_max_model_p_nt;
      
      double log_perrj=assignment_params.nerrorsj* assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
      double log_perrv=assignment_params.nerrorsv* assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
      double test_proba=
	assignment_params.log_probabase + log_perrj
	+ log_perrv+ log_max_pcutVDJ_loop_pV+
	assignment_params.log_max_model_pins + log_p_max_nt_VD_DJ_pV;
      if(test_proba< assignment_params.log_probability_threshold_factor+
	 assignment_params.log_highest_probability)
	{
	  assignment_params.skips++;
	  continue;
	}
	 
      //cout<<"start doing loop over j aplindrome"<<endl;
      //next start doing the loop over half-length of J palindrome
      for(assignment_params.npJ=npJ_start;;assignment_params.npJ+=np_step)
	{
	  cout<<"\t\t\t\t\t\t\tloop j palindrome: "<<assignment_params.npJ<<endl;
	  if(assignment_params.np_start_from_max&&(signed)assignment_params.npJ<npJ_end)
	    {
	      break; //we are done
	    }
	  if(!assignment_params.np_start_from_max&&(signed)assignment_params.npJ>npJ_end)
	    {
	      break;//we are done.
	    }
	  
	  if (assignment_params.skips > assignment_params.max_skips)
	    {
	      _assigns.n_assignments = assignment_params.in;
	      _assigns.skips =  assignment_params.skips;
	      //assignment_params.v_break_out=true;//===in this case, we want to jump all out.
	      return false;
	    }//       end
	  
	  assignment_params.J_start=_J.align_position[assignment_params.j][0]+assignment_params.ndJ1-assignment_params.npJ;

	  if(((unsigned)assignment_params.V_end)>=(unsigned)assignment_params.J_start)
	    //V with palindrome overlaps with J, not possible
	    {
	      continue;
	    }
	  
	  //now cut variable from J side
	  int ncutJ=assignment_params.ndJ-assignment_params.npJ;
	  double PcutJ=_model.PcutJ_given_J(assignment_params.j_g, ncutJ-_model.min_J_cut);
	  if(PcutV==0)
	    {
	      assignment_params.skips++;
	      continue;
	    }
	  assignment_params.log_PcutJ=log(PcutJ);
	  double log_max_pcutVDJ_loop_pJ=
	    assignment_params.log_PcutV+
	    assignment_params.log_PcutJ+assignment_params.log_max_pcutD_loop_d;
      
	  //start doing something new here for D 
	  //check for the valid D length and valid cut
	  assignment_params.n_D_aligns=_D.numOfAligned[assignment_params.d];
	  assignment_params.start_n_D_aligns=0;
	  
	  //find the first completely valid alignment
	  bool completely_valid_na=false;
	  unsigned c_na=0;
	  cout<<"starting doing the valid d na while loop........."<<endl;
	  while(!completely_valid_na&&c_na<assignment_params.n_D_aligns)
	    {
	      completely_valid_na=((signed)_D.align_position_left[assignment_params.d][c_na]>assignment_params.V_end) 
		&&((signed)_D.align_position_right[assignment_params.d][c_na]<assignment_params.J_start);
	      c_na++;
	    }//end of while for completely_valid_na;
	  unsigned first_valid_length=0;

	  if(completely_valid_na)
	    {
	      first_valid_length=_D.align_length[assignment_params.d][c_na-1];
	    }
	  else
	    {
	      //find the first partly valid alignment
	      bool partly_valid_na=false;
	      unsigned p_na=0;
	      while(~partly_valid_na&&p_na<assignment_params.n_D_aligns)
		{
		  partly_valid_na=((signed)_D.align_position_left[assignment_params.d][p_na]<assignment_params.J_start) 
		    ||((signed) _D.align_position_right[assignment_params.d][p_na]>assignment_params.V_end);
		  p_na++;
		}//end of partly_valid_na while loop

	      if(partly_valid_na)
		{
		  p_na--;//because it is increamented before exiting
		  int p_valid_start=((assignment_params.J_start-1)>((signed)_D.align_position_right[assignment_params.d] [p_na]))?_D.align_position_right[assignment_params.d] [p_na]:(assignment_params.J_start-1);
		  int p_valid_end=((assignment_params.V_end+1)>(signed)_D.align_position_left[assignment_params.d][ p_na])?(assignment_params.V_end+1):_D.align_position_left[assignment_params.d][p_na];
		  first_valid_length=p_valid_start-p_valid_end+1;
		  
		  assignment_params.start_n_D_aligns=p_na;
		}
	      else //no partially valid alignment
		{
		  first_valid_length=0;
		  //do only zero D. no even partly valid alignment.
		  assignment_params.start_n_D_aligns=assignment_params.n_D_aligns+1;//set it to be out of range, so for the later loop, we won't do anything, simply stop before the first iteration.
		}//end of partly valid else case;

	    }//end of completely_valid_na else case
	  
	  if(assignment_params.best_D_align_length<first_valid_length)
	    {
	      assignment_params.best_D_align_length=first_valid_length;
	    }
	  if(assignment_params.best_D_align_length<=8)
	    {
	      assignment_params.end_n_D_aligns=assignment_params.n_D_aligns+1;
	    }
	  else
	    {
	      assignment_params.end_n_D_aligns=assignment_params.n_D_aligns;
	    }

	  //NDN length
	  assignment_params.niVD_DJ_total=assignment_params.J_start-assignment_params.V_end-1;
	  //upper bound for nt bias factor
	  int niVD_DJ_min = (signed)assignment_params.niVD_DJ_total-first_valid_length;
	  if(assignment_params.start_n_D_aligns<=assignment_params.n_D_aligns)  //there is valid length
	    {
	      niVD_DJ_min=niVD_DJ_min-
	    _D.p_region_max_length_left[assignment_params.d][ assignment_params.start_n_D_aligns][ _D.deletions_left[assignment_params.d] [assignment_params.start_n_D_aligns]]-
		_D.p_region_max_length_right[assignment_params.d][ assignment_params.start_n_D_aligns] [_D.deletions_right[assignment_params.d][ assignment_params.start_n_D_aligns]];
	    }
	  if(niVD_DJ_min<0)
	    niVD_DJ_min=0;
			
	  double log_p_max_nt_VD_DJ_pJ=niVD_DJ_min*assignment_params.log_max_model_p_nt;
	  
	  test_proba=assignment_params.log_probabase+
	    log_perrv+log_perrj+log_max_pcutVDJ_loop_pJ+
	    assignment_params.log_max_model_pins+log_p_max_nt_VD_DJ_pJ;
	  
	  if(test_proba< assignment_params.log_probability_threshold_factor+
	     assignment_params.log_highest_probability)
	    {
	      assignment_params.skips++;
	      continue;
	    } 
	  //cout<<"doing assign D function..........."<<endl;
	  //===========here is page 12 on matlab code print out.
	  //---------> start doing d alignment from here
	  //====>will add the function to do the D alignment here
	  if(!assign_D(_model, _seq, _V, _D, _J, _genV, _numV, _genD, _numD, _genJ, 
		      _numJ, /*const double& _probability_threshold_factor,*/ _no_error,
		      /*const bool& _ignore_deep_error,*/ _do_smoothing, 
		      _force_all_alleles, _READ_LENGTH_CORRECTION,
		      /*output, input*/assignment_params,
		      /*output*/_assigns)
	     )
	    {
	      return false;//because in the case of returning false from assign_D
	      //we want to finish whatever we have for all loopings and don't
	      //do anything. this is the inner most level. the deeper level than
	      //this is the run_stat functions, where we don't do this checking (skips/assignments;
	    }
	  
	  
	  if(assignment_params.v_break_out||assignment_params.d_break_out||assignment_params.j_break_out)
	    {
	      break;
	    }
	}//end of npJ loop 
      if(assignment_params.v_break_out||assignment_params.d_break_out||assignment_params.j_break_out)
	{
	  break;
	}
    }//end of npV loop palindrome

  /*
  ===>>  //Note for myself, we probably we include the blow code inside the above blocks of for loops
  
  // ***************************
  //===========>Good code below from previous version, need to be carefully
  //in here we start doing the stats!! set things to the assigns.

  //Increment valid assignment number
  //assignment_params.in+=1;//disp(['assignment=', num2str(in)]) <==========here we move the increment to the end,
  //assignment_params.zeroD=false;
  unsigned v=assignment_params.v;
  unsigned j=assignment_params.j;
  unsigned d=assignment_params.d;
  unsigned na=0;//this is the n of aligned d seg index, pointing to the longest one for one
  //in the case of the best one is not valide we just go ahed to through it away for now

  //count<<"\t\t\tvdJ palindrome:J_1"<<endl;
  //insertions <--------added by Feng
  unsigned inVD=_D.align_position_left[d][na]-
    (_V.align_position[v][0]+_V.align_length[v]-1)-assignment_params.npV-assignment_params.npDl;
  //cout<<"\t\t\tvdj palidnrome:J_1.1"<<endl;
  if((signed)inVD<0)
    inVD=0;
  //cout<<"\t\t\tvdj palidnrome:J_1.2"<<endl;
  
  unsigned inDJ=_J.align_position[j][0]-_D.align_position_right[d][na]-assignment_params.npDr-assignment_params.npJ;
  if((signed)inDJ<0)
    {
      inDJ=0;
    }
//count<<"\t\t\tvdj palidnrome:J_1.3"<<endl;
 if(inVD>=_model.max_insertions||inDJ>=_model.max_insertions)
   return true;//continue;
  //% Compute final probability of assignment by multiplying base probability by factors for
  //% palindromes, deletions and the number of errors.
  //% multiply probability by factor for insertions of given lengths and sequence specific factors for
  //% nucleotide bias
 //cout<<"inVD:"<<inVD<<";inDJ:"<<inDJ<<endl;
 //cout<<"PinsVD:"<<_model.PinsVD.toString()<<endl;
 //cout<<"PinsVD(inVD):"<< _model.PinsVD(inVD)<<"PinsDJ(inDJ):"<<_model.PinsDJ(inDJ)<<endl;
 double log_pins;
 if(_model.PinsVD(inVD)*_model.PinsDJ(inDJ)==0)
   {
     log_pins=0;
   }
 else
   {
   log_pins= log(_model.PinsVD(inVD)*_model.PinsDJ(inDJ));//<=== be careful here, check matlab code
                             //cout<<"\t\t\tvdj palidnrome:J_2"<<endl;                       
   }
  //double log_pntbias_VD = log_RnucleotideVD_per_nucleotideVD_5prime*nucleotideVD(:);
  //double log_pntbias_DJ = log_RnucleotideDJ_per_nucleotideDJ_3prime*nucleotideDJ(:);
  assignment_params.nerrorsd=0;//set it to zero for now<-============
  //cout<<"e v:"<<assignment_params.nerrorsv<<"e j:"<< assignment_params.nerrorsj<<";e d:"<< assignment_params.nerrorsd<<endl;
   unsigned nerrors = assignment_params.nerrorsv + assignment_params.nerrorsj + assignment_params.nerrorsd;
   //% Cost for errors is (error_rate/3)^(n_errors) because we know the mistaken nucleotides.
   //count<<"nerrors:"<<nerrors<<endl;
   double log_perr;
   if (nerrors==0)
     log_perr = 0;
   else
     log_perr = nerrors*assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
   
                             ////cout<<"\t\t\tvdj palidnrome:J_3"<<endl;
   cout<<"log_probabase:"<<assignment_params.log_probabase<<";log_perr:"<< log_perr <<";log_pins:"<<    log_pins<<endl;        
   double log_proba = assignment_params.log_probabase + log_perr + / *log_PcutVDJ +* / 
    log_pins;// + log_pntbias_DJ + log_pntbias_VD;
                                                    
                                                    
  //% Update all highest probability so far values, if necessary
  
  
  if(log_proba > assignment_params.log_highest_probability_GIVEN_current_D_allele)
    assignment_params.log_highest_probability_GIVEN_current_D_allele = log_proba; //end;
  if(log_proba > assignment_params.log_highest_probability_GIVEN_current_V_allele)
    assignment_params.log_highest_probability_GIVEN_current_V_allele = log_proba; //end;
  if(log_proba > assignment_params.log_highest_probability_GIVEN_current_J_allele)
    assignment_params.log_highest_probability_GIVEN_current_J_allele = log_proba; //end;
  if(log_proba > assignment_params.log_highest_probability_GIVEN_current_V_deletions)
    assignment_params.log_highest_probability_GIVEN_current_V_deletions = log_proba; //end;
  if(log_proba > assignment_params.log_highest_probability_GIVEN_current_J_deletions)
    assignment_params.log_highest_probability_GIVEN_current_J_deletions = log_proba; //end;
  / *if(log_proba > log_highest_probability_GIVEN_current_Dl_deletions) 
    log_highest_probability_GIVEN_current_Dl_deletions = log_proba; //end;
  if(log_proba > log_highest_probability_GIVEN_current_Dr_deletions)
    log_highest_probability_GIVEN_current_Dr_deletions = log_proba; //end;
  * /
  //cout<<"\t\t\tvdj palidnrome:J_4"<<endl;
  //cout<<"log_proba"<<log_proba<<";log_probability_threshold_factor:"<<assignment_params.log_probability_threshold_factor<<";log_highest_probability:"<<assignment_params.log_highest_probability<<endl;
  if (assignment_params.D_align_length > assignment_params.best_D_align_length) 
    assignment_params.best_D_align_length = assignment_params.D_align_length; //end
                                                    
  //% if D length isvery short (<3) I let it count low prob events because the total prob of D short
  //% is spread out over many choices of D deletions. So if I skip them,
  //%if (~do_zeroD && (log_proba < log_probability_threshold_factor + log_highest_probability)) || (log_proba < log_probability_hopeless_threshold_factor + log_highest_probability)
  if (((log_proba < assignment_params.log_probability_threshold_factor + assignment_params.log_highest_probability)) || (log_proba < assignment_params.log_probability_hopeless_threshold_factor + assignment_params.log_highest_probability))
    {
      assignment_params.skips ++;
      //continue;
      return true;//go next round
    }//end
  unsigned ndDl1=0, ndDr1=0;
  //cout<<"\t\t\tvdj palidnrome:J_5"<<endl;
  if (~assignment_params.zeroD)
    {
      assignment_params.genic_length = 
	(_V.align_length[v] + _J.align_length[j] + 
	 assignment_params.npV + assignment_params.npJ -assignment_params.ndV1 
	 -assignment_params.ndJ1 + _D.align_length[d][na] + assignment_params.npDl + 
	 assignment_params.npDr - ndDl1 - ndDr1 
	 - assignment_params.high_error_region - _READ_LENGTH_CORRECTION);
    }
  else
    {
      assignment_params.genic_length = 
	(_V.align_length[v] + _J.align_length[j] 
	 +assignment_params.npV + assignment_params.npJ - 
	 assignment_params.ndV1 - assignment_params.ndJ1 - 
	 assignment_params.high_error_region - _READ_LENGTH_CORRECTION);
    }//end
  assignment_params.D_align_length=_D.align_length[d][na]-ndDl1-ndDr1;
  
  //cout<<"vdj palindrome "<<endl;
  assignment_params.n_assignments_v_gene ++;//  1; 
  assignment_params.n_assignments_d_gene ++;//= n_assignments_d_gene + 1;
  assignment_params.n_assignments_j_gene ++;//= n_assignments_j_gene + 1;
  if(log_proba > assignment_params.log_highest_probability)
    {
      assignment_params.log_highest_probability = log_proba; 
      _assigns.max_proba_index = assignment_params.in; //here because we increment the assignment_param.in before here 
    }//end;                                                    
  
  //v_err_pos_rel = zeros(max_V_depth,1);
  //j_err_pos_rel = zeros(max_J_depth,1);
  if (~_no_error) //we allow error in this case.
    {
      if(assignment_params.nerrorsv>0)
	{
	  //place holder
	}
      //% Set positions of errors in J to 1, accounting for deletions
      
      if(assignment_params.nerrorsj > 0)
	{
	  //plac holder
	}// end
      //% Set positions of errors in D to 1. All of these come from 'negative' deletions since D alignments do not
      //% allow errors.
                                                        
      if (~assignment_params.zeroD && assignment_params.nerrorsd > 0)
	{                           
	  / *if ndDl1 < 0 && d_ex_errs_left > 0
		     error_vs_position(d_err_excess_pos_left(d_ex_errs_left_i)) = 1;
	  end
	    
	    if ndDr1 < 0 && d_ex_errs_right > 0
		       error_vs_position(d_err_excess_pos_right(d_ex_errs_right_i)) = 1;
	  end
	    
	    if d_errs > 0
			  error_vs_position(d_err_pos(d_errs_i)) = 1;
			  end * /
	}//end
							  
    }//  end of _no_error if blaock
  
*/  
  return true;
}

//start doing assignment with D seg on both sides
//output: only two places that will give out false. 1)npDl and 2)npDr checking for 
//        maximum skips or maximum assignments
//        all other cases are returning true;

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
 )
{
  
  //start looping all the D alignments
  for(unsigned na=assignment_params.start_n_D_aligns;na<=assignment_params.end_n_D_aligns;
      na++)
    {
      cout<<"\t\t\t\t\t\t\tstart na loop for D:"<<na<<endl;
      bool zeroD=na>assignment_params.n_D_aligns;
      assignment_params.na=na;
      //this d_errs_interanal_max is a local variable used to set up break condistion
      //at the end of dr and dl deletions
      unsigned d_errs_internal_max=_D.n_errors[assignment_params.d][na];
      
      if(!zeroD && _D.align_length[assignment_params.d] [na]<assignment_params.best_D_align_length-8)
	return true;//???? this is same as break, since it is inside the outer most loop

      //if the whole D alignment lies within the V or J
      if(!zeroD && ((signed)_D.align_position_right[assignment_params.d][na]<= assignment_params.V_end || ((signed)_D.align_position_left[assignment_params.d][ na])>=assignment_params.J_start))
	{
	  continue;
	}
	
      //get ready for the deletion left loop
      int ndDl_min, ndDl_max;
      if(!zeroD) 
	{
	  ndDl_min=_D.deletions_left[assignment_params.d][na]+assignment_params.nd_start;
	  
	  int temp=_model.max_excess_D_deletions;
	  if(temp>(signed)(_model.max_D_deletions-_D.deletions_left[assignment_params.d][ na]))
	    temp=_model.max_D_deletions-_D.deletions_left[assignment_params.d][na];
	  if(temp> (signed)(_D.align_length[assignment_params.d] [na]-1))
	    temp=_D.align_length[assignment_params.d][na]-1;
	  
	  ndDl_max=_D.deletions_left[assignment_params.d][ na]+temp;
	  //d_errs_internal_max=_D.n_errors(na,d)
	}
      else
	{
	  ndDl_min=0;
	  
	  ndDl_max=_model.max_D_deletions;
	  if(ndDl_max>(signed)_genD[assignment_params.d_a].Get_Seq().GetLength())
	    {
	      ndDl_max=_genD[assignment_params.d_a].Get_Seq().GetLength();
	    }
	  d_errs_internal_max=0;
	}

      unsigned npDl_max, npDr_max; //local variables
      //start the deletion left loop
      for(assignment_params.ndDl=ndDl_min;assignment_params.ndDl<=ndDl_max;assignment_params.ndDl++)
	{
	  cout<<"\t\t\t\t\t\t\t\t\tstart ndDl loop :"<<assignment_params.ndDl<<endl;
	  if(assignment_params.ndDl<0 || assignment_params.ndDl>(signed)_model.max_D_deletions)
	    continue;

	  if(!zeroD)
	    {
	      assignment_params.ndDl1=assignment_params.ndDl-
		_D.deletions_left[assignment_params.d][na];
	      assignment_params.npDl_potential_max=_model.max_palindrome;
	      if(assignment_params.npDl_potential_max<
		 _D.p_region_max_length_left[assignment_params.d][assignment_params.ndDl][na])
		{
		  assignment_params.npDl_potential_max=
		    _D.p_region_max_length_left[assignment_params.d][assignment_params.ndDl][na];
		}
	    }
	  else
	    {
	      assignment_params.ndDl1=0;
	      assignment_params.npDl_potential_max=0;
	    }
		 
	  //only in this special case, we consider the p_nts
	  if(assignment_params.ndDl1==0)
	    {
	      npDl_max=assignment_params.npDl_potential_max;
	    }
	  else //in other case, no way to tell the p_nt from insertions
	    {
	      npDl_max=0;
	    }
	  
	  //errors
	  unsigned dim_size[]={3};
	  Matrix<unsigned> d_err_excess_pos_left
	    (1, dim_size, _D.excess_error_positions_left[assignment_params.d][na]);
	  Matrix<unsigned> d_err_excess_pos_right
	    (1, dim_size, _D.excess_error_positions_right[assignment_params.d][na]);
	  //now, we have the dim size about errors
	  if(assignment_params.ndDl1<0)  //negative deletions, where the deletion is "recorded" only because of sequencing error/mutations. 
	    {
	      //negative deletions
	      assignment_params.d_ex_errs_left_i=
		d_err_excess_pos_left>=(_D.align_position_left[assignment_params.d][na]+assignment_params.ndDl1);
	      assignment_params.d_ex_errs_left=sum_all_bool(assignment_params.d_ex_errs_left_i);
	      if(assignment_params.d_ex_errs_left>1)
		continue;
	      
	    }
	  else //positive deletions, deletions is actually longer, but because of insertion, it seems shorter.
	    {
	      assignment_params.d_ex_errs_left_i.clear();
	      assignment_params.d_ex_errs_left=0;
	    }
	  
	  //  }
	  //more !zero cases for doing deletions right
	  //and insertion nt distributions
	  unsigned ndDr_min, ndDr_max;
	  SequenceString insVD_nseq_min_loopdDl;
	  SequenceString insDJ_nseq_min_loopdDr;
	  double log_p_max_nt_VD_loop_Dl;
	  double log_p_max_nt_DJ_loop_Dr;
	  if(!zeroD)
	    {
	      //first get ready for right deletion on D
	      ndDr_min=_D.deletions_right[assignment_params.d][na]+assignment_params.nd_start;
	      ndDr_max=_model.max_excess_D_deletions;
	      if(ndDr_max>_model.max_D_deletions-_D.deletions_right[assignment_params.d][na])
		{
		  ndDr_max=_model.max_D_deletions-_D.deletions_right[assignment_params.d][na];
		}
	      if(ndDr_max>_D.align_length[assignment_params.d][na]-assignment_params.ndDl1-1)
		{
		  ndDr_max=_D.align_length[assignment_params.d][na]-assignment_params.ndDl1-1;
		}
	      ndDr_max+=_D.deletions_right[assignment_params.d][na];

	      //we need to do insertion nt calculation
	      //do some extra checking
	      if(assignment_params.V_end+1> (signed)(_D.align_position_left[assignment_params.d][na]+assignment_params.ndDl1-npDl_max-1))
		{
		  //it is possible in this case the ends cross, so we have to be careful
		  log_p_max_nt_VD_loop_Dl=0;
		  //within this "error" case, we did not set insVD_nseq_min_loopdDl sequence
		  //string, so it is a default empty one.
		}
	      else
		{
		  insVD_nseq_min_loopdDl=_seq.Sub(assignment_params.V_end+1, _D.align_position_left[assignment_params.d][na]+assignment_params.ndDl1-npDl_max-1);
		  if(insVD_nseq_min_loopdDl.GetLength()>_model.max_insertions)
		    {
		      continue;
		    }
		  //now get the statistics
		  log_p_max_nt_VD_loop_Dl=(insVD_nseq_min_loopdDl.GetLetterCount('A')+insVD_nseq_min_loopdDl.GetLetterCount('a'))*assignment_params.log_max_model_p_nt_VD(0);
		  log_p_max_nt_VD_loop_Dl+=(insVD_nseq_min_loopdDl.GetLetterCount('C')+insVD_nseq_min_loopdDl.GetLetterCount('c'))*assignment_params.log_max_model_p_nt_VD(1);
		  log_p_max_nt_VD_loop_Dl+=(insVD_nseq_min_loopdDl.GetLetterCount('G')+insVD_nseq_min_loopdDl.GetLetterCount('g'))*assignment_params.log_max_model_p_nt_VD(2);
		  log_p_max_nt_VD_loop_Dl+=(insVD_nseq_min_loopdDl.GetLetterCount('T')+insVD_nseq_min_loopdDl.GetLetterCount('t'))*assignment_params.log_max_model_p_nt_VD(3);
		}
	    }//end zeroD cases if loop
	  else
	    {
	      ndDr_min=_genD[assignment_params.d_a].Get_Seq().GetLength()-assignment_params.ndDl;
	      ndDr_max=ndDr_min;
	      log_p_max_nt_VD_loop_Dl=0;
	      
	    }//endo fo zeroD else loop
	  
	  assignment_params.log_highest_probability_GIVEN_current_Dl_deletions=-1000;
	  
	  //doing d deletions right
	  for(assignment_params.ndDr=ndDr_min;assignment_params.ndDr<=(signed)ndDr_max;assignment_params.ndDr++)
	    {
	      cout<<"\t\t\t\t\t\t\t\t\t\t\t\tstart doing ndDr loop:"<<assignment_params.ndDr<<endl;
	      if(assignment_params.ndDr<0 || assignment_params.ndDr>(signed)_model.max_D_deletions)
		{
		  continue;
		}
	      
	      if(!zeroD)
		{
		  assignment_params.npDr_potential_max=_model.max_palindrome;
		  if(assignment_params.npDr_potential_max<_D.p_region_max_length_right[assignment_params.d][na][assignment_params.ndDr])
		    {
		      assignment_params.npDr_potential_max=_D.p_region_max_length_right[assignment_params.d][na][assignment_params.ndDr];
		    }
		}
	      else
		{
		  assignment_params.npDr_potential_max=0;
		}

	      if(assignment_params.ndDr==0)
		{
		  npDr_max=assignment_params.npDr_potential_max;
		}
	      else
		{
		  npDr_max=0;
		}
	      
	      //doing nt distributions
	      double log_p_max_nt_VD_DJ_na_loop_Dr;
	      double log_max_pins_loop_dDr;
	      if(!zeroD)
		{
		  assignment_params.ndDr1=
		    assignment_params.ndDr-_D.deletions_right[assignment_params.d][na];
		  assignment_params.D_align_length=
		    _D.align_length[assignment_params.d][na]-assignment_params.ndDl1
		    -assignment_params.ndDr1;
		  
		  if((_D.align_position_right[assignment_params.d][na]-assignment_params.ndDr1+npDr_max+1) >= (unsigned)assignment_params.J_start-1)
		    {
		      //within this zero case, we set the sequencestring to be default empty
		      //string, so it has zero length
		      log_p_max_nt_DJ_loop_Dr=0;
		      log_max_pins_loop_dDr=log(_model.PinsDJ(0)*_model.PinsVD(insVD_nseq_min_loopdDl.GetLength()));
		    }
		  else
		    {

		      insDJ_nseq_min_loopdDr=_seq.Sub(_D.align_position_right[assignment_params.d][na]-assignment_params.ndDr1+npDr_max+1, assignment_params.J_start-1);
		      if(insDJ_nseq_min_loopdDr.GetLength()>_model.max_insertions)
			{
			  cout<<"BREAKOUT:maxi insertion:"<<_model.max_insertions<<";insDJ_nseq_min_loopDr:"<<insDJ_nseq_min_loopdDr.GetLength()<<endl;
			  continue;
			}
		  //now get the statistics
		      log_p_max_nt_DJ_loop_Dr=(insDJ_nseq_min_loopdDr.GetLetterCount('A')+insDJ_nseq_min_loopdDr.GetLetterCount('a'))*assignment_params.log_max_model_p_nt_DJ(0);
		      log_p_max_nt_DJ_loop_Dr+=(insDJ_nseq_min_loopdDr.GetLetterCount('C')+insDJ_nseq_min_loopdDr.GetLetterCount('c'))*assignment_params.log_max_model_p_nt_DJ(1);
		      log_p_max_nt_DJ_loop_Dr+=(insDJ_nseq_min_loopdDr.GetLetterCount('G')+insDJ_nseq_min_loopdDr.GetLetterCount('g'))*assignment_params.log_max_model_p_nt_DJ(2);
		      log_p_max_nt_DJ_loop_Dr+=(insDJ_nseq_min_loopdDr.GetLetterCount('T')+insDJ_nseq_min_loopdDr.GetLetterCount('t'))*assignment_params.log_max_model_p_nt_DJ(3);
		    }
		  log_p_max_nt_VD_DJ_na_loop_Dr=log_p_max_nt_DJ_loop_Dr+log_p_max_nt_VD_loop_Dl;

		  log_max_pins_loop_dDr=log(_model.PinsDJ(insDJ_nseq_min_loopdDr.GetLength())*_model.PinsVD(insVD_nseq_min_loopdDl.GetLength()));
		  
		  //calculate number of errors in D section, by leaving out errors that are in deleted region for this assignment.
		  dim_size[0]=_D.n_errors[assignment_params.d][na];
		  Matrix<unsigned> d_err_pos 
		    (1, dim_size, _D.error_positions[assignment_params.d][na]);
		  assignment_params.d_errs_i=
		    (d_err_pos>=_D.align_position_left[assignment_params.d][na]+assignment_params.ndDl1)&
		    (d_err_pos<=_D.align_position_right[assignment_params.d][na]-assignment_params.ndDr1);

		  assignment_params.d_errs=sum_all_bool(assignment_params.d_errs_i);
		}//end of non-zeroD case, for nt distribution 
	      else
		{
		  int niVD_DJ_min=assignment_params.niVD_DJ_total;
		  log_p_max_nt_VD_DJ_na_loop_Dr=niVD_DJ_min*assignment_params.log_max_model_p_nt;
		  log_max_pins_loop_dDr=assignment_params.log_max_model_pins;

		  assignment_params.ndDr1=0;
		  assignment_params.D_align_length=0;
		  assignment_params.d_errs=0;
		  assignment_params.d_errs_i.clear();
		  //d_err_pos.clear();
		}//end of zeroD else loop for nt distribution

	      //errors in the excess region
	      if(assignment_params.ndDr1<0)
		{
		  //error positions in extended j alignment for negative deletions
		  //are stored in d_err_excess_pos_right
		  //now get them for each 
		  assignment_params.d_ex_errs_right_i=
		    (d_err_excess_pos_right>0)&(d_err_excess_pos_right<=(_D.align_position_right[assignment_params.d][na]-assignment_params.ndDr1));
		  assignment_params.d_ex_errs_right=sum_all_bool(assignment_params.d_ex_errs_right_i);
		  if(assignment_params.d_ex_errs_right>1)
		    {
		      cout<<"here it is the negative cases:break"<<endl;
		      continue;
		    }
		}
	      else
		{
		  //case of positive deletions
		  assignment_params.d_ex_errs_right=0;
		  assignment_params.d_ex_errs_right_i.clear();
		  d_err_excess_pos_right.clear();

		}//end of error in excess

	      assignment_params.nerrorsd=assignment_params.d_errs+
		assignment_params.d_ex_errs_left+assignment_params.d_ex_errs_right;
	      if(assignment_params.nerrorsd>assignment_params.D_max_error)
		{
		  cout<<"BREKOUT:max d error:"<<assignment_params.D_max_error<<";nerrorsd:"<<assignment_params.nerrorsd<<endl;
		  continue;
		}
	      
	      assignment_params.log_highest_probability_GIVEN_current_Dr_deletions=-1000;
	      unsigned nerrors = assignment_params.nerrorsv+
		assignment_params.nerrorsd+ assignment_params.nerrorsj;
	      double log_perr;
	      if(nerrors==0)
		{
		  log_perr=0;
		}
	      else
		{
		  log_perr=nerrors*assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
		}

	    double log_max_pcutVDJ_loop_pJ=   assignment_params.log_PcutV+
	    assignment_params.log_PcutJ+assignment_params.log_max_pcutD_loop_d;

	    double test_proba=assignment_params.log_probabase+log_perr+
	      log_max_pcutVDJ_loop_pJ+log_max_pins_loop_dDr+log_p_max_nt_VD_DJ_na_loop_Dr;

	    if( (test_proba<assignment_params.log_probability_threshold_factor+assignment_params.log_highest_probability)||
		(test_proba<assignment_params.log_probability_hopeless_threshold_factor+assignment_params.log_highest_probability))
	      {
		cout<<"BREAKOUT due prob too low"<<endl;
		assignment_params.skips++;
		continue;
	      }
	    cout<<"ready to do the d palindrome........"<<endl;
	    //from this point on we are doing panlindrom
	    //define the loop iterator direction and steps
	    int npDl_start, npDl_end, npDr_start, npDr_end, np_step;
	    if(assignment_params.np_start_from_max)
	      {
		np_step=-1;
		npDl_start=npDl_max;
		npDl_end=0;
		npDr_start=npDr_max;
		npDr_end=0;
	      }
	    else
	      {
		np_step=1;
		npDl_start=0;
		npDl_end=npDl_max;
		npDr_start=0;
		npDr_end=npDr_max;
	      }
	    SequenceString insVD_nseq_min_looppDl;
	    SequenceString insDJ_nseq_min_looppDr;
	    double log_p_max_nt_VD_loop_pDl, log_p_max_nt_DJ_loop_pDr;
	    double log_p_max_nt_VD_DJ_na_loop_pDr;
	    //start the looping, over half-length of D left palindrome
	    for(assignment_params.npDl=npDl_start; ;assignment_params.npDl=+np_step)
	      {
		//check for the end of the loop
		if(assignment_params.np_start_from_max)//from big to small
		  {
		    if((signed)assignment_params.npDl<npDl_end)
		      break;
		  }
		else //from small to big
		  {
		    if((signed)assignment_params.npDl>npDl_end)
		      break;
		  }

		if(assignment_params.skips>assignment_params.max_skips)
		  {
		    _assigns.n_assignments = assignment_params.in;
		    _assigns.skips =  assignment_params.skips;
		    assignment_params.v_break_out=true;//===in this case, we want to jump all out. this is not necessary now. but just keep it here
		    return false;
		  }
		
		//start doing the loop 
		if(!zeroD)
		  {
		    if(assignment_params.V_end>=(signed)(_D.align_position_left[assignment_params.d][na]+assignment_params.ndDl1-assignment_params.npDl))
		      {
			//V with aplindrome overlaps .
			continue;
		      }
		    insVD_nseq_min_looppDl=
		      _seq.Sub(assignment_params.V_end+1, 
			       _D.align_position_left[assignment_params.d][na]+assignment_params.ndDl1-assignment_params.npDl-1);
		    if(insVD_nseq_min_looppDl.GetLength()>_model.max_insertions)
		      {
			continue;
		      }
		    //get sum stats for mono nts

		    log_p_max_nt_VD_loop_pDl=(insVD_nseq_min_looppDl.GetLetterCount('A')+insVD_nseq_min_looppDl.GetLetterCount('a'))*assignment_params.log_max_model_p_nt_VD(0);
		    log_p_max_nt_VD_loop_pDl+=(insVD_nseq_min_looppDl.GetLetterCount('C')+insVD_nseq_min_looppDl.GetLetterCount('c'))*assignment_params.log_max_model_p_nt_VD(1);
		    log_p_max_nt_VD_loop_pDl+=(insVD_nseq_min_looppDl.GetLetterCount('G')+insVD_nseq_min_looppDl.GetLetterCount('g'))*assignment_params.log_max_model_p_nt_VD(2);
		    log_p_max_nt_VD_loop_pDl+=(insVD_nseq_min_looppDl.GetLetterCount('T')+insVD_nseq_min_looppDl.GetLetterCount('t'))*assignment_params.log_max_model_p_nt_VD(3);
	 
		  }//zeroD case for the first loop 
		else //first zeroD case of else loop
		  {
		    log_p_max_nt_VD_loop_pDl=0;
		  }//end of first zeroD case of else loop
		
		//now loop over half-length of D right palindrome
		for(assignment_params.npDr=npDr_start; ; assignment_params.npDr+=np_step)
		  {
		    double log_max_pins_loop_pDr=0;
		    //check for the end of the loop
		    if(assignment_params.np_start_from_max)//from big to small
		      {
			if((signed)assignment_params.npDr<npDr_end)
			  break;
		      }
		    else //from small to big
		      {
			if((signed)assignment_params.npDr>npDr_end)
			  break;
		      }
		    //===>this following condition will not be used ever, need to remove later
		    //since this codition was checked above and there no update on skips in
		    //between till now.
		    if(assignment_params.skips>assignment_params.max_skips)
		      {
			_assigns.n_assignments = assignment_params.in;
			_assigns.skips =  assignment_params.skips;
			assignment_params.v_break_out=true;//===in this case, we want to jump all out.
			return false;
		      }

		    //if deletions is >0 and palindromes > 0 skip
		    if(!zeroD)
		      {
			if(assignment_params.J_start<=(signed)(_D.align_position_right[assignment_params.d][na]-assignment_params.ndDr1+assignment_params.npDr))
			  {
			    //J with palindrome overlaps D.
			    continue;
			  }
			
			//now doing the stats of mono-nt on DJ gap
			insDJ_nseq_min_looppDr=
			  _seq.Sub(_D.align_position_right[assignment_params.d][na]-assignment_params.ndDr1+assignment_params.npDr+1, assignment_params.J_start-1);
			if(insDJ_nseq_min_looppDr.GetLength()>_model.max_insertions)
			  {
			    continue;
			  }
			//get sum stats for mono nts
			
			log_p_max_nt_DJ_loop_pDr=(insDJ_nseq_min_looppDr.GetLetterCount('A')+insDJ_nseq_min_looppDr.GetLetterCount('a'))*assignment_params.log_max_model_p_nt_DJ(0);
			log_p_max_nt_DJ_loop_pDr+=(insDJ_nseq_min_looppDr.GetLetterCount('C')+insDJ_nseq_min_looppDr.GetLetterCount('c'))*assignment_params.log_max_model_p_nt_DJ(1);
			log_p_max_nt_DJ_loop_pDr+=(insDJ_nseq_min_looppDr.GetLetterCount('G')+insDJ_nseq_min_looppDr.GetLetterCount('g'))*assignment_params.log_max_model_p_nt_DJ(2);
			log_p_max_nt_DJ_loop_pDr+=(insDJ_nseq_min_looppDr.GetLetterCount('T')+insDJ_nseq_min_looppDr.GetLetterCount('t'))*assignment_params.log_max_model_p_nt_DJ(3);
			log_p_max_nt_VD_DJ_na_loop_pDr=log_p_max_nt_VD_loop_pDl+log_p_max_nt_DJ_loop_pDr;
		      }//npDr zeroD case 
		    else
		      {
			unsigned niVD_DJ_min=assignment_params.niVD_DJ_total;
			log_p_max_nt_VD_DJ_na_loop_pDr=niVD_DJ_min*assignment_params.log_max_model_p_nt;
		      }//npDr zeroD else case
		

		    //now start to figuring out the insertion on each side of D
		    //for this one, we need to consider the zeroD case, since it is more
		    //complicated. for the non-zero case, it is simple, because everything is 
		    //fixed. But zeroD case, we can allow one to vary, but the total is fixed
		    unsigned niVD_min, niVD_max;
		    unsigned niVD_nonzero=-1, niDJ_nonzero=-1;
		    
		    if(zeroD)
		      {//for zero D, we convolute over VD and DJ insertions
			if(assignment_params.niVD_DJ_total>2*_model.max_insertions)
			  {
			    continue;
			  }
			
			if(_do_smoothing)
			  {
			    niVD_min=0;
			    if(assignment_params.niVD_DJ_total>_model.max_insertions)
			      {
				niVD_min=assignment_params.niVD_DJ_total-_model.max_insertions;
			      }
			    niVD_max=_model.max_insertions;
			    if(assignment_params.niVD_DJ_total<_model.max_insertions)
			      {
				niVD_max=assignment_params.niVD_DJ_total;
			      }
			    
			    //check for sick cases
			    if(niVD_min>_model.max_insertions || 
			       ((unsigned)(assignment_params.niVD_DJ_total-niVD_max))>_model.max_insertions)
			      {
				continue;
			      }
			  }
			else //do smoothing else loop
			  {
			    //for strawman, we just split into two
			    niVD_min=assignment_params.niVD_DJ_total/2;
			    niVD_max=niVD_min;
			    //check for sick cases
			    if(niVD_min>_model.max_insertions || 
			       ((unsigned)(assignment_params.niVD_DJ_total-niVD_max))>_model.max_insertions)
			      {
				continue;
			      }
			  }//end of do smoothing else loop
			//doing stats
			log_max_pins_loop_pDr=-10000000;
			//go through a loop to get the biggest element in the combination
			//
			for(unsigned i=niVD_min;i<=niVD_max;i++)
			  {
			    if(log_max_pins_loop_pDr<log(_model.PinsVD(i)*_model.PinsDJ(assignment_params.niVD_DJ_total-i)))
			      {
				log_max_pins_loop_pDr=log(_model.PinsVD(i)*_model.PinsDJ(assignment_params.niVD_DJ_total-i));
				
			      }
			  }//for loop
			    
		      }
		    else //zero D else loop for insertion convolutions.
		      {
			niVD_nonzero=_D.align_position_left[assignment_params.d][na]+assignment_params.ndDl1-assignment_params.npDl-assignment_params.V_end-1;
			niDJ_nonzero=assignment_params.J_start -(_D.align_position_right[assignment_params.d][na]-assignment_params.ndDr1+assignment_params.npDr)-1;
			niVD_min=niVD_nonzero;
			niVD_max=niVD_nonzero;

			if(niVD_nonzero>_model.max_insertions || niDJ_nonzero> _model.max_insertions)
			  {
			    continue;
			  }
			
			log_max_pins_loop_pDr=log(_model.PinsDJ(niDJ_nonzero)*_model.PinsVD(niVD_nonzero));
						
		      }//end zero d else loop for insertion convolutions
		    
		    //start doing cut variables
		    assignment_params.ncutDl=assignment_params.ndDl-assignment_params.npDl;
		    assignment_params.ncutDr=assignment_params.ndDr-assignment_params.npDr;
		    //D deletions probability
		    double PcutDlDr=
		      _model.PcutDlcutDr_given_D(assignment_params.d,assignment_params.ncutDl-_model.min_D_cut,assignment_params.ncutDr-_model.min_D_cut);
		    
		    if(PcutDlDr==0)
		      {
			assignment_params.skips++;
			continue;
		      }
		    
		    assignment_params.log_PcutVDJ=
		      log(PcutDlDr)+assignment_params.log_PcutV
		      +assignment_params.log_PcutJ;
		    double log_perr=nerrors*assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
		    test_proba=assignment_params.log_probabase+
		      log_perr+assignment_params.log_PcutVDJ+
		      log_max_pins_loop_pDr+log_p_max_nt_VD_DJ_na_loop_pDr;

		    if((test_proba<assignment_params.log_probability_threshold_factor+assignment_params.log_highest_probability)||(test_proba<assignment_params.log_probability_hopeless_threshold_factor+assignment_params.log_highest_probability))
		      {
			assignment_params.skips++;
			continue;
		      }
		      
		    //loop over niVD for zero D case.
		    //for nonzero D, this is single instance loop
		    for(assignment_params.niVD=niVD_min;assignment_params.niVD<=niVD_max;assignment_params.niVD++)
		      {

			if(zeroD)
			  {
			    assignment_params.niDJ=assignment_params.niVD_DJ_total-assignment_params.niVD;
			  }
			else
			  {
			    assignment_params.niDJ=niDJ_nonzero;
			  }

			//check for other constraints
			//for instace, insertions can become negative 
			//with 'negative' deletions and palindromes together
			if((assignment_params.niDJ>_model.max_insertions)||
			   (assignment_params.niVD>_model.max_insertions))
			  {
			    continue;
			  }
			assignment_params.zeroD=zeroD;
			//====>>>>>>>from this point one, we generating all the stats
			//	     and ready to put them into assigns
			//by calling the function run_stats
			if(! run_stats_for_assignment(_model, _seq, _V, _D, _J,
						      _genV, _numV, _genD, _numD,
						      _genJ, _numJ, _no_error,
						      /*_ignore_deep_error,*/ _do_smoothing,
						      _force_all_alleles, _READ_LENGTH_CORRECTION,
						      assignment_params,/*output*/_assigns)
			   )
			  {
			    break;//because in matlab code it does this jump out too
			  }
			     
		      }//for loop of niVD for zeroD cases, nonzeroD case is single loop
		    if(assignment_params.v_break_out||assignment_params.j_break_out||assignment_params.d_break_out)
		      {
			break;
		      }
		  }//half-length of D right palindrome
	       	if(assignment_params.v_break_out||assignment_params.j_break_out||assignment_params.d_break_out)
		  {
		    break;
		  }
	      }//half-length of D left palindrome
	    
	    //here in ending the right deletions, we need to check for exit/breakingup cases
	    //==>unsigned d_errs_internal_max=_D.n_errors[assignment_params.d][na];
	    //here we define this as a local variable 
	    if(!zeroD && 
	       assignment_params.log_highest_probability_GIVEN_current_Dr_deletions>-900  &&
	       ( ( (assignment_params.ndDr1>2 && d_errs_internal_max==0)&&(assignment_params.log_highest_probability_GIVEN_current_Dr_deletions < (assignment_params.log_highest_probability+assignment_params.log_probability_threshold_factor))
		  )
		||
		 ( (assignment_params.ndDr1>0 && d_errs_internal_max==0)&&(assignment_params.log_highest_probability_GIVEN_current_Dr_deletions < (assignment_params.log_highest_probability+assignment_params.log_probability_hopeless_threshold_factor))
		  )
	       )
	      )
	      {
		break;
	      }
	    if(assignment_params.v_break_out||assignment_params.j_break_out||assignment_params.d_break_out)
	      {
		break;
	      }
	    
	    }//D right deletions, ndDr1
	  
	  //here in ending the right deletions, we need to check for exit/breakingup cases
	  
	  //here we define this as a local variable 
	  if(!zeroD && 
	     assignment_params.log_highest_probability_GIVEN_current_Dl_deletions>-900  &&
	     ( ( (assignment_params.ndDl1>2 && d_errs_internal_max==0)&&(assignment_params.log_highest_probability_GIVEN_current_Dl_deletions < (assignment_params.log_highest_probability+assignment_params.log_probability_threshold_factor))
		 )
	       ||
	       ( (assignment_params.ndDl1>0 && d_errs_internal_max==0)&&(assignment_params.log_highest_probability_GIVEN_current_Dl_deletions < (assignment_params.log_highest_probability+assignment_params.log_probability_hopeless_threshold_factor))
		 )
	       )
	     )
	    {
	      break;
	    }
	  if(assignment_params.v_break_out||assignment_params.j_break_out||assignment_params.d_break_out)
	    {
	      break;
	    }
	}//D left deletions, ndDl1
      if(assignment_params.v_break_out||assignment_params.d_break_out||assignment_params.j_break_out)
	{
	  break;
	}
    }//end of loop of D_alginments, out most.
     
  return true;
}//end of assign_D function

//start doing assignment with D seg on both sides
/*output: there is only one case under which we will get a false return
 *        it is in the very end of this function, we check for n_assignment_v/d/j
 *        so we will take care this on the outer loops for v/d/j_breakout
 *        Also, I think it is correct to set up to check for v/d/j_breakout at every loop
 *        because as I figured out, v/d/j_breakout is really jump_out_to_next_v/d/j_allele
 *        so in this sense it is jump to v/d/j loops.
 *     All other cases, a true is returned.
 */

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
 )
{
  cout<<"---------> start doing run stats:"<<endl;
  //first thing here is run stats for nt distribution
  SequenceString insVD_nseq=_seq.Sub(assignment_params.V_end+1, assignment_params.V_end+assignment_params.niVD);
  SequenceString insDJ_nseq=_seq.Sub(assignment_params.J_start-assignment_params.niDJ, assignment_params.J_start-1);

  SequenceString insVD_nseq_5prime_shift=_seq.Sub(assignment_params.V_end, assignment_params.V_end+assignment_params.niVD-1);
  SequenceString insDJ_nseq_3prime_shift=_seq.Sub(assignment_params.J_start-assignment_params.niDJ+1, assignment_params.J_start);

  //now doing mono nt distr first
  unsigned mononucleotideVD[]={0,0,0,0};
  unsigned mononucleotideDJ[]={0,0,0,0};

  mononucleotideVD[0]=insVD_nseq.GetLetterCount('A')+insVD_nseq.GetLetterCount('a');
  mononucleotideVD[1]=insVD_nseq.GetLetterCount('C')+insVD_nseq.GetLetterCount('c');
  mononucleotideVD[2]=insVD_nseq.GetLetterCount('G')+insVD_nseq.GetLetterCount('g');
  mononucleotideVD[3]=insVD_nseq.GetLetterCount('T')+insVD_nseq.GetLetterCount('t');
  
  mononucleotideDJ[0]=insDJ_nseq.GetLetterCount('A')+insDJ_nseq.GetLetterCount('a');
  mononucleotideDJ[1]=insDJ_nseq.GetLetterCount('C')+insDJ_nseq.GetLetterCount('c');
  mononucleotideDJ[2]=insDJ_nseq.GetLetterCount('G')+insDJ_nseq.GetLetterCount('g');
  mononucleotideDJ[3]=insDJ_nseq.GetLetterCount('T')+insDJ_nseq.GetLetterCount('t');
  
  //next the dinucleotide distribution
  unsigned matrix_dim[]={4,4,1,1}; //dimension of 4, but not necessarily all used, could be
  //only first say 2 to be used
  Matrix<unsigned> nucleotideVD(2, matrix_dim, (unsigned)0);
  Matrix<unsigned> nucleotideDJ(2, matrix_dim, (unsigned)0);
  
  unsigned cor_x, cor_y;
  //for VD
  for(unsigned nt_VD=0; nt_VD<assignment_params.niVD;nt_VD++)
    {
      switch (insVD_nseq.GetSequence().at(nt_VD))
	{
	case 'A':
	case 'a':
	  cor_x=0;
	  break;
	case 'C':
	case 'c':
	  cor_x=1;
	  break;
	 
	case 'G':
	case 'g':
	  cor_x=2;
	  break;
	case 'T':
	case 't':
	  cor_x=3;
	  break;
	default:
	  cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");
	}
      
      switch (insVD_nseq_5prime_shift.GetSequence().at(nt_VD))
	{
	case 'A':
	case 'a':
	  cor_y=0;
	  break;
	case 'C':
	case 'c':
	  cor_y=1;
	  break;
	 
	case 'G':
	case 'g':
	  cor_y=2;
	  break;
	case 'T':
	case 't':
	  cor_y=3;
	  break;
	default:
	  cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");
	}
      
      nucleotideVD(cor_x, cor_y)=nucleotideVD(cor_x, cor_y)+1;
    }//end of vd di nucleotide for loop
  
    //for DJ
  for(unsigned nt_DJ=0; nt_DJ<assignment_params.niDJ;nt_DJ++)
    {
      switch (insDJ_nseq.GetSequence().at(nt_DJ))
	{
	case 'A':
	case 'a':
	  cor_x=0;
	  break;
	case 'C':
	case 'c':
	  cor_x=1;
	  break;
	 
	case 'G':
	case 'g':
	  cor_x=2;
	  break;
	case 'T':
	case 't':
	  cor_x=3;
	  break;
	default:
	  cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");
	}
      
      switch (insDJ_nseq_3prime_shift.GetSequence().at(nt_DJ))
	{
	case 'A':
	case 'a':
	  cor_y=0;
	  break;
	case 'C':
	case 'c':
	  cor_y=1;
	  break;
	 
	case 'G':
	case 'g':
	  cor_y=2;
	  break;
	case 'T':
	case 't':
	  cor_y=3;
	  break;
	default:
	  cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");
	}
      nucleotideDJ(cor_x, cor_y)=nucleotideDJ(cor_x, cor_y)+1;
    }//end of vd di nucleotide for loop

  //% Compute final probability of assignment by multiplying base probability by factors for
  //% palindromes, deletions and the number of errors.
  //% multiply probability by factor for insertions of given lengths and sequence specific factors for
  //% nucleotide bias
  double log_pins=log(_model.PinsVD(assignment_params.niVD))*log(_model.PinsDJ(assignment_params.niDJ));
	
  double log_pntbias_VD=
    matrix_multiply_1D(assignment_params.log_RnucleotideVD_per_nucleotideVD_5prime,nucleotideVD.m2vec());
  double log_pntbias_DJ=
    matrix_multiply_1D(assignment_params.log_RnucleotideDJ_per_nucleotideDJ_3prime,nucleotideDJ.m2vec());

  //now finally, this is accurate!!!
  double log_proba=assignment_params.log_probabase+
    (assignment_params.nerrorsv+assignment_params.nerrorsj+assignment_params.nerrorsd)*assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3+
    assignment_params.log_PcutVDJ+log_pins+log_pntbias_VD+log_pntbias_DJ;

  //update all highest probability so far values, if necessary.
  if(log_proba>assignment_params.log_highest_probability_GIVEN_current_V_allele)
    {
      assignment_params.log_highest_probability_GIVEN_current_V_allele=log_proba;
    }
  if(log_proba>assignment_params.log_highest_probability_GIVEN_current_J_allele)
    {
      assignment_params.log_highest_probability_GIVEN_current_J_allele=log_proba;
    }
  if(log_proba>assignment_params.log_highest_probability_GIVEN_current_D_allele)
    {
      assignment_params.log_highest_probability_GIVEN_current_D_allele=log_proba;
    }
  if(log_proba>assignment_params.log_highest_probability_GIVEN_current_V_deletions)
    {
      assignment_params.log_highest_probability_GIVEN_current_V_deletions=log_proba;
    }
  if(log_proba>assignment_params.log_highest_probability_GIVEN_current_J_deletions)
    {
      assignment_params.log_highest_probability_GIVEN_current_J_deletions=log_proba;
    }
  if(log_proba>assignment_params.log_highest_probability_GIVEN_current_Dl_deletions)
    {
      assignment_params.log_highest_probability_GIVEN_current_Dl_deletions=log_proba;
    }
  if(log_proba>assignment_params.log_highest_probability_GIVEN_current_Dr_deletions)
    {
      assignment_params.log_highest_probability_GIVEN_current_Dr_deletions=log_proba;
    }

//the following two if statements are doing extra checks
  //==========> but probably unnecessary, since we have check similar conditions in above.
  //we have them here for now. might remove them in the future, 6/21/1015.
  if((signed)assignment_params.D_align_length<0)
    {
      assignment_params.D_align_length=0;
    }
  
  //update the best align length
  if(assignment_params.D_align_length> assignment_params.best_D_align_length)
    {
      assignment_params.best_D_align_length=assignment_params.D_align_length;
    }

  //% if D length isvery short (<3) I let it count low prob events because the total prob of D short
  //% is spread out over many choices of D deletions. So if I skip them,
  //%if (~do_zeroD && (log_proba < log_probability_threshold_factor + log_highest_probability)) || (log_proba < log_probability_hopeless_threshold_factor + log_highest_probability)
  if ((log_proba < assignment_params.log_probability_threshold_factor + assignment_params.log_highest_probability) || (log_proba < assignment_params.log_probability_hopeless_threshold_factor + assignment_params.log_highest_probability))
    {
      assignment_params.skips = assignment_params.skips + 1;
      return true;//continue; here, we skip the function. it is more like to skip the for loop
    }
  
  //dinucleotide???
  Matrix<unsigned> nucleotideVD_5prime(2, matrix_dim, (unsigned)0);
  Matrix<unsigned> nucleotideDJ_3prime(2, matrix_dim, (unsigned)0);
  
  Matrix<unsigned> sum_row=sum(nucleotideVD,0); //sum along the first dimension. get second dimension reserved
  nucleotideVD_5prime.SetSubMatrix(0,sum_row);
  nucleotideVD_5prime.SetSubMatrix(1,sum_row);
  nucleotideVD_5prime.SetSubMatrix(2,sum_row);
  nucleotideVD_5prime.SetSubMatrix(3,sum_row);
  
  Matrix<unsigned> sum_row2=sum(nucleotideDJ,0);
  nucleotideDJ_3prime.SetSubMatrix(0,sum_row2);
  nucleotideDJ_3prime.SetSubMatrix(1,sum_row2);
  nucleotideDJ_3prime.SetSubMatrix(2,sum_row2);
  nucleotideDJ_3prime.SetSubMatrix(3,sum_row2);
  
  sum_row.clear();
  sum_row2.clear();
  
  //edge dinucleotide distribution
  matrix_dim[0]=4;matrix_dim[1]=4;
  Matrix<unsigned> VD_left_edge_dinucleotide(2,matrix_dim,(unsigned)0);
  Matrix<unsigned> VD_right_edge_dinucleotide(2,matrix_dim,(unsigned)0);
  
  Matrix<unsigned> DJ_left_edge_dinucleotide(2,matrix_dim, (unsigned)0);
  Matrix<unsigned> DJ_right_edge_dinucleotide(2,matrix_dim, (unsigned)0);
  
  if(assignment_params.niVD >0)
    {
      switch (_seq.GetSequence().at(assignment_params.V_end))
	{
	case 'A':
	case 'a':
	  cor_x=0;
	  break;
	case 'C':
	case 'c':
	  cor_x=1;
	  break;
	case 'G':
	case 'g':
	  cor_x=2;
	  break;
	case 'T':
	case 't':
	  cor_x=3;
	  break;
	default:
	  cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	}
      switch (_seq.GetSequence().at(assignment_params.V_end+1))
	{
	case 'A':
	case 'a':
	  cor_y=0;
	  break;
	case 'C':
	case 'c':
	  cor_y=1;
	  break;
	case 'G':
	case 'g':
	  cor_y=2;
	  break;
	case 'T':
	case 't':
	  cor_y=3;
	  break;
	default:
	  cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	}
      VD_left_edge_dinucleotide(cor_x, cor_y)++;
      //now the right edge
      if(!assignment_params.zeroD)
	{
	  unsigned D_left_start=_D.align_position_left[assignment_params.d][assignment_params.na]+assignment_params.ndDl1-assignment_params.npDl;

	  switch (_seq.GetSequence().at(D_left_start-1))
	    {
	    case 'A':
	    case 'a':
	      cor_x=0;
	      break;
	    case 'C':
	    case 'c':
	      cor_x=1;
	      break;
	    case 'G':
	    case 'g':
	      cor_x=2;
	      break;
	    case 'T':
	    case 't':
	      cor_x=3;
	      break;
	    default:
	      cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	    }
	  switch (_seq.GetSequence().at(D_left_start))
	    {
	    case 'A':
	    case 'a':
	      cor_y=0;
	      break;
	    case 'C':
	    case 'c':
	      cor_y=1;
	      break;
	    case 'G':
	    case 'g':
	      cor_y=2;
	      break;
	    case 'T':
	    case 't':
	      cor_y=3;
	      break;
	    default:
	      cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	    }
	  VD_right_edge_dinucleotide(cor_x, cor_y)++;
	}//if loop, zeroD case for right edge VD

    }//niVD not zero if loop

  //now doing dj part, right edge
  if(assignment_params.niDJ>0)
    {
      
      switch (_seq.GetSequence().at(assignment_params.J_start-1))
	{
	case 'A':
	case 'a':
	  cor_x=0;
	  break;
	case 'C':
	case 'c':
	  cor_x=1;
	  break;
	case 'G':
	case 'g':
	  cor_x=2;
	  break;
	case 'T':
	case 't':
	  cor_x=3;
	  break;
	default:
	  cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	}
      switch (_seq.GetSequence().at(assignment_params.J_start))
	{
	case 'A':
	case 'a':
	  cor_y=0;
	  break;
	case 'C':
	case 'c':
	  cor_y=1;
	  break;
	case 'G':
	case 'g':
	  cor_y=2;
	  break;
	case 'T':
	case 't':
	  cor_y=3;
	  break;
	default:
	  cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	  throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	}
      DJ_right_edge_dinucleotide(cor_x, cor_y)++;
      //now the dj left edge
      if(!assignment_params.zeroD)
	{
	  unsigned D_right_end=_D.align_position_right[assignment_params.d][assignment_params.na]-assignment_params.ndDr1+assignment_params.npDr;

	  switch (_seq.GetSequence().at(D_right_end))
	    {
	    case 'A':
	    case 'a':
	      cor_x=0;
	      break;
	    case 'C':
	    case 'c':
	      cor_x=1;
	      break;
	    case 'G':
	    case 'g':
	      cor_x=2;
	      break;
	    case 'T':
	    case 't':
	      cor_x=3;
	      break;
	    default:
	      cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	    }
	  switch (_seq.GetSequence().at(D_right_end))
	    {
	    case 'A':
	    case 'a':
	      cor_y=0;
	      break;
	    case 'C':
	    case 'c':
	      cor_y=1;
	      break;
	    case 'G':
	    case 'g':
	      cor_y=2;
	      break;
	    case 'T':
	    case 't':
	      cor_y=3;
	      break;
	    default:
	      cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	    }
	  DJ_left_edge_dinucleotide(cor_x, cor_y)++;
	}//if loop, zeroD case for right edge VD
      
    }//if loop, niDJ not zero

  //now start doing the tri-nucleotide
  matrix_dim[2]=4;
  Matrix<unsigned> trinucleotideVD(3,matrix_dim, (unsigned)0);
  Matrix<unsigned> trinucleotideDJ(3,matrix_dim, (unsigned)0);

  unsigned cor_z;

  for(signed nt_VD=0;nt_VD<((signed)assignment_params.niVD-2);nt_VD++)
    {
      switch (insVD_nseq.GetSequence().at(nt_VD))
	    {
	    case 'A':
	    case 'a':
	      cor_x=0;
	      break;
	    case 'C':
	    case 'c':
	      cor_x=1;
	      break;
	    case 'G':
	    case 'g':
	      cor_x=2;
	      break;
	    case 'T':
	    case 't':
	      cor_x=3;
	      break;
	    default:
	      cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	    }
	  switch (insVD_nseq.GetSequence().at(nt_VD+1))
	    {
	    case 'A':
	    case 'a':
	      cor_y=0;
	      break;
	    case 'C':
	    case 'c':
	      cor_y=1;
	      break;
	    case 'G':
	    case 'g':
	      cor_y=2;
	      break;
	    case 'T':
	    case 't':
	      cor_y=3;
	      break;
	    default:
	      cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	    }
	  switch (insVD_nseq.GetSequence().at(nt_VD+2))
	    {
	    case 'A':
	    case 'a':
	      cor_z=0;
	      break;
	    case 'C':
	    case 'c':
	      cor_z=1;
	      break;
	    case 'G':
	    case 'g':
	      cor_z=2;
	      break;
	    case 'T':
	    case 't':
	      cor_z=3;
	      break;
	    default:
	      cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	    }
	  trinucleotideVD(cor_x, cor_y, cor_z)++;
    }//for loop for nt_VD trinucleotide;
  
  for(signed nt_DJ=0;nt_DJ<((signed)assignment_params.niDJ-2);nt_DJ++)
    {
      switch (insDJ_nseq.GetSequence().at(nt_DJ))
	    {
	    case 'A':
	    case 'a':
	      cor_x=0;
	      break;
	    case 'C':
	    case 'c':
	      cor_x=1;
	      break;
	    case 'G':
	    case 'g':
	      cor_x=2;
	      break;
	    case 'T':
	    case 't':
	      cor_x=3;
	      break;
	    default:
	      cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	    }
	  switch (insDJ_nseq.GetSequence().at(nt_DJ+1))
	    {
	    case 'A':
	    case 'a':
	      cor_y=0;
	      break;
	    case 'C':
	    case 'c':
	      cor_y=1;
	      break;
	    case 'G':
	    case 'g':
	      cor_y=2;
	      break;
	    case 'T':
	    case 't':
	      cor_y=3;
	      break;
	    default:
	      cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	    }
	  switch (insDJ_nseq.GetSequence().at(nt_DJ+2))
	    {
	    case 'A':
	    case 'a':
	      cor_z=0;
	      break;
	    case 'C':
	    case 'c':
	      cor_z=1;
	      break;
	    case 'G':
	    case 'g':
	      cor_z=2;
	      break;
	    case 'T':
	    case 't':
	      cor_z=3;
	      break;
	    default:
	      cout<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      cerr<<"unknown character in the sequence. please check in 'run stat()function in vdj model assignment";
	      throw runtime_error ("unknown character in the sequence. please check in 'run stat()function in vdj model assignment");	  
	    }
	  trinucleotideDJ(cor_x, cor_y, cor_z)++;
    }//for loop for nt_DJ trinucleotide;
  
  //start doing errors, this need to be reversed to starting from the end of J chain
  unsigned* error_vs_position=new unsigned[_model.model_params.maximum_read_length];
  memset(error_vs_position, 0, sizeof(unsigned)*_model.model_params.maximum_read_length);
  //to store positions with errors for this assignment
  
  //total nucleotides that are assigned to either V, D or J in this assignment
  //
  //unsigned genic_length;
  assignment_params.genic_length=_V.align_length[assignment_params.v]+_J.align_length[assignment_params.j]
    +assignment_params.npV+assignment_params.npJ-assignment_params.ndV1 - 
    assignment_params.ndJ1;
  if(!assignment_params.zeroD)
    {
      assignment_params.genic_length=assignment_params.genic_length
	+_D.align_length[assignment_params.d][assignment_params.na]
	+assignment_params.npDl+assignment_params.npDr-assignment_params.ndDl1-
	assignment_params.ndDr1;	
    }
  //in case of zeroD, the genic_length is good by now.

//the following two if statements are doing extra checks
  //==========> but probably unnecessary, since we have check similar conditions in above.
  //we have them here for now. might remove them in the future, 6/21/1015.  
  if((signed)assignment_params.genic_length<0)
    {
      //assignment_params.zeroD=true;
      //jump out this case;
      return true; //but don't affect the next case.
    }
  
  //here as in the matlab code, we skip some code for very short genic length. check matlab

  cout<<"check piont at L2543"<<endl;

  //store all positions that have been assigned to either v, d, or j.
  //this is used to estimate error rate vs psotion by dividing errors in a position by coverage
  unsigned v=assignment_params.v;
  unsigned j=assignment_params.j;
  unsigned d=assignment_params.d;
  unsigned na=assignment_params.na;

  unsigned* coverage =new unsigned[_model.model_params.maximum_read_length];
  memset(coverage, 0, sizeof(unsigned)*_model.model_params.maximum_read_length);
  unsigned V_start=_V.align_position[v][0];
  unsigned V_end=assignment_params.V_end;
  unsigned J_start=assignment_params.J_start;
  unsigned J_end=_J.align_position[j][0]+(_genJ[assignment_params.j_a].Get_Sequence().size()-_J.align_position[j][1])-1;
  
  unsigned starting_position=J_end;
  V_start=starting_position-V_start;
  V_end=starting_position-V_end;

  J_start=starting_position-J_start;
  J_end=starting_position-J_end;

  unsigned D_start, D_end;

  //set v segments to 1, reverse order
  for(unsigned i=V_end;i>=V_start;i++)
    {
      coverage[i]=1;
    }
  //set j segments to 1, reverse order
  for(unsigned i=J_end;i>=J_start;i++)
    {
      coverage[i]=1;
    }


  if(!assignment_params.zeroD)
    {
      D_start=_D.align_position_left[d][na]+
	assignment_params.ndDl1-assignment_params.npDl;
      D_end=_D.align_position_right[d][na]-
	assignment_params.ndDr1+assignment_params.npDr;
      
      D_start=starting_position-D_start;
      D_end=starting_position-D_end;

      for(unsigned i=D_end;i>=D_start;i++)
	{
	  coverage[i]=1;
	}
    }

  //store potential pal that is consistent with assignment.
  if(assignment_params.npV_potential_max >(assignment_params.niVD+assignment_params.npV))
    {
      assignment_params.npV_potential_max=assignment_params.niVD+assignment_params.npV;
    }
  if(assignment_params.npJ_potential_max >(assignment_params.niDJ+assignment_params.npJ))
    {
      assignment_params.npJ_potential_max=assignment_params.niDJ+assignment_params.npJ;
    }

  if(!assignment_params.zeroD)
    {
      if(assignment_params.npDl_potential_max>assignment_params.niVD+assignment_params.npDl)
	{
	  assignment_params.npDl_potential_max=assignment_params.niVD+assignment_params.npDl;
	}

      if(assignment_params.npDr_potential_max>assignment_params.niDJ+assignment_params.npDr)
	{
	  assignment_params.npDr_potential_max=assignment_params.niDJ+assignment_params.npDr;
	}
    }
  else
    {
      assignment_params.npDl_potential_max=0;
      assignment_params.npDr_potential_max=0;
    }
  
  //===========important node, this following incremental code are in the Matlab code.
  // we moved them to the end of the updating section for c++ code. be careful.
  //increment valid assignment number
  /*assignment_params.in++;
  assignment_params.n_assignments_v_gene++;
  assignment_params.n_assignments_d_gene++;
  assignment_params.n_assignments_j_gene++;
  */
  //doing errors now
  unsigned* v_err_pos_rel=new unsigned[assignment_params.max_V_depth];
  memset(v_err_pos_rel, 0, sizeof(unsigned)*assignment_params.max_V_depth);
  unsigned* j_err_pos_rel=new unsigned[assignment_params.max_J_depth];
  memset(j_err_pos_rel, 0, sizeof(unsigned)*assignment_params.max_J_depth);
  unsigned pos;
  //
  if(!_no_error)
    {
      matrix_dim[0]=3;
      Matrix<unsigned> v_err_excess_pos(1, matrix_dim,_V.excess_error_positions[v]); //% error positions in extended V alignment for 'negative' deleti
      matrix_dim[0]=_V.n_errors[v];
      Matrix<unsigned> v_err_pos(1, matrix_dim, _V.error_positions[v]);
      
      Matrix<unsigned> v_err_excess_pos_ok;
      Matrix<unsigned> v_err_pos_ok;
      //unsigned pos;
      if(assignment_params.nerrorsv>0)
	{
	  if(assignment_params.ndV1<0)//negative V deletions case
	    {
	      if(assignment_params.v_ex_errs>0)//that means we have negative errors, but can ONLY
		//be 1.
		{
		  v_err_excess_pos_ok=v_err_excess_pos.GetElements(assignment_params.v_ex_errs_i);
		  //so v_err_excess_pos_ok contains positions that we want to set to be in the error_vs_position array, this is contribution from negative case, but considered as a error in the sequences
		  //use for loop to set it
		  for(unsigned i=0;i<v_err_excess_pos_ok.size(0);i++)
		    {
		      pos=v_err_excess_pos_ok(i);
		      pos=starting_position-pos; //revers it from the j end;
		      error_vs_position[pos]=1;
		    }
		}//end of non zero negative v deletions case
	      v_err_pos_ok=v_err_pos;
	      for(unsigned i=0;i<v_err_pos_ok.size(0);i++)
		{
		  pos=v_err_pos_ok(i);
		  pos=starting_position-pos; //revers it from the j end;
		  error_vs_position[pos]=1;
		}
	    }
	  else //positive V deletions case
	    {
	      Matrix<bool> v_err_pos_ok_i= v_err_pos < assignment_params.V_end-assignment_params.npV+1; //this is equivalent to the matlab code, the idea is to remove npV in the end and only count the aligned errors
	      v_err_pos_ok=v_err_pos.GetElements(v_err_pos_ok_i);
	      for(unsigned i=0;i<v_err_pos_ok.size(0);i++)
		{
		  pos=v_err_pos_ok(i);
		  pos=starting_position-pos; //revers it from the j end;
		  error_vs_position[pos]=1;
		}	      
	    }//end of else loop that is positive V deletion case

	  //get relative err pos for v
	  for(unsigned i=0;i<_V.align_length[v];i++)
	    {
	      pos=i+assignment_params.ndV;
	      pos=starting_position-pos;
	      v_err_pos_rel[i+assignment_params.ndV]=error_vs_position[pos];
	    }
	}//end of nerrorsv if case

      //set positions of errors in J to accounting for deletions
      matrix_dim[0]=3;
      Matrix<unsigned> j_err_excess_pos(1, matrix_dim,_J.excess_error_positions[j]); //% error positions in extended V alignment for 'negative' deleti
      matrix_dim[0]=_J.n_errors[j];
      Matrix<unsigned> j_err_pos(1, matrix_dim, _J.error_positions[j]);
      
      Matrix<unsigned> j_err_excess_pos_ok;
      Matrix<unsigned> j_err_pos_ok;
      if(assignment_params.nerrorsj>0)
	{
	  if(assignment_params.ndJ1<0)//negative V deletions case
	    {
	      if(assignment_params.j_ex_errs>0)//that means we have negative errors, but can ONLY
		//be 1.
		{
		  j_err_excess_pos_ok=j_err_excess_pos.GetElements(assignment_params.j_ex_errs_i);
		  //so j_err_excess_pos_ok contains positions that we want to set to be in the error_vs_position array, this is contribution from negative case, but considered as a error in the sequences
		  //use for loop to set it
		  for(unsigned i=0;i<j_err_excess_pos_ok.size(0);i++)
		    {
		      pos=j_err_excess_pos_ok(i);
		      pos=starting_position-pos; //revers it from the j end;
		      error_vs_position[pos]=1;
		    }
		}//end of non zero negative v deletions case
	      j_err_pos_ok=j_err_pos;
	      for(unsigned i=0;i<j_err_pos_ok.size(0);i++)
		{
		  pos=j_err_pos_ok(i);
		  pos=starting_position-pos; //revers it from the j end;
		  error_vs_position[pos]=1;
		}
	    }
	  else //positive V deletions case
	    {
	      Matrix<bool> j_err_pos_ok_i= j_err_pos > assignment_params.J_start+assignment_params.npJ-1; //this is equivalent to the matlab code, the idea is to remove npV in the end and only count the aligned errors
	      j_err_pos_ok=j_err_pos.GetElements(j_err_pos_ok_i);
	      for(unsigned i=0;i<j_err_pos_ok.size(0);i++)
		{
		  pos=j_err_pos_ok(i);
		  pos=starting_position-pos; //revers it from the j end;
		  error_vs_position[pos]=1;
		}	      
	    }//end of else loop that is positive V deletion case

	  //get relative err pos for v
	  for(unsigned i=0;i<_J.align_length[j];i++)
	    {
	      pos=i+assignment_params.ndJ;
	      pos=starting_position-pos;
	      v_err_pos_rel[i+assignment_params.ndJ]=error_vs_position[pos];
	    }
	  
	}//end of nerrorsj if case

    }//end of if no error loop

  //now set up positions of errors in D to 1. all of these
  //come from 'negative' deleteions since D alignments do not allow error
  //the above is from Matlab, and it is true even in Matlab.
  //so be careful
  if(!assignment_params.zeroD&& (assignment_params.nerrorsd>0))
    {
      matrix_dim[0]=3;
      Matrix<unsigned> d_err_excess_pos_left(1, matrix_dim, _D.excess_error_positions_left[d][na]);
      Matrix<unsigned> d_err_excess_pos_right(1, matrix_dim, _D.excess_error_positions_right[d][na]);

      Matrix<unsigned> d_err_excess_pos_left_ok;
      Matrix<unsigned> d_err_excess_pos_right_ok;
      //for left
      if(assignment_params.ndDl1<0&&assignment_params.d_ex_errs_left>0)
	{
	  d_err_excess_pos_left_ok=d_err_excess_pos_left.GetElements(assignment_params.d_ex_errs_left_i);
	  for(unsigned i=0;i<d_err_excess_pos_left_ok.size(0);i++)
	    {
	      pos=d_err_excess_pos_left_ok(i);
	      pos=starting_position-pos;
	      error_vs_position[pos]=1;
	    }
	}

      //for right
      if(assignment_params.ndDr1<0&&assignment_params.d_ex_errs_right>0)
	{
	  d_err_excess_pos_right_ok=d_err_excess_pos_right.GetElements(assignment_params.d_ex_errs_right_i);
	  for(unsigned i=0;i<d_err_excess_pos_right_ok.size(0);i++)
	    {
	      pos=d_err_excess_pos_right_ok(i);
	      pos=starting_position-pos;
	      error_vs_position[pos]=1;
	    }
	}
      
      //for middle one according to 
      if(assignment_params.d_errs>0)
	{
	  unsigned temp_dim_size[]={_D.n_errors[d][na]};
	  Matrix<unsigned> d_err_pos(1, temp_dim_size, _D.error_positions[d][na]);
      
	  Matrix<unsigned> d_err_pos_ok=d_err_pos.GetElements(assignment_params.d_errs_i);
	  for(unsigned i=0;i<d_err_pos_ok.size(0);i++)
	    {
	      pos=d_err_pos_ok(i);
	      pos=starting_position-pos;
	      error_vs_position[pos]=1;
	    }
	}
    }//D error set up
  //finally we are done with assignments, next put things into assigns variable.

  cout<<"checkpoint L2805"<<endl;

  //=========================
  //% Store everything in assigns.
  //% Each variable in counter that begins with an 'nP' or 'nM' must be
  //% found in assigns without the 'nP' or 'nM'.
  //count<<"\t start updating >>>>>>>>>>>>"<<endl;
  //%%% The following are model variables
  unsigned in=assignment_params.in;                                                  
  _assigns.V(in) = assignment_params.v_g; //in-1, since we have incremented the 
  _assigns.DJ(in,0)=assignment_params.d_g;
  _assigns.DJ(in,1)=assignment_params.j_g; //% gene choices

  _assigns.Vallele_given_gene(in,0) =assignment_params.v_g;_assigns.Vallele_given_gene(in,1)= _genV[assignment_params.v_a].Get_Allele(); //% allele choice given gene
  _assigns.Dallele_given_gene(in,0) =assignment_params.d_g;_assigns.Dallele_given_gene(in,1)= _genD[assignment_params.d_a].Get_Allele();

  _assigns.cutV_given_V(in,0) = assignment_params.v_g;_assigns.cutV_given_V(in,1)= (assignment_params.ndV-assignment_params.npV)- _model.min_V_cut;
  _assigns.cutJ_given_J(in,0) = assignment_params.j_g; _assigns.cutJ_given_J(in,1) =(assignment_params.ndJ-assignment_params.npJ) - _model.min_J_cut ;
  _assigns.insVD(in)=assignment_params.niVD;// % insertions, niVD in matlab
  _assigns.insDJ(in)=assignment_params.niDJ;// niDJ in matlab

  _assigns.nucleotideVD.SetSubMatrix(in, nucleotideVD);
  _assigns.nucleotideVD_5prime.SetSubMatrix(in, nucleotideVD_5prime);
  _assigns.nucleotideDJ.SetSubMatrix(in,nucleotideDJ);
  _assigns.nucleotideDJ_3prime.SetSubMatrix(in, nucleotideDJ_3prime);

  cout<<"check point L2831"<<endl;

  _assigns.error(in) = assignment_params.nerrorsv+
    assignment_params.nerrorsd+ assignment_params.nerrorsj;//nerrors;//assignment_params.nerrorsv+assignment_params.nerrorsj+assignment_params.nerrorsd; //% numerator for error rate estimate
  unsigned nerrors=_assigns.error(in);
  _assigns.sequenced_nucleotide(in) = assignment_params.genic_length; //% denominator for error rate estimate
  //end of model variables........
  
  //% For tracking                                                  
  //% numbers of A,C,G and T in insertions: numerator for probabilities of A,C,G,T insertions
  cout<<"dimension of mono :"<<_assigns.mononucleotideVD.size().toString()<<endl;
  _assigns.mononucleotideVD(in,0) = mononucleotideVD[0];
  _assigns.mononucleotideVD(in,1) = mononucleotideVD[1];
  _assigns.mononucleotideVD(in,2) = mononucleotideVD[2];
  _assigns.mononucleotideVD(in,3) = mononucleotideVD[3];

  _assigns.insertionVD(in,0) =  assignment_params.niVD;
  _assigns.insertionVD(in,1) =  assignment_params.niVD;
  _assigns.insertionVD(in,2) =  assignment_params.niVD;
  _assigns.insertionVD(in,3) =  assignment_params.niVD;
  //, niVD, niVD, niVD]; //% denominator for probabilities of A,C,G,T insertions.
  _assigns.mononucleotideDJ(in,0) = mononucleotideDJ[0];
  _assigns.mononucleotideDJ(in,1) = mononucleotideDJ[1];
  _assigns.mononucleotideDJ(in,2) = mononucleotideDJ[2];
  _assigns.mononucleotideDJ(in,3) = mononucleotideDJ[3];

  _assigns.insertionDJ(in,0) = assignment_params.niDJ; 
  _assigns.insertionDJ(in,1) = assignment_params.niDJ;
  _assigns.insertionDJ(in,2) = assignment_params.niDJ;
  _assigns.insertionDJ(in,3) = assignment_params.niDJ; 
  //, niDJ, niDJ, niDJ]; //% denominator for probabilities of A,C,G,T insertions.
  //cout<<"check point, L2861"<<endl;
  //cout<<"size of 1:"<<_assigns.VD_left_edge_dinucleotide.dim()<<"dimension:"<<_assigns.VD_left_edge_dinucleotide.size().toString()<<endl;
  //cout<<"size of 2:"<<VD_left_edge_dinucleotide.dim()<<"dimension:"<<VD_left_edge_dinucleotide.size().toString()<<endl;
  _assigns.VD_left_edge_dinucleotide.SetSubMatrix(in, VD_left_edge_dinucleotide);
  _assigns.VD_right_edge_dinucleotide.SetSubMatrix(in, VD_right_edge_dinucleotide);
                                                    
  _assigns.DJ_left_edge_dinucleotide.SetSubMatrix(in, DJ_left_edge_dinucleotide);
  _assigns.DJ_right_edge_dinucleotide.SetSubMatrix(in, DJ_right_edge_dinucleotide);
                                                    
  _assigns.trinucleotideVD.SetSubMatrix(in,trinucleotideVD);
  _assigns.trinucleotideDJ.SetSubMatrix(in,trinucleotideDJ);

  cout<<"check point, L2871"<<endl;
  //                                             
  _assigns.VV_err_pos.SetSubMatrix(in,assignment_params.v_g, assignment_params.max_V_depth, v_err_pos_rel);
  _assigns.VV_align_length(in,0) = assignment_params.v_g;
  _assigns.VV_align_length(in,1) =assignment_params.V_align_length; //no need to be 1+V_align_length
                                                    
  _assigns.JJ_err_pos.SetSubMatrix(in,assignment_params.j_g, assignment_params.max_J_depth, j_err_pos_rel);
  _assigns.JJ_align_length(in,0) = assignment_params.j_g;
  _assigns.JJ_align_length(in,1)=assignment_params.J_align_length; //no need to be 1+J_align_length
  cout<<"check point, L2883"<<endl;                                                  
  _assigns.pVmax_delV_V(in,0) =assignment_params.v_g;
  _assigns.pVmax_delV_V(in,1) =assignment_params.ndV;
  _assigns.pVmax_delV_V(in,2) =assignment_params.npV_potential_max;
  //[1 + npV_potential_max, 1+ ndV, v_g];
  _assigns.pJmax_delJ_J(in,0) =assignment_params.j_g;
  _assigns.pJmax_delJ_J(in,1) =assignment_params.ndJ;
  _assigns.pJmax_delJ_J(in,2) =assignment_params.npJ_potential_max;
  //[1 + npJ_potential_max, 1+ ndJ, j_g];
                                                        
  _assigns.pDlmax_delDl_D(in,0) = assignment_params.d_g;
  _assigns.pDlmax_delDl_D(in,1) = assignment_params.ndDl;
  _assigns.pDlmax_delDl_D(in,2) = assignment_params.npDl_potential_max;
  //[1 + npDl_potential_max, 1+ ndDl, d_g];
  _assigns.pDrmax_delDr_D(in,0) = assignment_params.d_g;
  _assigns.pDrmax_delDr_D(in,1) = assignment_params.ndDr;
  _assigns.pDrmax_delDr_D(in,2) = assignment_params.npDr_potential_max;
  //[1 + npDr_potential_max, 1+ ndDr, d_g];
                                                    
  _assigns.zeroD(in) = assignment_params.zeroD;
  _assigns.V_align_length(in) = assignment_params.V_align_length;
  _assigns.D_align_length(in) = assignment_params.D_align_length;
  _assigns.J_align_length(in) = assignment_params.J_align_length;
  _assigns.VDJ(in,0)=assignment_params.v_g;
  _assigns.VDJ(in,1)=assignment_params.d_g;
  _assigns.VDJ(in,2)=assignment_params.j_g; //% gene choices

  cout<<"checkpoint L2904"<<endl;

  _assigns.pVdelV(in,0)=assignment_params.npV;
  _assigns.pVdelV(in,1)=assignment_params.ndV; //% palindromes and deletions
  _assigns.pJdelJ(in,0)=assignment_params.npJ;
  _assigns.pJdelJ(in,1)=assignment_params.ndJ;
  
  _assigns.delVinsVD(in,0) = assignment_params.ndV; _assigns.delVinsVD(in,1)= assignment_params.niVD;//in matlab it niVD
  _assigns.delVinsDJ(in,0) = assignment_params.ndV; _assigns.delVinsDJ(in,1)= assignment_params.niDJ;//in matlab it is niDJ
  _assigns.delVdelJ(in,0) =  assignment_params.ndV; _assigns.delVdelJ(in, 1)= assignment_params.ndJ;

  _assigns.delJinsVD(in,0) = assignment_params.ndJ; _assigns.delJinsVD(in,1)= assignment_params.niVD;//in matlab it is niVD
  _assigns.delJinsDJ(in,0) = assignment_params.ndJ; _assigns.delJinsDJ(in, 1)=assignment_params.niDJ;//in matlab it is niDJ
                                                    
  _assigns.insVDinsDJ(in,0) = assignment_params.niVD; _assigns.insVDinsDJ(in,1)= assignment_params.niDJ;//ni DJ in matlab
  
  _assigns.insDJ_D_align_length(in,0) = assignment_params.niDJ; //niDJ in matlab
  _assigns.insDJ_D_align_length(in,1)= assignment_params.D_align_length;
  
  _assigns.insVD_D_align_length(in,0) = assignment_params.niVD;//niVD in matlab
  _assigns.insVD_D_align_length(in,1)= assignment_params.D_align_length;
                                                    
  _assigns.insDJ_J_align_length(in,0) = assignment_params.niDJ;//ni in matlab
  _assigns.insDJ_J_align_length(in,1)= assignment_params.J_align_length;
  
  _assigns.insVD_J_align_length(in,0) = assignment_params.niVD; //niVD in matlab
  _assigns.insVD_J_align_length(in,1)= assignment_params.J_align_length;
                                                    
  _assigns.insDJ_V_align_length(in,0) = assignment_params.niDJ;//niDJ in matlab
  _assigns.insDJ_V_align_length(in,1) = assignment_params.V_align_length;
  
  _assigns.insVD_V_align_length(in,0) = assignment_params.niVD; //niVD in matlab
  _assigns.insVD_V_align_length(in, 1)= assignment_params.V_align_length;
                                                    
  _assigns.Dallele_D_align_length(in,0) = assignment_params.d_a;_assigns.Dallele_D_align_length(in, 1)=  assignment_params.D_align_length;
  //count<<"vdj palindrome 56"<<endl;                                                                          
  _assigns.VdelV(in,0) = assignment_params.v_g; _assigns.VdelV(in,1)=assignment_params.ndV ;
  _assigns.JdelJ(in,0) = assignment_params.j_g;_assigns.JdelJ(in,1)=assignment_params.ndJ;
  _assigns.VinsVD(in,0) =assignment_params.v_g; _assigns.VinsVD(in,1)=assignment_params.niVD;//niVd in matlab
  _assigns.DinsVD(in,0) = assignment_params.d_g; _assigns.DinsVD(in,1)=assignment_params.niVD;//niVD in matlab
  _assigns.DinsDJ(in,0) = assignment_params.d_g; _assigns.DinsDJ(in,1)= assignment_params.niDJ;//niDJ in matlab
  _assigns.JinsDJ(in,0) = assignment_params.j_g; _assigns.JinsDJ(in, 1)=assignment_params.niDJ;//niDJ in matlab
  _assigns.VdelJ(in,0) = assignment_params.v_g; _assigns.VdelJ(in,1)= assignment_params.ndJ;
  _assigns.JdelV(in,0) = assignment_params.j_g; _assigns.JdelJ(in,1)=assignment_params.ndV;
                                                    
  _assigns.DdelV(in,0) = assignment_params.d_g; _assigns.DdelV(in, 1)= assignment_params.ndV;
  _assigns.DdelJ(in,0) = assignment_params.d_g; _assigns.DdelJ(in,1)= assignment_params.ndJ;
  _assigns.VinsDJ(in,0) = assignment_params.v_g;_assigns.VinsDJ(in,1)=assignment_params.niDJ;//niDJ in matlab
  _assigns.JinsVD(in,0) = assignment_params.j_g;_assigns.JinsVD(in,1)= assignment_params.niVD;//niVD in matlab
  
  _assigns.pVinsVD(in,0) = assignment_params.npV; _assigns.pVinsVD(in,1) = assignment_params.niVD;
  _assigns.pVinsDJ(in,0) = assignment_params.npV; _assigns.pVinsDJ(in,1) = assignment_params.niDJ;
  _assigns.pVdelJ(in,0) = assignment_params.npV; _assigns.pVdelJ(in, 1)= assignment_params.ndJ;
  _assigns.VpV(in,0) = assignment_params.v_g; _assigns.VpV(in,1)=assignment_params.npV;
  _assigns.JpV(in,0) = assignment_params.j_g; _assigns.JpV(in,1) = assignment_params.npV;
  _assigns.DpV(in,0) = assignment_params.d_g; _assigns.DpV(in,1)= assignment_params.npV;
    
  _assigns.pJinsVD(in,0) = assignment_params.npJ;_assigns.pJinsVD(in,1)= assignment_params.niVD;
  _assigns.pJinsDJ(in,0) = assignment_params.npJ; _assigns.pJinsDJ(in,1) = assignment_params.niDJ;
  _assigns.pJdelV(in,0) = assignment_params.npJ; _assigns.pJdelV(in,1) = assignment_params.ndV;
  _assigns.VpJ(in,0) = assignment_params.v_g; _assigns.VpJ(in,1) = assignment_params.npJ;
  _assigns.JpJ(in,0) = assignment_params.j_g; _assigns.JpJ(in,1)=assignment_params.npJ;
  _assigns.DpJ(in,0) = assignment_params.d_g; _assigns.DpJ(in,1)=assignment_params.npJ;
  _assigns.pVpJ(in,0)  = assignment_params.npV; _assigns.pVpJ(in,1)=assignment_params.npJ;
  _assigns.VcutV(in,0) = assignment_params.v_g; _assigns.VcutV(in,1)=assignment_params.ndV-assignment_params.npV-_model.min_V_cut;//this is same as cutV_given_V;!!!
  
  _assigns.JcutJ(in,0) = assignment_params.j_g; _assigns.JcutJ(in,1)=assignment_params.ndJ-assignment_params.npJ-_model.min_J_cut;
  _assigns.DcutV(in,0) = assignment_params.d_g; _assigns.DcutV(in,1)=assignment_params.ndV-assignment_params.npV - _model.min_V_cut;
  _assigns.DcutJ(in,0) = assignment_params.d_g; _assigns.DcutJ(in,1) = assignment_params.ndJ-assignment_params.npJ - _model.min_J_cut ;
  _assigns.VcutJ(in,0) = assignment_params.v_g; _assigns.VcutJ(in, 1)= assignment_params.ndJ-assignment_params.npJ -_model.min_J_cut;
  _assigns.JcutV(in,0) = assignment_params.j_g; _assigns.JcutV(in,1)=assignment_params.ndV-assignment_params.npV - _model.min_V_cut;
  _assigns.insVDcutV(in,0) = assignment_params.niVD; _assigns.insVDcutV(in,1)=assignment_params.ndV-assignment_params.npV -_model.min_V_cut;
  _assigns.insDJcutV(in,0) = assignment_params.niDJ; _assigns.insDJcutV(in,1)=assignment_params.ndV-assignment_params.npV - _model.min_V_cut;
  _assigns.insVDcutJ(in,0) = assignment_params.niVD; _assigns.insVDcutJ(in,1)=assignment_params.ndJ-assignment_params.npJ - _model.min_J_cut;
  _assigns.insDJcutJ(in,0) = assignment_params.niDJ; _assigns.insDJcutJ(in,1)=assignment_params.ndJ-assignment_params.npJ - _model.min_J_cut;
                                                    
  _assigns.cutDlcutDr_given_D(in,0) =assignment_params.d_g;_assigns.cutDlcutDr_given_D(in,1)= assignment_params.ndDl-assignment_params.npDl -_model.min_D_cut; _assigns.cutDlcutDr_given_D(in,2)=assignment_params.ndDr-assignment_params.npDr  -_model.min_D_cut;
               
  //***need to be careful about above or overall assignment. Do I really need to change the matrix orientation to make it fit the p1_given_p2. is it necessary?                                   
  //% All possible pairwise joint distributions, just for tracking.
  
  _assigns.pDldelDl(in,0)=assignment_params.npDl; _assigns.pDldelDl(in,1)= assignment_params.ndDl;
  _assigns.pDrdelDr(in,0)=assignment_params.npDr; _assigns.pDrdelDr(in,1)=assignment_params.ndDr;
  _assigns.delVdelDl(in,0) = assignment_params.ndV;_assigns.delVdelDl(in,1) = assignment_params.ndDl;
  _assigns.delVdelDr(in,0) = assignment_params.ndV; _assigns.delVdelDr(in,1) = assignment_params.ndDr;
                                                    
  _assigns.delJdelDl(in,0) = assignment_params.ndJ; _assigns.delJdelDl(in,1) = assignment_params.ndDl;
  _assigns.delJdelDr(in,0) = assignment_params.ndJ; _assigns.delJdelDr(in,1) = assignment_params.ndDr;
  _assigns.delDlinsVD(in,0) = assignment_params.ndDl; _assigns.delDlinsVD(in,1) = assignment_params.niVD;
  _assigns.delDlinsDJ(in,0) = assignment_params.ndDl; _assigns.delDlinsDJ(in,1) = assignment_params.niDJ;
  _assigns.delDldelDr(in,0) = assignment_params.ndDl; _assigns.delDldelDr(in,1) = assignment_params.ndDr;
  _assigns.delDrinsVD(in,0) = assignment_params.ndDr; _assigns.delDrinsVD(in,1) = assignment_params.niVD;
  _assigns.delDrinsDJ(in,0) = assignment_params.ndDr; _assigns.delDrinsDJ(in,1)=assignment_params.niDJ;
  _assigns.DdelDl(in,0) = assignment_params.d_g; _assigns.DdelDl(in,1) = assignment_params.ndDl;
  _assigns.DdelDr(in,0) = assignment_params.d_g; _assigns.DdelDr(in,1) = assignment_params.ndDr ;

  _assigns.VdelDl(in,0) = assignment_params.v_g; _assigns.VdelDl(in,1) = assignment_params.ndDl;
  _assigns.VdelDr(in,0) = assignment_params.v_g; _assigns.VdelDr(in,1) = assignment_params.ndDr;
  _assigns.JdelDl(in,0) = assignment_params.j_g; _assigns.JdelDl(in,1) = assignment_params.ndDl;
  _assigns.JdelDr(in,0) = assignment_params.j_g; _assigns.JdelDr(in,1) = assignment_params.ndDr;
  _assigns.pVdelDl(in,0) = assignment_params.npV; _assigns.pVdelDl(in,1) = assignment_params.ndDl;
  _assigns.pVdelDr(in,0) = assignment_params.npV; _assigns.pVdelDr(in,1) = assignment_params.ndDr;
  _assigns.pJdelDl(in,0) = assignment_params.npJ;  _assigns.pJdelDl(in,1) = assignment_params.ndDl;
  _assigns.pJdelDr(in,0) = assignment_params.npJ; _assigns.pJdelDr(in,1) = assignment_params.ndDr;
  _assigns.pDlinsVD(in,0) = assignment_params.npDl;  _assigns.pDlinsVD(in,1) = assignment_params.niVD;
  _assigns.pDlinsDJ(in,0) = assignment_params.npDl; _assigns.pDlinsDJ(in,1) = assignment_params.niDJ;
  _assigns.pDldelV(in,0) = assignment_params.npDl; _assigns.pDldelV(in,1) = assignment_params.ndV;
  _assigns.pDldelJ(in,0)  = assignment_params.npDl; _assigns.pDldelJ(in,1)  = assignment_params.ndJ;
  _assigns.pDldelDr(in,0) = assignment_params.npDl; _assigns.pDldelDr(in,1) = assignment_params.ndDr;
  _assigns.VpDl(in,0) = assignment_params.v_g; _assigns.VpDl(in,1) = assignment_params.npDl;
  _assigns.JpDl(in,0) = assignment_params.j_g; _assigns.JpDl(in,1) = assignment_params.npDl;
  _assigns.DpDl(in,0) = assignment_params.d_g; _assigns.DpDl(in,1)=assignment_params.npDl;
                                                    
  _assigns.pDrinsVD(in,0) = assignment_params.npDr; _assigns.pDrinsVD(in,1) =assignment_params.niVD;
  _assigns.pDrinsDJ(in,0) = assignment_params.npDr; _assigns.pDrinsDJ(in,1) = assignment_params.niDJ;
  _assigns.pDrdelV(in,0)=  assignment_params.npDr; _assigns.pDrdelV(in,1)=  assignment_params.ndV;
  _assigns.pDrdelJ(in,0)=  assignment_params.npDr; _assigns.pDrdelJ(in,1)=  assignment_params.ndJ;
  _assigns.pDrdelDl(in,0) = assignment_params.npDr;  _assigns.pDrdelDl(in,1) = assignment_params.ndDl;
  _assigns.VpDr(in,0) = assignment_params.v_g; _assigns.VpDr(in,1) = assignment_params.npDr;
  _assigns.JpDr(in,0) = assignment_params.j_g; _assigns.JpDr(in,1) = assignment_params.npDr;
  _assigns.DpDr(in,0) = assignment_params.d_g; _assigns.DpDr(in,1) = assignment_params.npDr;
                                                    
  _assigns.pVpDl(in,0)  = assignment_params.npV; _assigns.pVpDl(in,1)  = assignment_params.npDl;
  _assigns.pVpDr(in,0)  = assignment_params.npV; _assigns.pVpDr(in,1)  = assignment_params.npDr;
  _assigns.pDlpDr(in,0) = assignment_params.npDl; _assigns.pDlpDr(in,1) = assignment_params.npDr;
  _assigns.pDlpJ(in,0) = assignment_params.npDl; _assigns.pDlpJ(in,1) = assignment_params.npJ;
  _assigns.pDrpJ(in,0) = assignment_params.npDr;  _assigns.pDrpJ(in,1) = assignment_params.npJ;
                                                    
  _assigns.DcutDl(in,0) = assignment_params.d_g; _assigns.DcutDl(in,1) = assignment_params.ndDl-assignment_params.npDl-_model.min_D_cut;
  _assigns.DcutDr(in,0) = assignment_params.d_g; _assigns.DcutDr(in,1) = assignment_params.ndDr-assignment_params.npDr - _model.min_D_cut;
  _assigns.VcutDl(in,0) = assignment_params.v_g; _assigns.VcutDl(in,1) = assignment_params.ndDl-assignment_params.npDl-_model.min_D_cut;
  _assigns.VcutDr(in,0) = assignment_params.v_g; _assigns.VcutDr(in,1) = assignment_params.ndDr-assignment_params.npDr - _model.min_D_cut;
  _assigns.JcutDl(in,0) = assignment_params.j_g; _assigns.JcutDl(in,1) = assignment_params.ndDl-assignment_params.npDl - _model.min_D_cut;
  _assigns.JcutDr(in,0) = assignment_params.j_g; _assigns.JcutDr(in,1) = assignment_params.ndDr-assignment_params.npDr - _model.min_D_cut;
  _assigns.insVDcutDl(in,0) = assignment_params.niVD;_assigns.insVDcutDl(in,1) = assignment_params.ndDl-assignment_params.npDl -_model.min_D_cut ;
  _assigns.insDJcutDl(in,0) = assignment_params.niDJ;  _assigns.insDJcutDl(in,1) = assignment_params.ndDl-assignment_params.npDl - _model.min_D_cut ;
  _assigns.insVDcutDr(in,0) = assignment_params.niVD;  _assigns.insVDcutDr(in,1) = assignment_params.ndDr - assignment_params.npDr - _model.min_D_cut;
  _assigns.insDJcutDr(in,0) = assignment_params.niDJ;  _assigns.insDJcutDr(in,0) = assignment_params.ndDr - assignment_params.npDr - _model.min_D_cut ;
                                                  
  _assigns.error_vs_position.SetSubMatrix(in, _model.model_params.maximum_read_length, error_vs_position);
  _assigns.coverage.SetSubMatrix(in, _model.model_params.maximum_read_length,  coverage);
  
  //%%% probability of this assignment
  cout<<"log_proba:"<<log_proba<<";log_proba_Rerror_normalization:"<<assignment_params.log_proba_Rerror_normalization<<endl;
  _assigns.proba(in)=exp(log_proba + assignment_params.log_proba_Rerror_normalization);
                                                    
  if (nerrors==0)
    {
      _assigns.event_probability(in) = exp(log_proba);
    }
  else
    _assigns.event_probability(in) = 0;
  
  cout<<"assigns.likelihood:"<<_assigns.likelihood<<endl;
  _assigns.generation_probability = _assigns.generation_probability + _assigns.event_probability(in);
  _assigns.likelihood = _assigns.likelihood + _assigns.proba(in);
  cout<<"assigns.likelihood:"<<_assigns.likelihood<<endl;
  
  //this following section for incrementing the counters were originally in Matlab somewhere
  //above, we moved them here to make it compatible for c++ code.
  in++;
  assignment_params.n_assignments_v_gene++;
  assignment_params.n_assignments_d_gene++;
  assignment_params.n_assignments_j_gene++;
  //in++;
  assignment_params.in=in;
  //?????need to check whether this is correct//_assigns.n_assignments=in; <==where we can set this????? ANSWER: we will set this up at the outermost loop, V loop

  //clean up
  delete[] v_err_pos_rel;
  delete[] j_err_pos_rel;
  delete[] error_vs_position;
  delete[] coverage;  
  
  //% If maximum number of valid assignments is reached, stop.
  //% If we are not smoothing, stop after 1 assignment
  //% Usually, we don't reach maximum number of assignments because we break out
  //% of all the loops due to the various thresholds below.
  if (in > _model.max_assignments || (~_do_smoothing && in>1))
    {
      _assigns.n_assignments = in;
      _assigns.skips = assignment_params.skips;
      
      assignment_params.v_break_out=true;//we add this because we are done with the assignment
      //need to jump all out.
      return false; //here we are done with overall assignment, 
    }//end
                                                    
if( _force_all_alleles && assignment_params.n_assignments_v_gene > ((double)_model.max_assignments)/_V.numOfAligned)
    {
      assignment_params.v_break_out = true;
      //break;<--break out since there are too many assignments for this segments
      return false;
    }
  //end
//cout<<"end here"<<endl;
if (_force_all_alleles && assignment_params.n_assignments_j_gene > ((double)_model.max_assignments)/_J.numOfAligned)
    {
      assignment_params.j_break_out = true;
      //break; <+++in matlab code this is break, but break is not working here, so we change it to return. 
      return false;
    }//end
//cout<<"in the middel"<<endl;
 if (_force_all_alleles && assignment_params.n_assignments_d_gene > ((double)_model.max_assignments)/_D.n_D_alleles)
    {
      assignment_params.d_break_out = true;
      //break; see above block
      return false;
    }//end
                                                    
 cout<<"exiting funcgtion"<<endl;


  return true;
}
