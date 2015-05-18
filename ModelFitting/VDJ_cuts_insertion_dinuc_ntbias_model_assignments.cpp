
#include <cmath>
#include "VDJ_cuts_insertion_dinuc_ntbias_model_assignments.hpp"

#include "VDJ_model_assignments_settings.hpp"
#include "../SIGPIG/AlignmentSettings.hpp"

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
  assignment_params.in=0;
  assignment_params.skips=0;
  //cout<<"&&&inside assignments:2"<<endl;
  assignment_params.max_J_depth=AlignmentSettings::max_J_length;
  assignment_params.max_V_depth=_seq.GetLength()+ AlignmentSettings::max_V_length;
  
  assignment_params.J_max_error=0;
  assignment_params.D_max_error=0;

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
  assignment_params.L_err=_seq.GetLength();
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
      assignment_params.nd_start=-1* _model.negative_excess_deletions_max;
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
  assign_VDJ_alleles(_model, _seq, _V, _D, _J, _genV, _numV, _genD, _numD,
		     _genJ, _numJ, _no_error, _ignore_deep_error, _do_smoothing,
		     _force_all_alleles, _READ_LENGTH_CORRECTION,
		     assignment_params, _assigns
		     );
    
  return true;

}
			   

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
  cout<<"_V.numOfAligned:"<<_V.numOfAligned<<endl;
  for(unsigned v=0;v<_V.numOfAligned;v++)
    {
      //cout<<"---->>>>loop v :"<<v<<endl;
      assignment_params.v=v;
      //check for validity
      if(assignment_params.skips>assignment_params.max_skips)
	{
	  _assigns.n_assignments=assignment_params.in;
	  _assigns.skips=assignment_params.skips;
	  assignment_params.v_break_out=true;
	  return false;
	}
      //cout<<"\tinside VDJ:1"<<endl;
      assignment_params.v_a=_V.alleles_all[v];
      assignment_params.v_g=_genV[assignment_params.v_a].Get_GeneIndex();

      //for high_error_region, we don't use in my code,
      //but keep it for now anyway
      if(_ignore_deep_error)
	assignment_params.high_error_region=_model.high_error_region;
      else
	assignment_params.high_error_region=_model.high_error_region;
//count<<"\tinside VDJ:2"<<endl;
      if(v==0)//for first round, set some starting point
	{
	  cout<<"in v=0"<<endl;
	  dim_size[0]=_V.n_errors[v];
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
	      //cout<<"too many errors compared with the best V"<<endl;
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
      for(unsigned j=0;j<_J.numOfAligned;j++)
	{
	  //cout<<"\tloop j:"<<j<<endl;
	  assignment_params.j=j;
	  if(assignment_params.skips>assignment_params.max_skips)
	    {
	      _assigns.n_assignments=assignment_params.in;
	      _assigns.skips=assignment_params.skips;
	      assignment_params.v_break_out=true;
	      return false; //=====is true good here
	    }
	  //cout<<"\tinside VDJ:j-1"<<endl;
	  if(_J.align_length[j]<_model.min_J_align_length)
	    continue;
	  assignment_params.j_a=_J.alleles_all[j];
	  assignment_params.j_g=_genJ[assignment_params.j_a].Get_GeneIndex();
	  
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
	      //cout<<"\t\tloop D:"<<d_i<<endl;
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
	      assignment_params.niVD_DJ_min=0;
	      //	      dim_pos[0]=v;
	      unsigned temp_int=assignment_params.niVD_DJ0- _D.align_length[d][0]-max_mf(_V.p_region_max_length[v], AlignmentSettings::V_maximum_deletion+1)-assignment_params.p_max_J - assignment_params.p_max_Dl- assignment_params.p_max_Dr;
	      if(temp_int>0)
		assignment_params.niVD_DJ_min=temp_int;
	      double log_p_max_nt_VD_DJ_d_allele=assignment_params.niVD_DJ_min*assignment_params.log_max_model_p_nt;
	      //cout<<"\t\tinside VDJ:D_4"<<endl;
	      //cout<<"log_probabase:"<<log(probabase)<<endl;
	      //cout<<"log_highest_probability:"<<assignment_params.log_highest_probability<<endl;
	      assignment_params.log_max_pcutD_loop_d=assignment_params.log_max_model_pcutD_given_gene(assignment_params.d_g);
	      //assignment_params.log_max_pcutVDJ_loop_d=assignment_params.log

	      if(probabase==0||
		 (log(probabase)+assignment_params.log_max_pcutV_loop_v+
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
		  if(assignment_params.v_break_out)//||assignment_param.j_break_out||assignment_param.d_break.out)
		    {
		      return false;
		    }
		  else
		    {
		      break;
		    }
		  //return false;
		  
		  }
	      if(assignment_params.v_break_out||assignment_params.j_break_out||assignment_params.d_break_out)
		{
		  break;
		}
	    }//end of d loop ====
	  if(assignment_params.v_break_out||assignment_params.j_break_out)
	    {
	      break;
	    }
	}//end of j loop =====
      if(assignment_params.v_break_out)
	{
	  break;
	}
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
  //count<<"<<<<inside DJ deletions:"<<endl;
  unsigned ndV_max=_model.max_excess_V_deletions;
  if(ndV_max>_V.align_length[assignment_params.v]-_model.min_V_length)
    {
      ndV_max=_V.align_length[assignment_params.v]-_model.min_V_length;
    }
  if(ndV_max>_model.max_V_deletions-_V.min_deletions[assignment_params.v])
    {
      ndV_max=_model.max_V_deletions-_V.min_deletions[assignment_params.v];
    }
  //loop over V deletions
  for(assignment_params.ndV1=assignment_params.nd_start;assignment_params.ndV1<(signed)ndV_max;assignment_params.ndV1++)
    {
      //cout<<"\t\t\tloop v_deletion:"<<assignment_params.ndV1<<endl;
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
	    assignment_params.nerrorsv = _V.n_errors[assignment_params.v] - sum_all_bool(assignment_params.v_err_pos > (_V.align_position[assignment_params.v][0]+_V.align_length[assignment_params.v]-1 - assignment_params.ndV1)) - sum_all_bool(assignment_params.v_err_pos > 0);
	  //            here change the very last part of the sum to not including
	  //the high error region. since we know we don't have this.
	}
                
      if( _no_error && assignment_params.nerrorsv > 0)
                    continue;
                           
      assignment_params.log_perrv = assignment_params.nerrorsv*assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
                
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
      assignment_params.niVD_DJ_min= assignment_params.niVD_DJ0+assignment_params.ndV1-_D.align_length[assignment_params.d][0] - assignment_params.npV_max - assignment_params.p_max_J-assignment_params.p_max_Dl-assignment_params.p_max_Dr;
      
      if(((signed)assignment_params.niVD_DJ_min) <0)
	{
	  assignment_params.niVD_DJ_min=0;
	}
      double log_p_max_nt_VD_DJ_dV = (assignment_params.niVD_DJ_min)*assignment_params.log_max_model_p_nt;
                //cout<<"\t\tV deletion:vd_3"<<endl;
      assignment_params.log_highest_probability_GIVEN_current_V_deletions = -1000;// % Highest probability so far GIVEN current V deletions AND outer loop variables {V, J, D alleles}                
      //cout<<"\t\tV deletion:vd_4"<<endl;
      double test_proba = assignment_params.log_probabase + assignment_params.log_perrv + (assignment_params.log_max_pcutV_loop_v+assignment_params.log_max_pcutJ_loop_j+assignment_params.log_max_pcutD_loop_d) + assignment_params.log_max_model_pins + log_p_max_nt_VD_DJ_dV;               
                
      if (test_proba < assignment_params.log_probability_threshold_factor + assignment_params.log_highest_probability)
	{
                    assignment_params.skips  ++;
                    continue;
	}
                
      assignment_params.V_align_length = _V.align_length[assignment_params.v] - assignment_params.ndV1 - _READ_LENGTH_CORRECTION;
      if(assignment_params.V_align_length>=300)
	{
	  assignment_params.v_break_out=true;
	  return false;
	}
//cout<<"\t\tV deletion:vd_6"<<endl;

      //======now start doing j deletions loop for
      //% Loop over J right deletions.
      //% ndJ1 is the deviation from the no. of deletions implied by the alignment (min_deletions).
      //% Lower bound is nd_start (usually -3).
      //% Upper bound is set by maximum deletions allowed and by minimum match length required
      unsigned temp_array[]={_model.max_excess_J_deletions,_J.align_length[assignment_params.j]-_model.min_J_assign_length , _model.max_J_deletions - _J.min_deletions[assignment_params.j]};
      unsigned ndJ_max = min_mf(temp_array, 3);
      for( assignment_params.ndJ1=assignment_params.nd_start;assignment_params.ndJ1<(signed)ndJ_max;assignment_params.ndJ1++)
	{	//		%ndJ1=nd_start
	  //count<<"\t\tloop j deletion"<<assignment_params.ndJ1<<endl;
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
	  if (assignment_params.ndJ < 0 || assignment_params.ndJ >(signed) _model.max_J_deletions)
	    {
	      //% Go to next iteration of loop.
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
	  //count<<"assignment_params.ndJ1:"<<assignment_params.ndJ1<<endl;
	  if (assignment_params.ndJ1 < 0)
	    {
	      //% Case of 'negative' deletions
	      assignment_params.j_ex_errs_i = j_err_excess_pos >= (_J.align_position[assignment_params.j][0] + assignment_params.ndJ1);
	      assignment_params.j_ex_errs = sum_all_bool(assignment_params.j_ex_errs_i);
	      if (assignment_params.j_ex_errs > 1)
		continue;
	      
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
                    
	  assignment_params.log_perrj = assignment_params.nerrorsj*assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
          
	  //%                     % Upper bound on insertions nt bias factor
	  assignment_params.niVD_DJ_min =  assignment_params.niVD_DJ0 + assignment_params.ndV1 + assignment_params.ndJ1 - _D.align_length[assignment_params.d][0] - _V.p_region_max_length[assignment_params.v][1 + assignment_params.ndV] - _J.p_region_max_length[assignment_params.j][1 + assignment_params.ndJ] - assignment_params.p_max_Dl - assignment_params.p_max_Dr;
	  if(assignment_params.niVD_DJ_min<0)
	    {
	      assignment_params.niVD_DJ_min=0;
	    }
	  double log_p_max_nt_VD_DJ_dJ = (assignment_params.niVD_DJ_min)*assignment_params.log_max_model_p_nt;
          
	  assignment_params.log_highest_probability_GIVEN_current_J_deletions = -1000;                   //cout<<"\t\t\tJ deletion:J_5"<<endl;
	  
	  double test_proba = assignment_params.log_probabase + assignment_params.log_perrv + assignment_params.log_perrj + (assignment_params.log_max_pcutV_loop_v+assignment_params.log_max_pcutJ_loop_j+assignment_params.log_max_pcutD_loop_d) + assignment_params.log_max_model_pins + log_p_max_nt_VD_DJ_dJ;
	            
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
				   _no_error, _ignore_deep_error, _do_smoothing,
				   _force_all_alleles, _READ_LENGTH_CORRECTION,
				   assignment_params, _assigns)
	     )
	    {
	      if(assignment_params.v_break_out||assignment_params.j_break_out||assignment_params.d_break_out)
		{
		  return false;
		}
	      else
		{
		  break;
		}
	      //return false;
	      }

	}//end of J deletion loop
    }//end of V deletions loop for
  
  return true;
}//end of assign_V_deletions

//=======================================
//***************************************
bool assign_VJ_palindrome
(VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 /*const double& _probability_threshold_factor,*/ const bool& _no_error,
 const bool& _ignore_deep_error, const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
 /*output, input*/VDJ_model_assignments_settings& assignment_params,
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns)
{
  //set the direction for VJ_palindrome search, IS THIS REALLY NECESSARY or IMPORTANT???
  //do it anyway for now
  //assignment_params.nd_start=0;
  assignment_params.np_start_from_max=true;
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
    
  //now looping the possible cases of npV first
  for(npV=npV_start;;npV+=np_step)
    {
      if(assignment_params.np_start_from_max&&npV>npV_end)
	{
	  break; //we are done
	}
      if(!assignment_params.np_start_from_max&&npV<npV_end)
	{
	  break;//we are done.
	}
      
      if (assignment_params.skips > assignment_params.max_skips)
	{
	  _assigns.n_assignments = assignment_params.in;
	  _assigns.skips =  assignment_params.skips;
	  assignment_params.v_break_out=true;//===in this case, we want to jump all out.
	  return false;
	}//       end
      
      int v_end=_V.align_positions(assignment_params.v, 0)+_V.align_length(assignment_params.v)-1-assignment_params.ndV1+npV;
      if(v_end>=_J.align_position(assignment_params.j, 0)+assignment_params.ndj1)
	//V with palindrome overlaps with J, not possible
	{
	  continue;
	}

      //now start doing the cut variable, Note: cut variable is the observed deletion since P_nt often is confused with negative deletions.
      int ncutV=assignment_params.ndV-npV;
      int PcutV=_model.PcutV_given_V(assignment_params.v_g, ncutV-_model.min_V_cut);
      if(PcutV==0)
	{
	  assignment_params.skips++;
	  continue;
	}
      assignment_params.log_PcutV=log(PcutV);
      double log_max_pcutVDJ_loop_pV=log_PcutV+log_max_pcutD_loop_d+log_max_pcutJ_loop_j;
      
      //upper bound on insertions nt bias factor
      int niVD_DJ_min=(signed)(assignment_params.niVD_DJ0)+assignment_params.ndV1+assignment_params.ndJ1-_D.align_length(d,0)-npV-_J.p_region_max_length(j, ndJ)-assignment_params.p_maxDl-assignment_params.p_maxDr;
      if( niVD_DJ_min<0)
	{
	  niVD_D=0;
	}
      double log_p_max_nt_VD_DJ_pV=niVD_DJ_min*assignment_params.log_max_model_p_nt;
      
      assignment_params.test_proba=
	assignment_params.log_probabase + assignment_params.log_perrorj
	+ assignment_params.log_perrorv+ log_max_pcutVDJ_loop_pV+
	assignment_params.log_max_model_pins + log_p_max_nt_VD_DJ_pV;
      if(assignment_params.test_proba< assignment_params.log_probability_threshold_factor+
	 log_highest_probability)
	{
	  assignment_params.skips++;
	  continue;
	}
	 
      //next start doing the loop over half-length of J palindrome
      for(npJ=npJ_start;;npJ+=np_step)
	{
	  if(assignment_params.np_start_from_max&&npJ>npJ_end)
	    {
	      break; //we are done
	    }
	  if(!assignment_params.np_start_from_max&&npJ<npJ_end)
	    {
	      break;//we are done.
	    }
	  
	  if (assignment_params.skips > assignment_params.max_skips)
	    {
	      _assigns.n_assignments = assignment_params.in;
	      _assigns.skips =  assignment_params.skips;
	      assignment_params.v_break_out=true;//===in this case, we want to jump all out.
	      return false;
	    }//       end
	  
	  int J_start=_J.align_position(assignment_params.j,0)+assignment_params.ndJ1-npJ;

	  if(v_end>=J_start)
	    //V with palindrome overlaps with J, not possible
	    {
	      continue;
	    }
	  
	  //now cut variable from J side
	  int ncutJ=assignment_params.ndJ-npJ;
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
	  assignment_params.n_D_aligns=_D.n_alignments(assignment_params.d);
	  assignment_params.start_n_D_aligns=0;
	  
	  //find the first completely valid alignment
	  bool completely_valid_na=false;
	  unsigned c_na=0;

	  while(!completely_valid&&c_na<assignment_params.n_D_aligns)
	    {
	      completely_valid_na=_D.align_position_lef(assignment_params.d, c_na)>V_end 
		&&_D.align_position_right(assignment_params.d, c_na)<J_start;
	      c_na++;
	    }//end of while for completely_valid_na;
	  unsigned first_valid_length=0;

	  if(completely_valid_na)
	    {
	      first_valid_length=_D.align_length(d, c_na_)
	    }
	  else
	    {
	      //find the first partly valid alignment
	      bool partly_valid_na=false;
	      unsigned p_na=0;
	      while(~partly_valid_na&&p_na<assignment_params.n_D_aligns)
		{
		  partly_valid_na=_D.align_position_left(d, p_na)<J_start 
		    || _D.align_position_right(d, p_na)>V_end;
		  p_na++;
		}//end of partly_valid_na while loop

	      if(partly_valid_na)
		{
		  int p_valid_start=((J_start-1)>(_D.align_position_right(d, p_na)))?_D.align_position_right(d, p_na):(J_start-1);
		  int p_valid_end=((V_end+1)>_D.align_position_left(d, p_na))?(V_end+1):_D.align_position_left(d,p_na);
		  first_valid_length=p_valid_start-p_valid_end+1;
		  
		  assignment_params.start_n_D_aligns=p_na;
		}
	      else //no partially valid alignment
		{
		  first_valid_length=0;
		  //do only zero D. no even partly valid alignment.
		  assignment_params.start_n_D_aligns=assignment_params.n_D_aligns+1;
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
	  assignment_params.niVD_DJ_total=J_start-V_end-1;
	  //upper bound for nt bias factor
	  int niVD_DJ_min= (signed)niVD_D_total-frist_valid_length-
	    _D.p_region_max_length_left(d, assignment_params.start_n_D_aligns, _D.deletions_left(d, assignment_params.start_n_D_aligns))-
	    _D.p_region_max_length_right(d, assignment_params.start_n_D_aligns, _D.deletions_right(d, assignment_params.start_n_D_aligns));
	  if(niVD_DJ_min<0)
	    niVD_DJ_min=0;
			
	  double log_p_max_nt_VD_DJ_pJ=niVD_DJ_min*assignment_params.log_max_model_p_nt;
	  
	  assignment_params.test_proba=assigmnent_params.log_probabase+
	    assigmnent_params.log_perrv+log_max_pcutVDJ_loop_pJ+
	    assigmnent_params.log_max_model_pins+log_p_max_nt_VD_DJ_pJ;
	  
	  if(assignment_params.test_proba< assignment_params.log_probability_threshold_factor+
	     log_highest_probability)
	    {
	      assignment_params.skips++;
	      continue;
	    } 
	  
	  //---------> start doing d alignment from here

	}//end of npJ loop 
      
    }//end of npV loop palindrome
  //Note for myself, we probably we include the blow code inside the above blocks of for loops
  
  //***************************
  //===========>Good code below from previous version, need to be carefully
  //in here we start doing the stats!! set things to the assigns.

  //Increment valid assignment number
  //assignment_params.in+=1;//disp(['assignment=', num2str(in)]) <==========here we move the increment to the end,
  assignment_params.zeroD=false;
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
   double log_proba = assignment_params.log_probabase + log_perr + /*log_PcutVDJ +*/ 
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
  /*if(log_proba > log_highest_probability_GIVEN_current_Dl_deletions) 
    log_highest_probability_GIVEN_current_Dl_deletions = log_proba; //end;
  if(log_proba > log_highest_probability_GIVEN_current_Dr_deletions)
    log_highest_probability_GIVEN_current_Dr_deletions = log_proba; //end;
  */
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
  if((signed)assignment_params.D_align_length<0)
    {
      assignment_params.D_align_length=0;
    }
  if((signed)assignment_params.genic_length<0)
    {
      //assignment_params.zeroD=true;
      //jump out this case;
      return true; //but don't affect the next case.
    }
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
	  /*if ndDl1 < 0 && d_ex_errs_left > 0
		     error_vs_position(d_err_excess_pos_left(d_ex_errs_left_i)) = 1;
	  end
	    
	    if ndDr1 < 0 && d_ex_errs_right > 0
		       error_vs_position(d_err_excess_pos_right(d_ex_errs_right_i)) = 1;
	  end
	    
	    if d_errs > 0
			  error_vs_position(d_err_pos(d_errs_i)) = 1;
			  end*/
	}//end
							  
    }//  end of _no_error if blaock
  
  //% Store everything in assigns.
  //% Each variable in counter that begins with an 'nP' or 'nM' must be
  //% found in assigns without the 'nP' or 'nM'.
  //count<<"\t start updating >>>>>>>>>>>>"<<endl;
  //%%% The following are model variables
  unsigned in=assignment_params.in;                                                  
  _assigns.V(in) = assignment_params.v_g; //in-1, since we have incremented the 
  _assigns.DJ(in,0)=assignment_params.d_g;_assigns.DJ(in,1)=assignment_params.j_g; //% gene choices
  _assigns.Vallele_given_gene(in,0) =assignment_params.v_g;_assigns.Vallele_given_gene(in,1)= _genV[assignment_params.v_a].Get_Allele(); //% allele choice given gene
  _assigns.Dallele_given_gene(in,0) =assignment_params.d_g;_assigns.Dallele_given_gene(in,1)= _genD[assignment_params.d_a].Get_Allele();
  //_assigns.cutV_given_V(in,:) = [ 1 - model.min_V_cut + ncutV, v_g];
  //assigns.cutJ_given_J(in,:) = [1 - model.min_J_cut + ncutJ, j_g];
  _assigns.insVD(in)=inVD;// % insertions, niVD in matlab
  _assigns.insDJ(in)=inDJ;// niDJ in matlab
  //assigns.nucleotideVD(in,:,:) = nucleotideVD;
  //assigns.nucleotideVD_5prime(in,:,:) = nucleotideVD_5prime;
  //assigns.nucleotideDJ(in,:,:) = nucleotideDJ;
  //assigns.nucleotideDJ_3prime(in,:,:) = nucleotideDJ_3prime;
  
  _assigns.error(in) = nerrors;//assignment_params.nerrorsv+assignment_params.nerrorsj+assignment_params.nerrorsd; //% numerator for error rate estimate
  //nerrors=_assigns.error(in);
  _assigns.sequenced_nucleotide(in) = assignment_params.genic_length; //% denominator for error rate estimate
  //% For tracking                                                  
  //% numbers of A,C,G and T in insertions: numerator for probabilities of A,C,G,T insertions
  //assigns.mononucleotideVD(in,:) = mononucleotideVD;
  //assigns.insertionVD(in,:) = [ niVD , niVD, niVD, niVD]; //% denominator for probabilities of A,C,G,T insertions.
  //assigns.mononucleotideDJ(in,:) = mononucleotideDJ;
  //assigns.insertionDJ(in,:) = [niDJ , niDJ, niDJ, niDJ]; % denominator for probabilities of A,C,G,T insertions.
                                                    
  //assigns.VD_left_edge_dinucleotide(in,:) = VD_left_edge_dinucleotide;
  //assigns.VD_right_edge_dinucleotide(in,:) = VD_right_edge_dinucleotide;
                                                    
  //assigns.DJ_left_edge_dinucleotide(in,:) = DJ_left_edge_dinucleotide;
  //assigns.DJ_right_edge_dinucleotide(in,:) = DJ_right_edge_dinucleotide;
                                                    
  //assigns.trinucleotideVD(in,:,:,:) = trinucleotideVD;
  //assigns.trinucleotideDJ(in,:,:,:) = trinucleotideDJ;
  //????????????????iiiii                                             
  //assigns.VV_err_pos(in,v_g,:) = v_err_pos_rel;
  _assigns.VV_align_length(in,0) = assignment_params.v_g;
  _assigns.VV_align_length(in,1) =assignment_params.V_align_length; //no need to be 1+V_align_length
                                                    
  //_assigns.JJ_err_pos(in,j_g,:) = j_err_pos_rel;
  _assigns.JJ_align_length(in,0) = assignment_params.j_g;
  _assigns.JJ_align_length(in,1)=assignment_params.J_align_length; //no need to be 1+J_align_length
                                                    
  
  //assigns.pVmax_delV_V(in,:) = [1 + npV_potential_max, 1+ ndV, v_g];
  //assigns.pJmax_delJ_J(in,:) = [1 + npJ_potential_max, 1+ ndJ, j_g];
                                                    
                                                    
  //assigns.pDlmax_delDl_D(in,:) = [1 + npDl_potential_max, 1+ ndDl, d_g];
  //assigns.pDrmax_delDr_D(in,:) = [1 + npDr_potential_max, 1+ ndDr, d_g];
                                                    
  _assigns.zeroD(in) = assignment_params.zeroD;
  _assigns.V_align_length(in) = assignment_params.V_align_length;
  _assigns.D_align_length(in) = assignment_params.D_align_length;
  _assigns.J_align_length(in) = assignment_params.J_align_length;
  _assigns.VDJ(in,0)=assignment_params.v_g;
  _assigns.VDJ(in,1)=assignment_params.d_g;
  _assigns.VDJ(in,2)=assignment_params.j_g; //% gene choices
  //assigns.pVdelV(in,:)=[ 1 + npV , 1 + ndV]; % palindromes and deletions
  //assigns.pJdelJ(in,:)=[ 1 + npJ , 1 + ndJ];
  _assigns.delVinsVD(in,0) = assignment_params.ndV; _assigns.delVinsVD(in,1)= inVD;//in matlab it niVD
  _assigns.delVinsDJ(in,0) = assignment_params.ndV; _assigns.delVinsDJ(in,1)=  inDJ;//in matlab it is niDJ
  _assigns.delVdelJ(in,0) =  assignment_params.ndV; _assigns.delVdelJ(in, 1)= assignment_params.ndJ;
  _assigns.delJinsVD(in,0) = assignment_params.ndJ; _assigns.delJinsVD(in,1)=  inVD;//in matlab it is niVD
  _assigns.delJinsDJ(in,0) = assignment_params.ndJ; _assigns.delJinsDJ(in, 1)=  inDJ;//in matlab it is niDJ
                                                    
  _assigns.insVDinsDJ(in,0) =  inVD; _assigns.insVDinsDJ(in,1)= inDJ;//ni DJ in matlab
  
  _assigns.insDJ_D_align_length(in,0) =  inDJ; //niDJ in matlab
  _assigns.insDJ_D_align_length(in,1)= assignment_params.D_align_length;
  
  _assigns.insVD_D_align_length(in,0) = inVD;//niVD in matlab
  _assigns.insVD_D_align_length(in,1)= assignment_params.D_align_length;
                                                    
  _assigns.insDJ_J_align_length(in,0) =  inDJ;//ni in matlab
  _assigns.insDJ_J_align_length(in,1)= assignment_params.J_align_length;
  
  _assigns.insVD_J_align_length(in,0) = inVD; //niVD in matlab
  _assigns.insVD_J_align_length(in,1)= assignment_params.J_align_length;
                                                    
  _assigns.insDJ_V_align_length(in,0) = inDJ;//niDJ in matlab
  _assigns.insDJ_V_align_length(in,1) = assignment_params.V_align_length;
  
  _assigns.insVD_V_align_length(in,0) = inVD; //niVD in matlab
  _assigns.insVD_V_align_length(in, 1)= assignment_params.V_align_length;
                                                    
  _assigns.Dallele_D_align_length(in,0) = assignment_params.d_g;_assigns.Dallele_D_align_length(in, 1)=  assignment_params.D_align_length;
  //count<<"vdj palindrome 56"<<endl;                                                                          
  _assigns.VdelV(in,0) = assignment_params.v_g; _assigns.VdelV(in,1)=assignment_params.ndV ;
  _assigns.JdelJ(in,0) = assignment_params.j_g;_assigns.JdelJ(in,1)=assignment_params.ndJ;
  _assigns.VinsVD(in,0) =assignment_params.v_g; _assigns.VinsVD(in,1)= inVD;//niVd in matlab
  _assigns.DinsVD(in,0) = assignment_params.d_g; _assigns.DinsVD(in,1)=inVD;//niVD in matlab
  _assigns.DinsDJ(in,0) = assignment_params.d_g; _assigns.DinsDJ(in,1)= inDJ;//niDJ in matlab
  _assigns.JinsDJ(in,0) = assignment_params.j_g; _assigns.JinsDJ(in, 1)= inDJ;//niDJ in matlab
  _assigns.VdelJ(in,0) = assignment_params.v_g; _assigns.VdelJ(in,1)= assignment_params.ndJ;
  _assigns.JdelV(in,0) = assignment_params.j_g; _assigns.JdelJ(in,1)=assignment_params.ndV;
                                                    
  _assigns.DdelV(in,0) = assignment_params.d_g; _assigns.DdelV(in, 1)= assignment_params.ndV;
  _assigns.DdelJ(in,0) = assignment_params.d_g; _assigns.DdelJ(in,1)= assignment_params.ndJ;
  _assigns.VinsDJ(in,0) = assignment_params.v_g;_assigns.VinsDJ(in,1)= inDJ;//niDJ in matlab
  _assigns.JinsVD(in,0) = assignment_params.j_g;_assigns.JinsVD(in,1)= inVD;//niVD in matlab
  
  //assigns.pVinsVD(in,:) = [1 + npV, 1 + niVD];
  //assigns.pVinsDJ(in,:) = [1 + npV, 1 + niDJ];
  //assigns.pVdelJ(in,:) = [1 + npV, 1 + ndJ];
  //assigns.VpV(in,:) = [v_g, 1 + npV];
  //assigns.JpV(in,:) = [j_g, 1 + npV];
  /*assigns.DpV(in,:) = [d_g, 1 + npV];
    
                                                    assigns.pJinsVD(in,:) = [1 + npJ, 1 + niVD];
                                                    assigns.pJinsDJ(in,:) = [1 + npJ, 1 + niDJ];
                                                    assigns.pJdelV(in,:) = [1 + npJ, 1 + ndV];
                                                    assigns.VpJ(in,:) = [v_g, 1 + npJ];
                                                    assigns.JpJ(in,:) = [j_g, 1 + npJ];
                                                    assigns.DpJ(in,:) = [d_g, 1 + npJ];
                                                    assigns.pVpJ(in,:)  = [1 + npV, 1+ npJ];
                                                    assigns.VcutV(in,:) = [ v_g, 1 - model.min_V_cut + ncutV];
                                                    assigns.JcutJ(in,:) = [ j_g, 1 - model.min_J_cut + ncutJ];
                                                    assigns.DcutV(in,:) = [ d_g, 1 - model.min_V_cut + ncutV];
                                                    assigns.DcutJ(in,:) = [ d_g, 1 - model.min_J_cut + ncutJ];
                                                    assigns.VcutJ(in,:) = [ v_g, 1 - model.min_J_cut + ncutJ];
                                                    assigns.JcutV(in,:) = [ j_g, 1 - model.min_V_cut + ncutV];
                                                    assigns.insVDcutV(in,:) = [ 1 + niVD, 1 - model.min_V_cut + ncutV];
                                                    assigns.insDJcutV(in,:) = [ 1 + niDJ, 1 - model.min_V_cut + ncutV];
                                                    assigns.insVDcutJ(in,:) = [ 1 + niVD, 1 - model.min_J_cut + ncutJ];
                                                    assigns.insDJcutJ(in,:) = [ 1 + niDJ, 1 - model.min_J_cut + ncutJ];
                                                    
                                                    assigns.cutDlcutDr_given_D(in,:) = [ 1 - model.min_D_cut + ncutDl, 1 - model.min_D_cut + ncutDr, d_g];
  */                                                
  //% All possible pairwise joint distributions, just for tracking.
  /*                                                
                                                    assigns.pDldelDl(in,:)=[ 1 + npDl , 1 + ndDl];
                                                    assigns.pDrdelDr(in,:)=[ 1 + npDr , 1 + ndDr];
                                                    assigns.delVdelDl(in,:) = [1 + ndV, 1 + ndDl];
                                                    assigns.delVdelDr(in,:) = [1 + ndV, 1 + ndDr];
                                                    
                                                    assigns.delJdelDl(in,:) = [1 + ndJ, 1 + ndDl];
                                                    assigns.delJdelDr(in,:) = [1 + ndJ, 1 + ndDr];
                                                    assigns.delDlinsVD(in,:) = [1 + ndDl, 1 + niVD];
                                                    assigns.delDlinsDJ(in,:) = [1 + ndDl, 1 + niDJ];
                                                    assigns.delDldelDr(in,:) = [1 + ndDl, 1 + ndDr];
                                                    assigns.delDrinsVD(in,:) = [1 + ndDr, 1 + niVD];
                                                    assigns.delDrinsDJ(in,:) = [1 + ndDr, 1 + niDJ];
                                                    assigns.DdelDl(in,:) = [d_g, 1 + ndDl];
                                                    assigns.DdelDr(in,:) = [d_g, 1 + ndDr ];

                                                    assigns.VdelDl(in,:) = [v_g, 1 + ndDl];
                                                    assigns.VdelDr(in,:) = [v_g, 1 + ndDr];
                                                    assigns.JdelDl(in,:) = [j_g, 1 + ndDl];
                                                    assigns.JdelDr(in,:) = [j_g, 1 + ndDr];
                                                    assigns.pVdelDl(in,:) = [1 + npV, 1 + ndDl];
                                                    assigns.pVdelDr(in,:) = [1 + npV, 1 + ndDr];
                                                    assigns.pJdelDl(in,:) = [1 + npJ, 1 + ndDl];
                                                    assigns.pJdelDr(in,:) = [1 + npJ, 1 + ndDr];
                                                    assigns.pDlinsVD(in,:) = [1 + npDl, 1+niVD];
                                                    assigns.pDlinsDJ(in,:) = [1 + npDl, 1+niDJ];
                                                    assigns.pDldelV(in,:) = [1 + npDl, 1+ndV];
                                                    assigns.pDldelJ(in,:)  = [1 + npDl, 1+ndJ];
                                                    assigns.pDldelDr(in,:) = [1 + npDl, 1+ndDr];
                                                    assigns.VpDl(in,:) = [v_g, 1+ npDl];
                                                    assigns.JpDl(in,:) = [j_g, 1+ npDl];
                                                    assigns.DpDl(in,:) = [d_g, 1+ npDl];
                                                    
                                                    assigns.pDrinsVD(in,:) = [1 + npDr, 1 + niVD];
                                                    assigns.pDrinsDJ(in,:) = [1 + npDr, 1 + niDJ];
                                                    assigns.pDrdelV(in,:)=  [1 + npDr, 1 + ndV];
                                                    assigns.pDrdelJ(in,:)=  [1 + npDr, 1 + ndJ];
                                                    assigns.pDrdelDl(in,:) = [1 + npDr, 1 + ndDl];
                                                    assigns.VpDr(in,:) = [v_g, 1 + npDr];
                                                    assigns.JpDr(in,:) = [j_g, 1 + npDr];
                                                    assigns.DpDr(in,:) = [d_g, 1 + npDr];
                                                    
                                                    assigns.pVpDl(in,:)  = [1 + npV, 1+ npDl];
                                                    assigns.pVpDr(in,:)  = [1 + npV, 1+ npDr];
                                                    assigns.pDlpDr(in,:) = [1 + npDl, 1+ npDr];
                                                    assigns.pDlpJ(in,:) = [1 + npDl, 1+ npJ];
                                                    assigns.pDrpJ(in,:) = [1 + npDr, 1+ npJ];
                                                    
                                                    assigns.DcutDl(in,:) = [ d_g, 1 - model.min_D_cut + ncutDl];
                                                    assigns.DcutDr(in,:) = [ d_g, 1 - model.min_D_cut + ncutDr];
                                                    assigns.VcutDl(in,:) = [ v_g, 1 - model.min_D_cut + ncutDl];
                                                    assigns.VcutDr(in,:) = [ v_g, 1 - model.min_D_cut + ncutDr];
                                                    assigns.JcutDl(in,:) = [ j_g, 1 - model.min_D_cut + ncutDl];
                                                    assigns.JcutDr(in,:) = [ j_g, 1 - model.min_D_cut + ncutDr];
                                                    assigns.insVDcutDl(in,:) = [1 + niVD,  1 - model.min_D_cut + ncutDl];
                                                    assigns.insDJcutDl(in,:) = [1 + niDJ,  1 - model.min_D_cut + ncutDl];
                                                    assigns.insVDcutDr(in,:) = [1 + niVD,  1 - model.min_D_cut + ncutDr];
                                                    assigns.insDJcutDr(in,:) = [1 + niDJ,  1 - model.min_D_cut + ncutDr];
  */                                                
  //_assigns.error_vs_position(in) = error_vs_position;
  //assigns.coverage(in,:) = coverage;
  
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
  in++;
  assignment_params.in=in;
  //_assigns.n_assignments=in;
  
  
  //% If maximum number of valid assignments is reached, stop.
  //% If we are not smoothing, stop after 1 assignment
  //% Usually, we don't reach maximum number of assignments because we break out
  //% of all the loops due to the various thresholds below.
  if (in > _model.max_assignments || (~_do_smoothing && in>1))
    {
      _assigns.n_assignments = in-1;
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



