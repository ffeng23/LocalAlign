
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
  //start preparing for parameters
  VDJ_model_assignments_settings assignment_params;

  assignment_params.log_probability_threshold_factor=log(_probability_threshold_factor);
  assignment_params.log_probability_hopeless_threshold_factor=assignment_params.log_probability_threshold_factor+log(1E-6);
  assignment_params.in=0;
  assignment_params.skips=0;
 
  assignment_params.max_J_depth=AlignmentSettings::max_J_length;
  assignment_params.max_V_depth=_seq.GetLength()+ AlignmentSettings::max_V_length;
  
  assignment_params.J_max_error=0;
  assignment_params.D_max_error=0;

  assignment_params.assume_palindrome_negative_deletions=true;

  assignment_params.log_max_model_p_nt_DJ=matrix_log(max(_model.RnucleotideDJ_per_nucleotideDJ_3prime,2));
  assignment_params.log_max_model_p_nt_VD=matrix_log(max(_model.RnucleotideVD_per_nucleotideVD_5prime,2));
  assignment_params.log_max_model_p_nt=max(assignment_params.log_max_model_p_nt_DJ);
  double temp=max(assignment_params.log_max_model_p_nt_VD);
  if(assignment_params.log_max_model_p_nt<temp)
    {
      assignment_params.log_max_model_p_nt=temp;
    }

  assignment_params.log_RnucleotideDJ_per_nucleotideDJ_3prime=matrix_log(_model.RnucleotideDJ_per_nucleotideDJ_3prime.m2vec());
  assignment_params.log_RnucleotideVD_per_nucleotideVD_5prime=matrix_log(_model.RnucleotideVD_per_nucleotideVD_5prime.m2vec());

  if(_model.Rerror_per_sequenced_nucleotide>0)
    {
      assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3=
	log(_model.Rerror_per_sequenced_nucleotide/3);
    }
  else
    {
      assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3=-1000.0;
    }
  
  assignment_params.L_err=_seq.GetLength();
  assignment_params.log_proba_Rerror_normalization=-1.0*assignment_params.L_err*log(1+_model.Rerror_per_sequenced_nucleotide);
  assignment_params.log_max_model_pcutV_given_gene=matrix_log(max(_model.PcutV_given_V,1));
  assignment_params.log_max_model_pcutJ_given_gene=matrix_log(max(_model.PcutJ_given_J,1));
  Matrix<double> temp_m=max(_model.PcutDlcutDr_given_D,2);
  assignment_params.log_max_model_pcutD_given_gene=matrix_log(max(temp_m,1));
  
  assignment_params.log_max_model_pinsVD=log(max(_model.PinsVD));
  assignment_params.log_max_model_pinsDJ=log(max(_model.PinsDJ));
  assignment_params.log_max_model_pins=assignment_params.log_max_model_pinsVD+assignment_params.log_max_model_pinsDJ;
  
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
  assignment_params.log_highest_probability=-1000.0;
  assignment_params.best_D_align_length=0;
  //assignment_params.max_nerrorsv_1=0;
  //====>??????call vdj assign

  return true;

}
			   

bool assign_VDJ_alleles
(VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 const double& _probability_threshold_factor, const bool& _no_error,
 const bool& _ignore_deep_error, const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
 /*output, input*/VDJ_model_assignments_settings& assignment_params,
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns )
{
  unsigned dim_size[4]={0,0,0,0};
  unsigned deletable_errors=0;
  unsigned max_nerrorsv_1=0;
  unsigned min_nerrorsv, max_nerrorsv;

  //double log_max_pcutV_loop_v;
  for(unsigned v=0;v<_V.numOfAligned;v++)
    {
      assignment_params.v=v;
      //check for validity
      if(assignment_params.skips>assignment_params.max_skips)
	{
	  _assigns.n_assignments=assignment_params.in;
	  _assigns.skips=assignment_params.skips;
	  return true;
	}
      assignment_params.v_a=_V.alleles_all[v];
      assignment_params.v_g=_genV[assignment_params.v_a].Get_GeneIndex();

      //for high_error_region, we don't use in my code,
      //but keep it for now anyway
      if(_ignore_deep_error)
	assignment_params.high_error_region=_model.high_error_region;
      else
	assignment_params.high_error_region=_model.high_error_region;

      if(v==0)//for first round, set some starting point
	{
	  dim_size[0]=_V.n_errors[v];
	  assignment_params.v_err_pos.initialize(1, dim_size, _V.error_positions[v]);
	  deletable_errors=sum_all_bool(assignment_params.v_err_pos > (_V.align_length[v]-_model.max_excess_V_deletions));
	  min_nerrorsv=_V.n_errors[v]-deletable_errors;
	  //here it is different from Matlab code,
	  //since we don't have high_error region anyway.
	  //We simply don't account for it
	  max_nerrorsv_1=_V.n_errors[v];//again, this is different from Matlab code.
	  //since we don't have high error region
	  if(_no_error && min_nerrorsv>0)
	    {
	      //we don't allow errors , but do have errors in here, so do next
	      continue;
	    }
	}
      else  //for other cases not v==1
	{
	  //check if this v has too many more eerors than the best V
	  assignment_params.v_err_pos.clear();
	  dim_size[0]=_V.n_errors[v];
	  assignment_params.v_err_pos.initialize(1, dim_size, _V.error_positions[v]);
	  deletable_errors=sum_all_bool(assignment_params.v_err_pos > (_V.align_length[v]-_model.max_excess_V_deletions));
	  min_nerrorsv=_V.n_errors[v]-deletable_errors;
	  //here it is different from Matlab code,
	  //since we don't have high_error region anyway.
	  //We simply don't account for it
	  max_nerrorsv=_V.n_errors[v];//again, this is different from Matlab code.
	  //since we don't have high error region
	  if(_no_error && min_nerrorsv>0)
	    {
	      //we don't allow errors , but do have errors in here, so do next
	      continue;
	    }
	  //check for too many more errors compared with best V
	  double log_min_diff_perrv=(max_nerrorsv-max_nerrorsv_1) * assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
	  if(log_min_diff_perrv<assignment_params.log_probability_threshold_factor)
	    {
	      assignment_params.skips+=1;
	      continue;//why did not count skips in the above no_error position
	    }
	  
	}//end of else
      assignment_params.log_max_pcutV_loop_v=assignment_params.log_max_model_pcutV_given_gene(assignment_params.v_g);

      assignment_params.log_highest_probability_GIVEN_current_V_allele=-1000.0;
      assignment_params.v_break_out=false;
      assignment_params.n_assignments_v_gene=0;

      //========loop over J alleles
      for(unsigned j=0;j<_J.numOfAligned;j++)
	{
	  assignment_params.j=j;
	  if(assignment_params.skips>assignment_params.max_skips)
	    {
	      _assigns.n_assignments=assignment_params.in;
	      _assigns.skips=assignment_params.skips;
	      return true; //=====is true good here
	    }
	  if(_J.align_length[j]<_model.min_J_align_length)
	    continue;
	  assignment_params.j_a=_J.alleles_all[j];
	  assignment_params.j_g=_genJ[assignment_params.j_a].Get_GeneIndex();
	  
	  assignment_params.niVD_DJ0=_J.align_position[j][0] -(_V.align_position[v][0]+ _V.align_length[v]-1)-1;
	  //NOTE:::here the code is different, since the alignment of V not starting from zero. By adding _V.align_position[v][0]

	  assignment_params.log_max_pcutJ_loop_j=assignment_params.log_max_model_pcutJ_given_gene(assignment_params.j_g);

	  assignment_params.log_highest_probability_GIVEN_current_J_allele=-1000.0;
	  assignment_params.n_assignment_j_gene=0;
	  assignment_params.j_break_out=false;

	  //=======loop over D alleles
	  for(unsigned d_i=0;d_i<_numD;d_i++)
	    {
	      if(assignment_params.skips>assignment_params.max_skips)
		{
		  _assigns.n_assignments=assignment_params.in;
		  _assigns.skips=assignment_params.skips;
		  return true; //====is true good here?
		}
	      unsigned d=_D.allele_order[d_i];
	      assignment_params.d=d;
	      assignment_params.d_a=d;
	      assignment_params.d_g=_genD[assignment_params.d_a].Get_GeneIndex();

	      assignment_params.n_assignments_d_gene=0;
	      assignment_params.d_break_out=false;

	      //base probability of gene choices, multiplied by conditional probability for allele choices, given genes
	      //J have no distinugishable alles for original code, but not for project
	      double probabase=_model.PV(assignment_params.v_g);
	      probabase=probabase*_model.PDJ(assignment_params.d_g, assignment_params.j_g);
	      probabase*=_model.PVallele_given_gene(assignment_params.v_g, _genV[assignment_params.v_a].Get_Allele() );
	      probabase*=_model.PDallele_given_gene(assignment_params.d_g, _genD[assignment_params.d_a].Get_Allele());
	      probabase*=_model.PJallele_given_gene( assignment_params.j_g, _genJ[assignment_params.j_a].Get_Allele());
	      //code change*****, different from matlab, since we add PJallele_given_gene. 
	      //since we can distinguish J allele now.
	      //unsigned dim_pos[]={AlignmentSettings::D_maximum_deletion+1,AlignmentSettings::D_maximum_deletion+1,0};
	      assignment_params.p_max_Dl=max_mf2(_D.p_region_max_length_left[d],_D.numOfAligned[d],AlignmentSettings::D_maximum_deletion+1);
	      assignment_params.p_max_Dr=max_mf2(_D.p_region_max_length_left[d],_D.numOfAligned[d],AlignmentSettings::D_maximum_deletion+1);
	      //dim_pos[0]=j;
	      assignment_params.p_max_J=max_mf(_J.p_region_max_length[j],AlignmentSettings::J_maximum_deletion+1);
	      
	      //upper boun on insertion nt bias factor
	      assignment_params.niVD_DJ_min=0;
	      //	      dim_pos[0]=v;
	      unsigned temp_int=assignment_params.niVD_DJ0- _D.align_length[d][0]-max_mf(_V.p_region_max_length[v], AlignmentSettings::V_maximum_deletion+1)-assignment_params.p_max_J - assignment_params.p_max_Dl- assignment_params.p_max_Dr;
	      if(temp_int>0)
		assignment_params.niVD_DJ_min=temp_int;
	      double log_p_max_nt_VD_DJ_d_allele=assignment_params.niVD_DJ_min*assignment_params.log_max_model_p_nt;
	      
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
		
	      assignment_params.log_highest_probability_GIVEN_D_allele=-1000.0;

	      //=======now loop over V deletions
	      //???assign_V_deletions();<==========

	    }//end of d loop ====
	  
	}//end of j loop =====
      
    }//end of for outer for v loop  ========
  return true;
}//end of function of assign VDJ alleles

//return false if something not right, like too many skips and so we
//want the outer caller finish too. just stop doing this alignment.
bool assign_VJ_deletions
(VDJ_cuts_insertion_dinuc_ntbias_model& _model, const SequenceString& _seq,
 const Alignment_Object& _V, const Alignment_D& _D, const Alignment_Object& _J,
 const GenomicV* _genV, const unsigned& _numV,
 const GenomicD* _genD, const unsigned& _numD,
 const GenomicJ* _genJ, const unsigned& _numJ,
 const double& _probability_threshold_factor, const bool& _no_error,
 const bool& _ignore_deep_error, const bool& _do_smoothing,
 const bool& _force_all_alleles, const unsigned& _READ_LENGTH_CORRECTION,
 /*output, input*/VDJ_model_assignments_settings& assignment_params,
 /*output*/VDJ_cuts_insertion_dinuc_ntbias_assigns& _assigns )

{
  //% Loop over V deletions.
  //% ndV1 is the deviation from the no. of deletions implied by the alignment (min_deletions).
  //% Lower bound is nd_start (usually -3).
  //% Upper bound is set by maximum deletions allowed and by minimum match length required.
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
      if (assignment_params.skips > assignment_params.max_skips)
	{
	  _assigns.n_assignments = assignment_params.in;
	  _assigns.skips =  assignment_params.skips;
	  return false;
	}//       end
                
      //% Actual number of deletions
      assignment_params.ndV = _V.min_deletions[assignment_params.v] + assignment_params.ndV1;

      //% This is a violation of a constraint.
      //% Happens for eg. when min_deletions is 0 and I consider 'negative' deletions.
      if( assignment_params.ndV < 0 || assignment_params.ndV > ((signed)_model.max_V_deletions))
	//                    % Go to next iteration of loop.
                    continue;
      //        end

      //% Calculate number of errors in V section, by leaving out errors that are in deleted region for this assignment.
      //% Also leaves out errors in left most 'high_error_region'.
                
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
	  if(assignment_params.npV_max > _V.p_region_max_length[assignment_params.v][ 1+assignment_params.ndV] )
	    {
	      assignment_params.npV_max=_V.p_region_max_length[assignment_params.v][ 1+assignment_params.ndV]; 
	    }
	  assignment_params.npV_potential_max = assignment_params.npV_max;
	}
      else
	{
	  assignment_params.npV_max = 0;
	  assignment_params.npV_potential_max = _model.max_palindrome;
	  if(assignment_params.npV_potential_max> _V.p_region_max_length[assignment_params.v] [1+assignment_params.ndV]);
	  {
	    assignment_params.npV_potential_max=_V.p_region_max_length[assignment_params.v][1+assignment_params.ndV];
	  }
	}                
                
      //Upper bound on insertions nt bias factor
      assignment_params.niVD_DJ_min= assignment_params.niVD_DJ0+assignment_params.ndV1-_D.align_length[assignment_params.d][0] - assignment_params.npV_max - assignment_params.p_max_J-assignment_params.p_max_Dl-assignment_params.p_max_Dr;
      
      if(((signed)assignment_params.niVD_DJ_min) <0)
	{
	  assignment_params.niVD_DJ_min=0;
	}
      double log_p_max_nt_VD_DJ_dV = (assignment_params.niVD_DJ_min)*assignment_params.log_max_model_p_nt;
                
      assignment_params.log_highest_probability_GIVEN_current_V_deletions = -1000;// % Highest probability so far GIVEN current V deletions AND outer loop variables {V, J, D alleles}                
      
      double test_proba = assignment_params.log_probabase + assignment_params.log_perrv + (assignment_params.log_max_pcutV_loop_v+assignment_params.log_max_pcutJ_loop_j+assignment_params.log_max_pcutD_loop_d) + assignment_params.log_max_model_pins + log_p_max_nt_VD_DJ_dV;               
                
      if (test_proba < assignment_params.log_probability_threshold_factor + assignment_params.log_highest_probability)
	{
                    assignment_params.skips  ++;
                    continue;
	}
                
      assignment_params.V_align_length = _V.align_length[assignment_params.v] - assignment_params.ndV1 - assignment_params.READ_LENGTH_CORRECTION;


      //======now start doing j deletions loop for
      //% Loop over J right deletions.
      //% ndJ1 is the deviation from the no. of deletions implied by the alignment (min_deletions).
      //% Lower bound is nd_start (usually -3).
      //% Upper bound is set by maximum deletions allowed and by minimum match length required
      unsigned temp_array[]={_model.max_excess_J_deletions,_J.align_length[assignment_params.j]-_model.min_J_assign_length , _model.max_J_deletions - _J.min_deletions[assignment_params.j]};
      unsigned ndJ_max = min_mf(temp_array, 3);
      for( assignment_params.ndJ1=assignment_params.nd_start;assignment_params.ndJ1<(signed)ndJ_max;assignment_params.ndJ1++)
	{	//		%ndJ1=nd_start
	  if (assignment_params.skips > assignment_params.max_skips)
	    {
	      _assigns.n_assignments = assignment_params.in;
	      _assigns.skips = assignment_params.skips;
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
                                 
	  //niVD_DJ_total_with_ps = niVD_DJ0 + ndV1 + ndJ1;
          //          % Violation of constraint. V and J overlap in this assignment. So move on until this is not true.
	  if (assignment_params.niVD_DJ0+ assignment_params.ndV1+assignment_params.ndJ1 < 0)
	    continue;
                                      
	  //% Calculate number of errors in J section, by leaving out errors that are in deleted region for this assignment.
	  dim_size[0]=_J.n_errors[assignment_params.j];
	  Matrix<unsigned> j_err_pos(1, dim_size, _J.error_positions[assignment_params.j]); //% error positions in J alignment
	  dim_size[0]=3;
	  Matrix<unsigned> j_err_excess_pos(1, dim_size, _J.excess_error_positions[assignment_params.j]);// % error positions in extended J alignment for 'negative' deletions
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
	      assignment_params.nerrorsj = sum_all_bool( j_err_pos >= (_J.align_position[assignment_params.j][0] + assignment_params.ndJ1));
	    }//      end
                    
	  if (_no_error && assignment_params.nerrorsj > 0)
	    continue;
	                      
	  if( assignment_params.nerrorsj > assignment_params.J_max_error || ((double)assignment_params.nerrorsj)/(_J.align_length[assignment_params.j]-assignment_params.ndJ1) > 0.3)
	    {
	      //					disp('in the condition')
	      continue;
	    }
                    
                    
	  assignment_params.log_perrj = assignment_params.nerrorsj*assignment_params.log_Rerror_per_sequenced_nucleotide_divided_by_3;
          
	  //%                     % Upper bound on insertions nt bias factor
	  assignment_params.niVD_DJ_min =  assignment_params.niVD_DJ0 + assignment_params.ndV1 + assignment_params.ndJ1 - _D.align_length[assignment_params.d][0] - _V.p_region_max_length[assignment_params.v][1 + assignment_params.ndV] - _J.p_region_max_length[assignment_params.j][1 + assignment_params.ndJ] - assignment_params.p_max_Dl - assignment_params.p_max_Dr;
	  if(assignment_params.niVD_DJ_min<0)
	    {
	      assignment_params.niVD_DJ_min=0;
	    }
	  double log_p_max_nt_VD_DJ_dJ = (assignment_params.niVD_DJ_min)*assignment_params.log_max_model_p_nt;
          
	  assignment_params.log_highest_probability_GIVEN_current_J_deletions = -1000;                   
	  
	  double test_proba = assignment_params.log_probabase + assignment_params.log_perrv + assignment_params.log_perrj + (assignment_params.log_max_pcutV_loop_v+assignment_params.log_max_pcutJ_loop_j+assignment_params.log_max_pcutD_loop_d) + assignment_params.log_max_model_pins + log_p_max_nt_VD_DJ_dJ;
	            
	  if( test_proba < assignment_params.log_probability_threshold_factor + assignment_params.log_highest_probability)
	    continue;
                   
                    
	  // % Maximum possible half-lengths of palindromes, given current deletions values
	  if (assignment_params.ndJ==0 || ~assignment_params.assume_palindrome_negative_deletions)
	    {
	      assignment_params.npJ_max = _model.max_palindrome;
	      if(assignment_params.npJ_max> _J.p_region_max_length[assignment_params.j][1+assignment_params.ndJ])
		
		assignment_params.npJ_potential_max = assignment_params.npJ_max;
	    }
	  else
	    {
              assignment_params.npJ_max = 0;
	      assignment_params.npJ_potential_max = _model.max_palindrome; 
	      if(assignment_params.npJ_potential_max>_J.p_region_max_length[assignment_params.j][ 1+assignment_params.ndJ])
		assignment_params.npJ_potential_max=_J.p_region_max_length[assignment_params.j][ 1+assignment_params.ndJ];
		
	    }//             end
                    
	  //-------->NOTE: the following part is need to be revised for the future-------
             /*       
	  if (np_start_from_max)
	    {
	      np_step = -1;
	      npV_start = npV_max;
	      npV_end = 0;
	      npJ_start = npJ_max;
	      npJ_end = 0;
	    }
	  else
	    {
	      np_step = 1;
	      npV_start = 0;
	      npV_end = npV_max;
	      npJ_start = 0;
	      npJ_end = npJ_max;
	    }//       end
	     */     
	  assignment_params.J_align_length = _J.align_length[assignment_params.j] - assignment_params.ndJ1;                                        
	  //% Loop over half-length of V palindrome
	  //for now it is just a place holder for the following code
	  //it simply return true and hopefully this will work for now.
	  if(!assign_VJ_palindrome())
	    {
	      return false;
	    }
	}//end of J deletion loop
    }//end of V deletions loop for
  
  return true;
}//end of assign_V_deletions

bool assign_VJ_palindrome()
{
  return true;
}



