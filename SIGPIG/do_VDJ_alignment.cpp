#include "do_VDJ_alignment.hpp"

//to determine v_end, we need to get right most position on the seq
//align_position+align_length
static unsigned determine_v_end(const Alignment_Object& _v_align)
{
  unsigned v_end=0;
  //go hrough to find the right most postion
  for(unsigned i=0;i<_v_align.numOfAligned;i++)
    {
      if(v_end<(_v_align.align_position[i][0]+_v_align.align_length[i]))
	{
	  v_end=_v_align.align_position[i][0]+_v_align.align_length[i];
	}
    }
  return v_end;
}
//go throught 
static unsigned determine_j_start(const Alignment_Object& _j_align)
{
  unsigned j_start=_j_align.align_position[0][0];
  //go hrough to find the right most postion
  for(unsigned i=0;i<_j_align.numOfAligned;i++)
    {
      if(j_start<(_j_align.align_position[i][0]))
	{
	  j_start=_j_align.align_position[i][0];
	}
    }
  
  return j_start;
}


/*do vdj alignment for each individual input sequence
 *
 */
static bool do_VDJ_alignment_single
(const SequenceString& _seq, const GenomicV* _genVs, const unsigned& _numOfVs,
 const GenomicD* _genDs, const unsigned& _numOfDs,
 const GenomicJ* _genJs, const unsigned& _numOfJs, 
 const double& _error_cost, const ScoreMatrix* _sm,
 const unsigned& _max_D_align,
 /*output*/Alignment_Object& _v_align, Alignment_D& _d_align, Alignment_Object& _j_align)
{
  unsigned v_end=0; //by default(when v is failed), we the start of the sequence
  unsigned j_start=_seq.GetLength()-1; //by default (when match_j is failed), 
       //we use the end of the sequence.
  
//do match J first
  bool seq_j_ok=match_J
	(_seq, _genJs, _numOfJs, 
	 AlignmentSettings::J_minimum_alignment_length, 
	 AlignmentSettings::J_maximum_deletion, 
	 AlignmentSettings::negative_excess_deletions_max, 
	 AlignmentSettings::J_allowed_errors, _error_cost, 
	 _j_align);
  if(!seq_j_ok)
    {
      cout<<"====>match J failed"<<endl;
    }
  else
    {
      j_start=determine_j_start(_j_align);
    }

  //do match V second
  bool seq_v_ok=match_V
	(_seq, _genVs, _numOfVs, 
	 AlignmentSettings::V_minimum_alignment_length, 
	 AlignmentSettings::V_maximum_deletion, 
	 AlignmentSettings::negative_excess_deletions_max, 
	 AlignmentSettings::V_allowed_errors, _error_cost, _v_align);
  if(!seq_v_ok)
    {
      cout<<"====>match v failed"<<endl;
    }
  else
    {
      v_end=determine_v_end(_v_align);
    }
  
  //do match d last
  bool seq_d_ok=match_D
	(_seq, _genDs,  _numOfDs, 
	 v_end, j_start, AlignmentSettings::flank_length, _sm, 
	 AlignmentSettings::D_maximum_deletion, 
	 AlignmentSettings::negative_excess_deletions_max,
	 _max_D_align, _d_align);
  if(!seq_d_ok)
    {
      cout<<"====>match D failed"<<endl;
    }
  return seq_j_ok&&seq_v_ok&&seq_d_ok;
}
/*return alignment objects have been intialized in the caller. the size is
 *_numOfSeq elements
 * this is the function to be called by each working thread.
 */
unsigned do_VDJ_alignment
( /*vector of const sequenstring*/ vector<SequenceString>::const_iterator _it_seq,
  const unsigned& _numOfSeq, 
  const GenomicV* _genVs, const unsigned& _numOfVs,
  const GenomicD* _genDs, const unsigned& _numOfDs,
  const GenomicJ* _genJs, const unsigned& _numOfJs, 
  const double& _error_cost, const ScoreMatrix* _sm,
  const unsigned& _max_D_align,
  /*output*/Alignment_Object* _v_align, Alignment_D* _d_align, 
  Alignment_Object* _j_align
		       )
{
  unsigned numOfAligned=0;
  bool align_ok;
  for(unsigned i=0;i<_numOfSeq;i++)
    {
      //do we need to check for end????
      //calling
      align_ok=do_VDJ_alignment_single(*_it_seq, _genVs, _numOfVs, _genDs, _numOfDs,
				       _genJs, _numOfJs,
			     _error_cost, _sm, _max_D_align, _v_align[numOfAligned],
			     _d_align[numOfAligned], _j_align[numOfAligned]);
      //check for the failure
      if(align_ok)
	{
	  numOfAligned++;
	}
      else
	{
	  _v_align[numOfAligned].ResetData();
	  _d_align[numOfAligned].ResetData();
	  _j_align[numOfAligned].ResetData();
	}
      _it_seq++;
    }

  return numOfAligned;
}
