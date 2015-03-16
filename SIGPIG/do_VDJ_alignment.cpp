#include <pthread.h>
#include "do_VDJ_alignment.hpp"
extern pthread_mutex_t progressMutex;
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
  cout<<"------------>start j alignment"<<endl;
  bool seq_j_ok=match_J
	(_seq, _genJs, _numOfJs, 
	 AlignmentSettings::J_minimum_alignment_length, 
	 AlignmentSettings::J_maximum_deletion, 
	 AlignmentSettings::negative_excess_deletions_max, 
	 AlignmentSettings::J_allowed_errors, _error_cost, 
	 _j_align);
  //cout<<"done!!!"<<endl;
  if(!seq_j_ok)
    {
      cout<<"====>match J failed"<<endl;
    }
  else
    {
      j_start=determine_j_start(_j_align);
    }
  cout<<"------------>start v alignment"<<endl;
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
  cout<<"---------->start d alignment"<<endl;
  cout<<"j_start:"<<j_start<<";v_end:"<<v_end<<endl;
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
  Alignment_Object* _j_align, bool* _vdj_align_ok
		       )
{
  unsigned numOfAligned=0;
  bool align_ok;
  pthread_mutex_lock(&progressMutex);
  cout<<"<<)))))))_numOfSeq:"<<_numOfSeq<<endl;
  pthread_mutex_unlock(&progressMutex);
  for(unsigned i=0;i<_numOfSeq;i++)
    {
      //do we need to check for end????
      //calling
      pthread_mutex_lock(&progressMutex);
      cout<<"numOfAligned:"<<numOfAligned<<endl;
      cout<<"before doing the alignment:"<<endl;
      cout<<"\t_j_align[]:"<<_j_align[numOfAligned].align_length<<endl;
      pthread_mutex_unlock(&progressMutex);
      align_ok=do_VDJ_alignment_single(*_it_seq, _genVs, _numOfVs, _genDs, _numOfDs,
				       _genJs, _numOfJs,
			     _error_cost, _sm, _max_D_align, _v_align[i],
			     _d_align[i], _j_align[i]);
      pthread_mutex_lock(&progressMutex);
      cout<<"after doing the alignment:"<<endl;
      cout<<"\t_j_align[]"<<_j_align[i].align_length<<endl;
      
      cout<<"====doing work i:"<<i<<"/"<<_numOfSeq<<"!!!"<<endl;
      cout<<"\talign_ok:"<<align_ok<<endl;
      pthread_mutex_unlock(&progressMutex);
      //check for the failure
      if(align_ok)
	{
	  cout<<" A good one with numOfAligned:"<<numOfAligned<<endl;
	  numOfAligned++;
	  _vdj_align_ok[i]=align_ok;
	}
      else
	{//in this new version, we want to save everything in case we want to 
	  //check for debugging, but keep another indicator array whethter this
	  //one is a good. or not.
	  _vdj_align_ok[i]=align_ok;
	  /*
	  pthread_mutex_lock(&progressMutex);
	  cout<<"\t*****resetting !!!"<<endl;
	  pthread_mutex_unlock(&progressMutex);
	  _v_align[numOfAligned].ResetData();
	  pthread_mutex_lock(&progressMutex);
	  cout<<"\t*****resetting 2!!!"<<endl;
	  pthread_mutex_unlock(&progressMutex);
	  _d_align[numOfAligned].ResetData();
	  pthread_mutex_lock(&progressMutex);
	  cout<<"\t*****resetting 3!!!"<<endl;
	  cout<<_j_align[numOfAligned].toString()<<endl;
	  pthread_mutex_unlock(&progressMutex);
	  _j_align[numOfAligned].ResetData();

	  pthread_mutex_lock(&progressMutex);
	  cout<<"\t*****resetting DONE !!!"<<endl;
	  pthread_mutex_unlock(&progressMutex);
	  */
	}
      _it_seq++;
      pthread_mutex_lock(&progressMutex);
      cout<<"\tend of loop i:"<<i<<" !!!"<<endl;
      pthread_mutex_unlock(&progressMutex);
    }
  pthread_mutex_lock(&progressMutex);
  cout<<"\t====done in pthread claling functin"<<endl;
  //cout<<"\talign_ok:"<<align_ok<<endl;
  pthread_mutex_unlock(&progressMutex);
  return numOfAligned;
}

//p_thread function wrapper.
void * do_VDJ_align_pthread(void * param_thread)
{
  //now decoding the param input
  param_alignment_pthread* p_a_p=(param_alignment_pthread*) param_thread;
  
  //now calling it
  pthread_mutex_lock(&progressMutex);
  cout<<"==>thread::["<<p_a_p->p_thread_id<<"]:"<<endl;
   pthread_mutex_unlock(&progressMutex);

  *p_a_p->p_numOfAligned =do_VDJ_alignment(p_a_p->p_it, p_a_p->p_numOfSeq, 
		   p_a_p->p_genVs, p_a_p->p_numOfVs,
		   p_a_p->p_genDs, p_a_p->p_numOfDs,
		   p_a_p->p_genJs, p_a_p->p_numOfJs,
		   p_a_p->p_error_cost, p_a_p->p_sm,
		   p_a_p->p_max_D_align, p_a_p->p_v_align,
					   p_a_p->p_d_align, p_a_p->p_j_align, 
					   p_a_p->p_vdj_align_ok);
  
  pthread_mutex_lock(&progressMutex);
  cout<<"====end of thread!!!"<<endl;
  pthread_mutex_unlock(&progressMutex);
  pthread_exit(NULL);//return numOfAligned as a pointer
}

void PackUpAlignmentParameterForThread
  (
   /*params for do_VDJ_alignment function*/ 
   vector<SequenceString>::const_iterator _it, 
   const unsigned& _numOfSeq, 
   const GenomicV* _genVs, const unsigned& _numOfVs,
   const GenomicD* _genDs, const unsigned& _numOfDs,
   const GenomicJ* _genJs, const unsigned& _numOfJs, 
   const double& _error_cost, const ScoreMatrix* _sm,
   const unsigned& _max_D_align,
   Alignment_Object* _v_align, Alignment_D* _d_align, 
   Alignment_Object* _j_align,
   unsigned* _numOfAligned,
   /*pointer to param to use*/
   param_alignment_pthread* p_a_p,
   const unsigned _thread_id,
   bool* _vdj_align_ok
   )
{
  p_a_p->p_it=_it;
  p_a_p->p_numOfSeq=_numOfSeq;
  p_a_p->p_genVs=_genVs; p_a_p->p_numOfVs=_numOfVs;
  p_a_p->p_genDs=_genDs; p_a_p->p_numOfDs=_numOfDs;
  p_a_p->p_genJs=_genJs; p_a_p->p_numOfJs=_numOfJs;
  p_a_p->p_error_cost=_error_cost; p_a_p->p_sm=_sm;
  p_a_p->p_max_D_align=_max_D_align; 
  p_a_p->p_v_align=_v_align;
  p_a_p->p_d_align=_d_align; 
  p_a_p->p_j_align=_j_align;
  p_a_p->p_numOfAligned= _numOfAligned;
  p_a_p->p_thread_id=_thread_id;
  p_a_p->p_vdj_align_ok=_vdj_align_ok;
  //done;

}

