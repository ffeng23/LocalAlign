#ifndef DO_VDJ_ALIGNMENT_HPP
#define DO_VDJ_ALIGNMENT_HPP

#include <vector>

#include "MatrixFunctions.hpp"
#include "Alignment.hpp"
#include "Alignment_V.hpp"
#include "Alignment_D.hpp"
#include "genomicSegments.hpp"
#include "../SequenceString.hpp"


/*return alignment objects have been intialized in the caller. the size is
 *_numOfSeq elements
 */
unsigned do_VDJ_alignment
( /*vector of const sequenstring*/ vector<SequenceString>::const_iterator _it, 
  const unsigned& _numOfSeq, 
  const GenomicV* _genVs, const unsigned& _numOfVs,
  const GenomicD* _genDs, const unsigned& _numOfDs,
  const GenomicJ* _genJs, const unsigned& _numOfJs, 
  const double& _error_cost, const ScoreMatrix* _sm,
  const unsigned& _max_D_align,
  /*output*/ Alignment_Object* _v_align, Alignment_D* _d_align, 
  Alignment_Object* _j_align, bool* _vdj_align_ok
  );


//================following code to define things for multi-threading 
//this is the struct wrapper for the parameter the alignment
struct param_alignment_pthread
{
  vector<SequenceString>::const_iterator p_it;
  unsigned p_numOfSeq;
  const GenomicV* p_genVs; unsigned p_numOfVs;
  const GenomicD* p_genDs; unsigned p_numOfDs;
  const GenomicJ* p_genJs; unsigned p_numOfJs;
  double p_error_cost; 
  const ScoreMatrix* p_sm;
  unsigned p_max_D_align;
  Alignment_Object* p_v_align;
  Alignment_D* p_d_align;
  Alignment_Object* p_j_align;
  //aboved are the parameters to do_vdj_alignment function
  unsigned* p_numOfAligned;/*return parameter from vdj alignment function*/
  
  unsigned p_thread_id;
  bool* p_vdj_align_ok;
};

void * do_VDJ_align_pthread(void * param_thread);

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
   );


#endif

