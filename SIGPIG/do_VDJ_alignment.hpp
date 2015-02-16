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
  Alignment_Object* _j_align
		  );




#endif
