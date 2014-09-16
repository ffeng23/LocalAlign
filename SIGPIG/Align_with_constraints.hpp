#ifndef ALIGN_WITH_CONSTRAINTS_HPP
#define ALIGN_WITH_CONSTRAINTS_HPP

#include "../score.hpp"
#include "../SequenceString.hpp"
#include "../AlignmentString.hpp"
#include "../AffineGapModel.hpp"
#include "../MarkovChainGapModel_454.hpp"
#include "../TracebackTable.hpp"


class Align_with_constraints:public PairwiseAlignment
{
 public:
 Align_with_constraints(SequenceString* _pattern, SequenceString* _subject, 
		 ScoreMatrix* _m=&nuc44, const double& _gopen=-8, 
		 const double& _gextension=-5, const double& _scale=1);

 virtual ~Align_with_constraints();

protected:
  virtual void align();
  virtual void traceBack();
public: //this following public members is in order to be compatible to the matlab code
  unsigned Align_with_constraints_fast_left(SequenceString const* _pattern, SequenceString const* _subject, const unsigned int& maximum_errors
					    , const double& error_cost);
  unsigned Align_with_constraints_fixed_left(SequenceString const* _pattern, SequenceString const* _subject, const unsigned int& maximum_errors
					    , const double& error_cost);
  
  unsigned 
  

};
#endif
