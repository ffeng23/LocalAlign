#ifndef GLOBALALIGNMENT_HPP
#define GLOBALALIGNMENT_HPP
#include "pairwiseAlignment.hpp"
#include <vector>
using namespace std;


class GlobalAlignment: public PairwiseAlignment
{
public:
  GlobalAlignment(SequenceString* _pattern, SequenceString* _subject, 
		 ScoreMatrix* _m=&nuc44, const double& _gopen=-8, 
		 const double& _gextension=-5, const double& _scale=1);
  
  virtual ~GlobalAlignment();

  //const double* GetScoreArr();
  //const AlignmentString* GetAlignmentArr();
  //void alignLM();
protected:
  virtual void align();
  virtual void traceBack();
};


#endif
