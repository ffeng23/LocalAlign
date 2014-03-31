#ifndef LOCALALIGNMENT_HPP
#define LOCALALIGNMENT_HPP
#include "pairwiseAlignment.hpp"

class LocalAlignment: public PairwiseAlignment
{
public:
  LocalAlignment(SequenceString* _pattern, SequenceString* _subject, 
		 ScoreMatrix* _m=&nuc44, const double& _gopen=-8, 
		 const double& _gextension=-5, const double& _scale=1,const int& _numOfAlignments=1);
  
  virtual ~LocalAlignment();

  const double* GetScoreArr();
  const AlignmentString* GetAlignmentArr();
  //void alignLM();
protected:
  virtual void align();
  virtual void traceBack();

  //**********************************
  //the following are the ones used to return and keep track of all non-intersect alignments
  //this is the array holding numOfAlignments requested
  //the one defined by base class only holds the optimal one
  int c_numOfAlignments;
  AlignmentString* c_alignmentArr;
  double* c_scoreArr;
  //in this class we define the alignmentstring and score string and
  //then delete it upon destruction. so the outside caller need to take care(copy)
  //if they need to use the score or alignment after the alignment scope expires.

  
  

};
#endif
