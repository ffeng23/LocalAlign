#ifndef PAIRWISEALIGNMENT_HPP
#define PAIRWISEALIGNMENT_HPP

#include "score.hpp"
#include "SequenceString.hpp"
#include "AlignmentString.hpp"
#include "AffineGapModel.hpp"
#include "MarkovChainGapModel_454.hpp"
#include "TracebackTable.hpp"

//this file is used to define a interface/abstract class for pairwise alignment classes
//LocalAlign
//overlapAlign
//globalAlign will be implemented in the future


//Abstract class
class PairwiseAlignment
{
  //public
public:
  //pairwiseAlignment();
  PairwiseAlignment(SequenceString* _pattern, SequenceString* _subject, 
		    ScoreMatrix* _m=&nuc44, const double& _gopen=-8, 
		    const double& _gextension=-5, const double& _scale=1,const short& _typeOfGapModel=1);
  
  virtual ~PairwiseAlignment()=0;

  double GetScore();
  AlignmentString GetAlignment();

protected:
  virtual void align()=0;
  virtual void traceBack();

  SequenceString* c_pattern;//this is follwing R style pairwiseAlignment
  SequenceString* c_subject;//this is following R style pairwiseAlignment
  ScoreMatrix* c_sm;
  double c_gapOpen;
  double c_gapExtension;
  double c_scale;
  
  AlignmentString c_alignment;
  double c_score;
  unsigned int c_optimalIndex[2];
  
  //double* c_dp_table;
  TracebackTable* c_traceback_table;
  GapModel* c_gm;
};

#endif
