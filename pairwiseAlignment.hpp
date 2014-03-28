#ifndef PAIRWISEALIGNMENT_HPP
#define PAIRWISEALIGNMENT_HPP

#include "score.hpp"
#include "SequenceString.hpp"

//this file is used to define a interface/abstract class for pairwise alignment classes
//LocalAlign
//overlapAlign
//globalAlign will be implemented in the future

//Abstract class
class pairwiseAlign
{
  //public
public:
  //pairwiseAlignment();
  pairwiseAlignment(const SequenceString* _pattern, const SequenceSting* _sujbect, 
		    const ScoreMatrix* _m=&nuc44, const double& _gopen=-8, 
		    const double& _gextension=-5, const double& _scale=1);
  
  ~pairwiseAlignment();

  virtual double Score()=0;
  virtual AlignmentString GetAlignment()=0;

protected:
  virtual align()=0;
  ScoreMatrix* c_sm;
  double c_gapOpen;
  double c_gapExtension;
  double c_scale;
  
  const SequenceString* c_pattern;//this is follwing R style pairwiseAlignment
  const SequenceString* c_subject;//this is following R style pairwiseAlignment
  
  AlignmentString* c_alignment;
  double score;
}

#endif
