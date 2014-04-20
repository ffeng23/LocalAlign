#ifndef PAIRWISEALIGNMENT_HPP
#define PAIRWISEALIGNMENT_HPP

#include "score.hpp"
#include "SequenceString.hpp"
#include "AlignmentString.hpp"

//this file is used to define a interface/abstract class for pairwise alignment classes
//LocalAlign
//overlapAlign
//globalAlign will be implemented in the future

//this is struct to record the link back for tracing
//zero is for local align to indicate this is a zero node for terminating a path
enum LinkBack
  {
    UP, LEFT, UPLEFT, ZERO
  };

class TracebackTableEntry
{
public:
  TracebackTableEntry();
  TracebackTableEntry(const LinkBack& _link, const unsigned int& _numOfIndels=0); 
  
  LinkBack GetLinks();
  unsigned int GetNumOfIndels();
  
  void SetLinks(const LinkBack& _link);
  void SetNumOfIndels(const unsigned int _numOfIndels);

private:
  LinkBack c_link;//this is pointer pointing to the one leading to this, could be left, up, upleft, zero
  unsigned int c_numOfIndels;//this is only works for left or up(indels),using to indicate how indels leads to this,not necessarily only 1.
};

class TracebackTable
{
public:
  //typeOfAlignment, 0, global;1, overlap; 2, local.
  TracebackTable(SequenceString* _pattern, SequenceString* _subject, const short& _typeOfAlignment=0);
  ~TracebackTable();
  void SetTableEntry(const unsigned int& _patternIndex, const unsigned int& _subjectIndex, const LinkBack& _link, const unsigned int& _numOfIndels);

  LinkBack GetLink(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const;
  unsigned int GetNumOfIndels(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const ;
private:
  TracebackTableEntry* c_tbt;
  unsigned int c_lenOfPattern;
  unsigned int c_lenOfSubject;
};


//Abstract class
class PairwiseAlignment
{
  //public
public:
  //pairwiseAlignment();
  PairwiseAlignment(SequenceString* _pattern, SequenceString* _subject, 
		    ScoreMatrix* _m=&nuc44, const double& _gopen=-8, 
		    const double& _gextension=-5, const double& _scale=1);
  
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
};

#endif
