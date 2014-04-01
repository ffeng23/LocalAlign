#ifndef LOCALALIGNMENT_HPP
#define LOCALALIGNMENT_HPP
#include "pairwiseAlignment.hpp"
#include <vector>
using namespace std;

//this entry class is used to define each local alignment entry/path
class Path
{
public:
  Path();
  Path( unsigned int _OptimalIndex[2], unsigned int _startIndex [2], const double& _optimalValue );
  void SetOptimalIndex(unsigned int _OptimalIndex[2]);
  void SetStartIndex(unsigned int _startIndex [2]);
  void SetOptimalValue(const double& _optimalValue);
  int* GetOptimalIndex();
  int* GetStartIndex();
  double GetOptimalValue();
private:
  unsigned int c_optimalIndex[2];
  unsigned int c_startIndex[2];
  double c_optimalValue;
};



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

  vector<Path*> c_path_vec;

};


#endif
