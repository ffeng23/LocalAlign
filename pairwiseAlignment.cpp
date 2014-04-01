#include "pairwiseAlignment.hpp"
#include <iostream>
using namespace std;


//pairwiseAlignment();
PairwiseAlignment::PairwiseAlignment(SequenceString* _pattern, SequenceString* _subject, 
		    ScoreMatrix* _m, const double& _gopen, 
				     const double& _gextension, const double& _scale):
  c_pattern( _pattern), c_subject (_subject), c_sm( _m), c_gapOpen(_gopen), 
  c_gapExtension(_gextension), c_scale( _scale), 
  c_traceback_table(NULL)
{
  //need to allocate the 
  //this->align();
  c_optimalIndex[0]=0;
  c_optimalIndex[1]=0;
}
  
PairwiseAlignment::~PairwiseAlignment()
{
  c_pattern=NULL;
  c_subject=NULL;
  c_sm=NULL;
  
  if(c_traceback_table)
    delete[] c_traceback_table;
  
}

double PairwiseAlignment::GetScore()
{
  return c_score;
}
AlignmentString PairwiseAlignment::GetAlignment()
{
  return c_alignment;
}

void PairwiseAlignment::traceBack()
{
  //emptyhere need to fill
  cout<<"doing traceback in pairwiseAlignment"<<endl;
}



TracebackTableEntry::TracebackTableEntry():c_link(ZERO), c_numOfIndels(0)
{
  //empty
}

TracebackTableEntry::TracebackTableEntry(const LinkBack& _link, const unsigned int& _numOfIndels):
  c_link(_link), c_numOfIndels(_numOfIndels)
{
  //empty
} 
  
LinkBack TracebackTableEntry::GetLinks()
{
  return c_link;
}
unsigned int TracebackTableEntry::GetNumOfIndels()
{
  return c_numOfIndels;
}
  
void TracebackTableEntry::SetLinks(const LinkBack& _link)
{
  c_link=_link;
}
void TracebackTableEntry::SetNumOfIndels(const unsigned int _numOfIndels)
{
  c_numOfIndels=_numOfIndels;
}
