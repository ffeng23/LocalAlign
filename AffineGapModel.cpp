#include "AffineGapModel.hpp"

AffineGapModel::AffineGapModel(const double& _gopen, const double& _gextension):
  c_gopen(_gopen), c_gextension(_gextension)
{

  //empty
}

~AffineGapModel()
{
  //empty
}

//input
//_tbTable--the tracebackTable so far.
//_patternIndex, _colIndex -- this is the "coordinate" of current entry,
//_rowGap -- this is the flag to indicate whether we are doing row gap (true) or col gap (false).
//_prevEntryValue -- this is the score of the previous entry in the same row (for row gap) or in the same col (for col gap)
//_MaxGapValue -- the running Max Gap value
//_MaxGapValue-- the index of the Max Gap starting point
//
//Output
//_MaxGapValue
//_MaxGapIndex
virtual double AffineGapModel::GapValue(TracebackTable* _tbTable, const unsigned int& _rowIndex, const unsigned int& _colIndex, const bool& _rowGap,
			const double& _prevEntryValue,
			  double& _MaxGapValue, unsigned int& _MaxGapIndex) const;
{
  //are we doing it by row of by or by col
  //for affine gap value, we have to know this, we only need to 
  //get the maxiGapvalue and compare with the newly opened one to get 
  //the best ones
  double tempMaxGapValue=_MaxGapValue;
  

  //newly opened gap
  double newGapOpenValue=_prevEntryValue+c_gopen+c_gextension;
  double extendedGapValue=_MaxGapValue+c_gextension;

  if(newGapOpenValue>=extendedGapValue)
    {
      _MaxGapValue=newGapOpenValue;
      if(_rowGap)//we are doing row gaps
	{
	  //_colIndex--;
	  _MaxGapIndex=_colIndex-1;
	  
	}
      else
	{
	  //_rowIndex--;
	  _MaxGapIndex=_rowIndex-1;
	}
      
    }
  else
    {
      _MaxGapValue=extendedGapValue;
    }
  //_tbTable[_rowIndex+_colIndex*
  return _MaxGapValue;
}
double AffineGapModel::GetGapOpenValue() const
{
  return c_gopen;
}
double AffineGapModel::GetGapExtensionValue() const
{
  return c_gextension;
}
void AffineGapModel::SetGapOpenValue(const double& _gopen)
{
  c_gopen=_gopen;
}
void AffineGapModel::SetGapExtensionValue(const double& _gextension)
{
  c_gextension=_gextension;
}




