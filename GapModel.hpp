#ifndef GAPMODEL_HPP
#define GAPMODEL_HPP

#include "pairwiseAlignment.hpp"

//this is the abstract class base class for gap models

class GapModel
{
public:
  GapModel();
  virtual ~GapModel()=0;
  
  virtual double GapValue(TracebackTable* _tbTable, unsigned int _rowIndex, unsigned int _colIndex, bool _rowGap,
			  double& _MaxGapValue, unsigned int& _MaxGapIndex)=0;
  
};


#endif
