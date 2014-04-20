#ifndef AFFINEGAPMODEL_HPP
#define AFFINEGAPMODEL_HPP

#include "pairwiseAlignment.hpp"

class AffineGapModel:public GapModel
{

public:
  AffineGapModel(const double& _gopen, const double& _gextension);

  virtual ~AffineGapModel();
  virtual double GapValue(TracebackTable* _tbTable, const unsigned int& _rowIndex, const unsigned int& _colIndex, const bool& _rowGap,
			  double& _MaxGapValue, unsigned int& _MaxGapIndex) const;

  double GetGapOpenValue() const;
  double GetGapExtensionValue() const;
  void SetGapOpenValue(const double& _gopen);
  void SetGapExtensionValue(const double& _gextension);

protected:
  double c_gopen;
  double c_gextension;
  
}
#endif 
