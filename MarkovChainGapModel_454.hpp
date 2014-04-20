#ifndef MARKOVCHAINGAPMODEL_HPP
#define MARKOVCHAINGAPMODEL_HPP

#include <cmath>

class 454_MarkovChainGapModel:public AffineGapModel
{
public:
  454_MarkovChainGapModel(const double& _gopen, const double& _gextension,const string& _patternStr, const string& _subjectStr);
  virtual ~454_MarkovChainGapModel();
  virtual double GapValue(TracebackTable* _tbTable, const unsigned int _patternIndex, const unsigned int _subjectIndex, const bool& _patternGap,
			  double& _MaxGapValue, unsigned int& _MaxGapIndex) const;
protected:
  string c_patternString;
  string c_subjectString;

  double GetMarkovChainGapValue(const double& _regularGapExtension, const unsigned int& _runLen, const unsigned int& _gapLen) const; 
};
#endif
