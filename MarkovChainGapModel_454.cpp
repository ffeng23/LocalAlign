#include "454_MarkovChainGapModel.hpp"

//this model is for 454 sequencing reading, to model the umbiguity of stretch of identical nucleotides run. it is using affine model too

454_MarkovChainGapModel::454_MarkovChainGapModel(const double& _gopen, const double& _gextension, 
						 const string& _patternStr, const string& _subjectStr)
:AffineGapModel(_gopen, _gextension), c_patternString(_patternStr), c_subjectString(_subjectStr)
{
  //empty
}

454_MarkovChainGapModel::~454_MarkovChainGapModel()
{
  //empty
}
double 454_MarkovChainGapModel::GapValue(const TracebackTable* _tbTable, const unsigned int& _patternIndex, const unsigned int& _subjectIndex, const bool& _patternGap,
			  double& _MaxGapValue, unsigned int& _MaxGapIndex)
{
  //to start, there are only two possible cases: 1) the extended from the MaxGapIndex;2)a new open gap
  //for both cases, we need to check whether this is a stretch of identical nucleotides and then decide how we 
  //are going to apply cost.
  //first take care of extended from the MaxGapIndex
  //check to see whether we need to use the Markov Chain Model
  char curr_char;
  bool regular_cost_flag=true;
  double extendedGapValue;
  if(_patternGap)//doing row gap now
    {
      curr_char=c_subjectString[_subjectIndex-1];
      for(unsigned int i=1;_patternIndex-i<=_MaxGapIndex&&_patternIndex-i!=0&&_subjectIndex-i!=0;i++)
	{
	  if(c_subjectString[_subjectIndex-i-1]!=curr_char||c_patternString[_patternIndex-i-1]!=curr_char)
	    {
	      regular_cost_flag=false;
	      break;
	    }
	  
	}
      if(regular_cost_flag) //regular
	{
	  extendedGapValue=_MaxGapValue+c_gextension;
	}
      else //Markov chain
	{
	  //need to check how long the run is
	  unsigned int runLen=0;
	  //if we are here, _patternIndex-_MaxGapIndex won't be zero, and _subjectIndex-(_patternIndex-_MaxGapIndex)won't either
	  /*while(c_subjectString[_subjectIndex-(_patternIndex-_MaxGapIndex)-runLen-1]!=curr_char||c_patternString[_MaxGapIndex-runLen-1]!=curr_char)
	    {
	      
	      runLen++;
	      if(_subjectIndex-(_patternIndex-_MaxGapIndex)-runLen-1<0&&_MaxGapIndex-runLen-1<0)
		{
		  break;
		}
		}*/
	  //try to trace back using trace back table in order to figure out the length of the run of the identical nucleotide
	  while(_tbTable->GetLink(_MaxGapIndex-runLen, _subjectIndex-(_patternIndex-_MaxGapIndex)-runLen)==UPLEFT)
	    {//here we don't have to check whether the index will get negative, since the zero ones won't pass the check of while clause.
	      if(c_subjectString[_subjectIndex-(_patternIndex-_MaxGapIndex)-runLen-1]==curr_char||c_patternString[_MaxGapIndex-runLen-1]==curr_char)
		runLen++;
	      else
		break;
	    }
	  //get the gap len
	  unsigned int gapLen=_patternIndex-_MaxGapIndex;

	  //get the cost
	  extendedGapValue=GetMarkovChainGapValue(c_gextension, runLen, _gapLen);	  
	}
    }
  else //doing subject gap
    {
      curr_char=c_patternString[_patternIndex-1];
      for(unsigned int i=1;_subjectIndex-i<=_MaxGapIndex&&_subjectIndex-i!=0&&_patternIndex-i!=0;i++)
	{
	  if(c_patternString[_patternIndex-i-1]!=curr_char||c_subjectString[_subjectIndex-i-1]!=curr_char)
	    {
	      regular_cost_flag=false;
	      break;
	    }
	  
	}
      if(regular_cost_flag) //regular
	{
	  extendedGapValue=_MaxGapValue+c_gextension;
	}
      else //Markov chain
	{
	  //need to check how long the run is
	  unsigned int runLen=0;
	  //if we are here, _patternIndex-_MaxGapIndex won't be zero, and _subjectIndex-(_patternIndex-_MaxGapIndex)won't either
	  /*while(c_patternString[_patternIndex-(_subjectIndex-_MaxGapIndex)-runLen-1]!=curr_char||c_subjectString[_MaxGapIndex-runLen-1]!=curr_char)
	    {
	      
	      runLen++;
	      if(_patternIndex-(_subjectIndex-_MaxGapIndex)-runLen-1<0&&_MaxGapIndex-runLen-1<0)
		{
		  break;
		}
	    }
	  */
	  while(_tbTable->GetLink(_patternIndex-(_subjectIndex-_MaxGapIndex)-runLen,_MaxGapIndex-runLen )==UPLEFT)
	    {//here we don't have to check whether the index will get negative, since the zero ones won't pass the check of while clause.
	      if(c_patternString[_patternIndex-(_subjectIndex-_MaxGapIndex)-runLen-1]==curr_char||c_subjectString[_MaxGapIndex-runLen-1]==curr_char)
		runLen++;
	      else
		break;
	    }
	  //get the gap len
	  unsigned int gapLen=_subjectIndex-_MaxGapIndex;

	  //get the cost
	  extendedGapValue=GetMarkovChainGapValue(c_gextension, runLen, _gapLen);	  
	}
    }//patterngap or subjectgap

  //now we want to solve the problem of figure out the open new gap value
  runLen=0;
  if(_patternGap)//pattern gap
    {
      while(_tbTable->GetLink(_patternIndex-runLen,_subjectIndex-runLen-1 )==UPLEFT)
	    {//here we don't have to check whether the index will get negative, since the zero ones won't pass the check of while clause.
	      if(c_patternString[_patternIndex-runLen]==curr_char||c_subjectString[_subjectIndex-runLen-1]==curr_char)
		runLen++;
	      else
		break;
	    }
    }
  else
    {
      while(_tbTable->GetLink(_patternIndex-runLen-1,_subjectIndex-runLen )==UPLEFT)
	    {//here we don't have to check whether the index will get negative, since the zero ones won't pass the check of while clause.
	      if(c_patternString[_patternIndex-runLen-1]==curr_char||c_subjectString[_subjectIndex-runLen]==curr_char)
		runLen++;
	      else
		break;
	    }
    }
  double newOpenGapValue;
  if(runLen==0) //regular one
    {
      newOpenGapValue=c_gopen+c_extension;
    }
  else
    {
      newOpenGapValue =GetMarkovChainGapValue(c_gopen, c_gextension, runLen, 1);
    }
  
  
}

double 454_MarkovChainGapModel::GetMarkovChainGapValue(const double& _regularGapOpen, const double& _regularGapExtension, const unsigned int& _runLen, const unsigned int& _gapLen) const
{
  return (_regularGapOpen+_regularGapExtension)*(1/(1+20*exp(-20.0*(((double)_gapLen)/(_gapLen+_runLen)-0.5))));
} 
