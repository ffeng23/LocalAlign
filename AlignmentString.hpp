#ifndef ALIGNMENTSTRING_HPP
#define ALIGNMENTSTRING_HPP
#include <string>

using namespace std;

class AlignmentString
{
public:
  AlignmentString();
  AlignmentString(const string& _pattern, const int& p_start, const int& p_end,
		  const string& _subject, const int& c_start, const int& c_end,
		  const string& _pattern_wg, const string& _subject_wg, const double& _score);
  
  //true to set the sub string with gap (true) or without gap(fasle)
  void SetPattern(const string& _pattern, const bool& gapFlag=false);
  void SetSubject(const string& _subject, const bool& gapFlag=false);
  
  string GetPattern(const bool& gapFlag=false);
  string GetSubject(const bool& gapFlag=false);
  
  void SetPatternIndex(const int& _start, const int& _end);
  void SetSubjectIndex(const int& _start, const int& _end);

  unsigned int GetPatternIndexStart();
  unsigned int GetPatternIndexEnd();
  unsigned int GetSubjectIndexStart();
  unsigned int GetSubjectIndexEnd();

  //this is the method 
  //unsigned int FoundAlignedIndexOnOppositeStrand(const byte& _strand, const unsigned int& _index);

  string toString();
  void SetScore(const double& _s);
  double GetScore();
private:
  string c_pattern;//the original aligned substring 
  string c_subject;//the original aligned substring
  
  //starting and ending index of the aligned substring in the full string
  //starting at zero
  int c_pattern_start;
  int c_pattern_end;
  int c_subject_start;  
  int c_subject_end;

  //the substring aligned with gaps
  string c_pattern_w_gap;
  string c_subject_w_gap;

  double c_score;

};
#endif
