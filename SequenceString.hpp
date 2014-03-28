#ifndef SEQUENCESTRING_HPP
#define SEQUENCESTRING_HPP

#include <string>
using namespace std;

class SequenceString
{
public:
  SequenceString();

  SequenceString(string _name, string _seq);

  void SetName(string _name);
  void SetSequence(string _seq);

  string GetName();
  string GetSequence();

private:
  string c_name;
  string c_seq;

};
#endif
