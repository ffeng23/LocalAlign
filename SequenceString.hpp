#ifndef SEQUENCESTRING_HPP
#define SEQUENCESTRING_HPP

#include <string>
using namespace std;

class SequenceString
{
public:
  SequenceString();

  SequenceString(const string& _name, const string& _seq);

  void SetName(const string& _name);
  void SetSequence(const string& _seq);

  const string GetName();
  const string GetSequence();
  const unsigned GetLength();

  string toString(bool _fasta=false);

private:
  string c_name;
  string c_seq;
  //unsigned int c_len;

};
#endif
