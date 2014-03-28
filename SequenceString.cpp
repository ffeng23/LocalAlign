
#include "SequenceString.hpp"

SequenceString::SequenceString():c_name(""),c_seq("")
{
  //empty
}

SequenceString::SequenceString(string _name, string _seq) : c_name(_name),c_seq(_seq)
{
  //empty
}
void SequenceString::SetName(string _name)
{
  c_name=_name;
}
void SequenceString::SetSequence(string _seq)
{
  c_seq=_seq;
}

string SequenceString::GetName()
{
  return c_name;
}

string SequenceString::GetSequence()
{
  return c_seq;
}
