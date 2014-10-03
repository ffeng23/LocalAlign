
#include "SequenceString.hpp"
#include<iostream>
#include <string>
#include <sstream>

using namespace std;
SequenceString::SequenceString():c_name(""),c_seq("")
{
  //empty
}

SequenceString::SequenceString(const string& _name, const string& _seq) : c_name(_name),c_seq(_seq)
{
  //empty
}
void SequenceString::SetName(const string& _name)
{
  c_name=_name;
}
void SequenceString::SetSequence(const string& _seq)
{
  c_seq=_seq;
  //c_len=_seq.length()
}

const string SequenceString::GetName() const
{
  return c_name;
}

const string SequenceString::GetSequence() const
{
  return c_seq;
}

const unsigned int SequenceString::GetLength() const
{
  return c_seq.length();
}

string SequenceString::toString(bool _fasta) const
{
  ostringstream ss("");
  if(!_fasta)
    {
       ss<<c_seq.length() <<"-character long SequenceString instance\n"<<c_name<<":"<<c_seq<<"\n";
    }
  else
    {
      ss<< ">"<<c_name<<"\n"<<c_seq<<"\n";
    }
  return ss.str();
}

bool SequenceString::operator < (const SequenceString& other) const
{
	int ret=c_seq.compare(other.c_seq);
	//cout<<"compare in the function:"<<ret<<endl;
	if(ret<=0)
		return true;
	else
		return false;
}
