#ifndef SEQUENCESTRING_HPP
#define SEQUENCESTRING_HPP

#include <string>
#include <fstream>
using namespace std;

class SequenceString
{
public:
  SequenceString();
  
  ~SequenceString();//destructor empty one

  //no copy constructor since it is not necessary. since there is no pointer or dynamically allocated mem

  SequenceString(const string& _name, const string& _seq);

  void SetName(const string& _name);
  void SetSequence(const string& _seq);

  const string GetName()const;
  const string GetSequence()const;
  const unsigned GetLength()const;

  string toString(bool _fasta=false) const;

  bool operator < (const SequenceString& other) const;

  void Serialize(ofstream& _ofs)const;
  void Deserialize(ifstream& _ifs);
private:
  string c_name;
  string c_seq;
  //unsigned int c_len;

};
#endif
