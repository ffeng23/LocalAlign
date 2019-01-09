#ifndef TRACEBACKTABLE_HPP
#define TRACEBACKTABLE_HPP

#include "Accessory/SequenceString.hpp"

//this is struct to record the link back for tracing
//zero is for local align to indicate this is a zero node for terminating a path
enum LinkBack
  {
    UP, LEFT, UPLEFT, ZERO
  };


class TracebackTableEntry
{
public:
  TracebackTableEntry();
  TracebackTableEntry(const LinkBack& _link, const unsigned int& _numOfIndels=0); 
  
  LinkBack GetLinks();
  unsigned int GetNumOfIndels();
  
  void SetLinks(const LinkBack& _link);
  void SetNumOfIndels(const unsigned int _numOfIndels);

private:
  LinkBack c_link;//this is pointer pointing to the one leading to this, could be left, up, upleft, zero
  unsigned int c_numOfIndels;//this is only works for left or up(indels),using to indicate how indels leads to this,not necessarily only 1.
};

class TracebackTable
{
public:
  //typeOfAlignment, 0, global;1, overlap; 2, local.
  TracebackTable(SequenceString* _pattern, SequenceString* _subject, const short& _typeOfAlignment=0);
  ~TracebackTable();
  void SetTableEntry(const unsigned int& _patternIndex, const unsigned int& _subjectIndex, const LinkBack& _link, const unsigned int& _numOfIndels);

  LinkBack GetLink(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const;
  unsigned int GetNumOfIndels(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const ;
private:
  TracebackTableEntry* c_tbt;
  unsigned int c_lenOfPattern;
  unsigned int c_lenOfSubject;
};

#endif
