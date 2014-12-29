#ifndef SEQUENCEHANDLERCOMMON_HPP
#define SEQUENCEHANDLERCOMMON_HPP
#include "SequenceString.hpp"
#include "score.hpp"


SequenceString ReverseComplement(SequenceString& seq);

//compare two strings character by character, return # of chars that are different. if two strings are of different size, 
//the ones longer are also counted as different chars
//for example: "abc" vs. "acc" return 1;
//"abc" vs "acb" return 2
//"abc" vs "ab" reurn 1
//"abc" vs "cab" return 3
unsigned int CompareStrings(const string& str1, const string& str2);

SequenceString Reverse(SequenceString& seq);

#endif