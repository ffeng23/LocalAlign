#include "SequenceHandlerCommon.hpp"

using namespace std;


unsigned int CompareStrings(const string& str1, const string& str2)
{
  unsigned int ret=0;
  unsigned int len1, len2, len;
  len1=str1.length();len2=str2.length();
  len=len1;

  
  if(len>len2)
    {
      len=len2;
      ret=len1-len2;
    }
  else
    {
      ret=len2-len1;
    }
  for(unsigned int i=0;i<len;i++)
    {
      
      if(str1.at(i)!=str2.at(i))
	{
	  ret++;
	}
    }
  
  return ret;
}


SequenceString ReverseComplement(SequenceString& seq)
{
  string tempStr=seq.GetSequence();
  SequenceString temp(seq.GetName(), "");
  string tempStrReturn("");
  //cout<<"temStr (seq get sequence):"<<tempStr<<endl;
  //cout<<"length:"<<tempStr.length()<<endl;
  for(unsigned int i=tempStr.length();i>0;i--)
    {
      //cout<<"\tloop i="<<i<<endl;
      switch(tempStr.at(i-1))
	{
	case 'A':
	case 'a':
	  tempStrReturn.push_back('T');
	  break;
	case 'T':
	case 't':
	  tempStrReturn.push_back('A');
	  break;
	case 'U':
	case 'u':
	  tempStrReturn.push_back('A');
	  break;
	case 'G':
	case 'g':
	  tempStrReturn.push_back('C');
	  break;
	case 'C':
	case 'c':
	  tempStrReturn.push_back('G');
	  break;
	case 'Y':
	case 'y':
	  tempStrReturn.push_back('R');
	  break;
	case 'R':
	case 'r':
	  tempStrReturn.push_back('Y');
	  break;
	case 'S':
	case 's':
	  tempStrReturn.push_back('S');
	  break;
	case 'W':
	case 'w':
	  tempStrReturn.push_back('W');
	  break;
	case 'K':
	case 'k':
	  tempStrReturn.push_back('M');
	  break;
	case 'M':
	case 'm':
	  tempStrReturn.push_back('K');
	  break;
	case 'B':
	case 'b':
	  tempStrReturn.push_back('V');
	  break;
	case 'D':
	case 'd':
	  tempStrReturn.push_back('H');
	  break;
	case 'H':
	case 'h':
	  tempStrReturn.push_back('D');
	  break;
	case 'V':
	case 'v':
	  tempStrReturn.push_back('B');
	  break;
	case 'N':
	case 'n':
	default:
	  tempStrReturn.push_back('N');
	  break;
	}
      
    }

  temp.SetSequence(tempStrReturn);

  return temp;
}
