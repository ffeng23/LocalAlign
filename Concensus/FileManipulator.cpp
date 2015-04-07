#include <iostream>
#include "FileManipulator.hpp"
#include <stdlib.h>

using namespace std;
//the caller need to declare, but don't initialize, since we don't know how many,
//but the caller has to delete/clean up memory
void GetFileNames(const char* _path, /*output*/std::string** fileNames, unsigned& _numOfFiles)
{
  DIR* d = opendir(_path);
  _numOfFiles=0;
  if (d == NULL) 
    {
      cout<<"ERROR: can not open the directory, quit"<<endl;
      cerr<<"ERROR: can not open the directory, quit"<<endl;
      exit(-1);
    }
  for(struct dirent *de = NULL; (de = readdir(d)) != NULL; )
    {
      string fn(de->d_name);
      
      if(fn.find(".fasta")==std::string::npos)
	{
	  cout<<"WARNING::find one irrelevant file "<< de->d_name <<endl;
	  continue;
	}
      //else we are good
      _numOfFiles++;
      
    }
  //do the 
  (*fileNames)=new string[_numOfFiles];
  rewinddir(d);
  _numOfFiles=0;
  for(struct dirent *de = NULL; (de = readdir(d)) != NULL; )
    {
      string fn(de->d_name);
      
      if(fn.find(".fasta")==std::string::npos)
	{
	  cout<<"WARNING::find one irrelevant file "<< de->d_name <<endl;
	  continue;
	}
      //else we are good
      
      (*fileNames)[_numOfFiles]=de->d_name;
      _numOfFiles++;
    }
  //done
  closedir(d);
}
