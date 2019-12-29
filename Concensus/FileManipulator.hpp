#ifndef FILE_MANIPULATOR_HPP
#define FILE_MANIPULATOR_HPP

#include <dirent.h>
#include <string>

//the caller need to declare, but don't initialize, since we don't know how many,
//but the caller has to delete/clean up memory
void GetFileNames(const char* _path, /*output*/std::string** fileNames, unsigned& _numOfFiles);


#endif
