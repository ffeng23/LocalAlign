#ifndef FILE_MANIPULATOR_HPP
#define FILE_MANIPULATOR_HPP

#include <dirent.h>
//the caller need to declare, but don't initialize, since we don't know how many,
//but the caller has to delete/clean up memory
void GetFileNames(const char* _path, /*output*/string* fileNames, unsigned& _numOfFiles);


#endif
