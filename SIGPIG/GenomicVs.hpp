#ifndef GENOMICVS_HPP
#define GENOMICVS_HPP

#include "genomicSegments.hpp"

//inherited from Genomic_Segments class
class GenomicVs : public Genomic_Segments
{
public:
  GenomicVs(const string& _fastaFile, const string& _infoFile);
  ~GenomicVs();
  

  //for this one, we have to have two files, but not one file name
  virtual bool ReadGenomicSegments(const string& _fastFileName, const string& _infoFileName="");

};
#endif
