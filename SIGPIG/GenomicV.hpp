#ifndef GENOMICV_HPP
#define GENOMICV_HPP

#include "genomicSegment.hpp"

//inherited from Genomic_Segment class

//this is exactly like the base class, just for consistence purpuse, we make it a inherited class
class GenomicV : public Genomic_Segment
{
public:
  GenomicV(/*const string& _fastaFile, const string& _infoFile*/);
  ~GenomicV();
  
  /*
  //for this one, we have to have two files, but not one file name
  virtual bool ReadGenomicSegments(const string& _fastFileName, const string& _infoFileName="");
  */
};
#endif
