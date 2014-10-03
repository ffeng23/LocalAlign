#ifndef GENOMICJ_HPP
#define GENOMICJ_HPP

#include "genomicSegment.hpp"

//inherited from Genomic_Segment class
//this is the individual J segment class
//compared to the base class, there is one more field in this J seg
class GenomicJ : public Genomic_Segment
{
public:
  GenomicJ(/*const string& _fastaFile, const string& _infoFile*/);
  ~GenomicJ();
  
  /*
  //for this one, we have to have two files, but not one file name
  virtual bool ReadGenomicSegments(const string& _fastFileName, const string& _infoFileName="");
  */
private:
  SequenceString c_seq_untill_primer;
};
#endif
