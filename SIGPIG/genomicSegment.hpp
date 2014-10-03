#ifndef GENOMIC_SEGMENT_HPP
#define GENOMIC_SEGMENT_HPP

//this is the class defining the template for genomic gene segaments from the library.
//Will be instantiated into 3 different segments V, D, J

#include <string>

#include "../SequenceString.hpp"

using namespace std;

//abstract class
class Genomic_Segment
{
public:
  Genomic_Segment();
  virtual ~Genomic_Segment()=0;

  //accessors
  SequenceString Get_Seq() const;
  string Get_Sequence() const;
  string Get_Name() const;
  unsigned Get_GeneIndex() const;
  unsigned Get_TotalNumberOfAlleles() const;
  unsigned Get_AlleleNum() const;
  

  /*
  //for v, we need two files, but D, J, we only need the first file containing segments
  virtual bool ReadGenomicSegments(const string& _fastaFileName, const string& _infoFileName="")=0;
  */
//private
private: 
  SequenceString c_seq;
  unsigned c_gene_index;
  unsigned c_n_alleles;  //this is the total number of alleles for this gene segment
  unsigned c_allele; //which allele in terms of number is this current one
  
};
#endif
