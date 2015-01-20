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
  unsigned Get_n_alleles() const;
  unsigned Get_Allele() const;
  string Get_Gene() const;//return the name of the gene.

  void Set_Seq(const SequenceString& _seq);
  void Set_GeneIndex(const unsigned int& _gene_index);
  void Set_n_alleles(const unsigned int& _n_alleles);
  void Set_Allele(const unsigned int& _allele);
  void Set_Gene(const string& _gene); 
  /*
  //for v, we need two files, but D, J, we only need the first file containing segments
  virtual bool ReadGenomicSegments(const string& _fastaFileName, const string& _infoFileName="")=0;
  */
//private
private: 
  SequenceString c_seq;
  string c_gene;//this is the gene name without allele number
  unsigned c_gene_index;//index of the genes in the gene array, not sure where the array is now (?)
  //but probably will use this array later.
  unsigned c_n_alleles;  //this is the total number of alleles for this gene segment
  unsigned c_allele; //which allele in terms of number is this current one
  
};
#endif
