#include "genomicSegment.hpp"

Genomic_Segment::Genomic_Segment(): c_seq(),c_gene_index(0),c_n_alleles(0), c_allele(0)
{
  //empty
}
Genomic_Segment::~Genomic_Segment()
{
  //empty
}

SequenceString Genomic_Segment::Get_Seq() const
{
  return c_seq;
}

string Genomic_Segment::Get_Sequence() const
{
  return c_seq.GetSequence();
}

string Genomic_Segment::Get_Name() const
{
  return c_seq.GetName();
}
unsigned Genomic_Segment::Get_TotalNumberOfAlleles() const
{
  return c_n_alleles;
}
unsigned Genomic_Segment::Get_AlleleNum() const
{
  return c_allele;
}

//bool Genomic_Segments::ReadGenomicSegments(const string& _fileName)
