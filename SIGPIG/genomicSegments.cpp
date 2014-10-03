#include "genomicSegments.hpp"

Genomic_Segments::Genomic_Segments(): c_seq(),c_gene_index(0),c_n_alleles(0), c_allele(0)
{
  //empty
}
Genomic_Segments::~Genomic_Segments()
{
  //empty
}

SequenceString Genomic_Segments::Get_Seq() const
{
  return c_seq;
}

string Genomic_Segments::Get_Sequence() const
{
  return c_seq.GetSequence();
}

string Genomic_Segments::Get_Name() const
{
  return c_seq.GetName();
}
unsigned Genomic_Segments::Get_TotalNumberOfAlleles() const
{
  return c_n_alleles;
}
unsigned Genomic_Segments::Get_AlleleNum() const
{
  return c_allele;
}

//bool Genomic_Segments::ReadGenomicSegments(const string& _fileName)
