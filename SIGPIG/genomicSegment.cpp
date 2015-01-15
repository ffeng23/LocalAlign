#include "genomicSegment.hpp"

Genomic_Segment::Genomic_Segment(): c_seq(),c_gene(),c_gene_index(0),c_n_alleles(0), c_allele(0)
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
unsigned Genomic_Segment::Get_n_alleles() const
{
  return c_n_alleles;
}
unsigned Genomic_Segment::Get_Allele() const
{
  return c_allele;
}
unsigned Genomic_Segment::Get_GeneIndex() const
{
  return c_gene_index;
}

//bool Genomic_Segments::ReadGenomicSegments(const string& _fileName)


void Genomic_Segment::Set_Seq(const SequenceString& _seq)
{
  c_seq=_seq;
}
void Genomic_Segment::Set_GeneIndex(const unsigned int& _gene_index)
{
  c_gene_index=_gene_index;
}
void Genomic_Segment::Set_n_alleles(const unsigned int& _n_alleles)
{
  c_n_alleles=_n_alleles;
}
void Genomic_Segment::Set_Allele(const unsigned int& _allele)
{
  c_allele=_allele;
}

string Genomic_Segment::Get_Gene() const
{
  return c_gene;
}

void Genomic_Segment::Set_Gene(const string& _gene)
{
  c_gene=_gene;
}
