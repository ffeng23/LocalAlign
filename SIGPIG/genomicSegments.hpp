#ifndef GENOMIC_SEGMENTS_HPP
#define GENOMIC_SEGMENTS_HPP

//this is the class defining the template for genomic gene segaments from the library.
//Will be instantiated into 3 different segments V, D, J

#include <string>
#include <sstream>
#include <vector>

#include "../SequenceString.hpp"
#include "genomicSegment.hpp"
#include "GenomicJ.hpp"
#include "GenomicV.hpp"
#include "GenomicD.hpp"



using namespace std;

/*
//abstract class
class Genomic_Segments
{
public:
  Genomic_Segments();
  virtual ~Genomic_Segments()=0;

  //accessors
  SequenceString Get_Seq() const;
  string Get_Sequence() const;
  string Get_Name() const;
  unsigned Get_GeneIndex() const;
  unsigned Get_TotalNumberOfAlleles() const;
  unsigned Get_AlleleNum() const;

  //for v, we need two files, but D, J, we only need the first file containing segments
  virtual bool ReadGenomicSegments(const string& _fastaFileName, const string& _infoFileName="")=0;
  //private
private: 
  SequenceString c_seq;
  unsigned c_gene_index;
  unsigned c_n_alleles;  //this is the total number of alleles for this gene segment
  unsigned c_allele; //which allele in terms of number is this current one
  
  };*/
/******here in this file, we are not doing any class, instead we simply define some functions to read the file and populate the genomic segments objects*/
//this is the function to read all the similar things from
//input:
//   file name, fasta file name for gene segment sequence
//   info filename, we need this only for geneV read, not others

//Output: the pointer to the array of the segments
//the caller doesn't initialize the array, because the caller doesn't know 
//how many are there, so the callee will initialize the array,
//but the caller has to clean up the memeory afterwards
unsigned ReadGenomicV(const string& _fastaFileName, GenomicV** _gseg);

//see the definition above

unsigned ReadGenomicD(const string& _fastaFileName, GenomicD** _gseg);

unsigned ReadGenomicJ(const string& _fastaFileName, GenomicJ** _gseg);

//here we assume the name is arranged with a fix format
//for example >Jxxxx|IGHJ1*01|Homo Sapiens|F|......
string ParseName(const string& _seqName);

string ParseGeneName(const string& _seqName);

string ParseSequenceName(const string& _seqName);

//only flip the sequence of the object, but do not change the name
//return a new object and leave the input object no change
SequenceString FlipSequenceString(const SequenceString& _ss);

vector<string> DetermineOutputFileNames(const string& _outFileNameBase, const unsigned& _NPerFile, const unsigned& _totalNSeq);

bool GenomicVCompare_bySequenceName(const GenomicV& ss1,const GenomicV& ss2);

#endif
