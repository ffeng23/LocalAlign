
#include <iostream>
#include <fstream>
#include "GenomicD.hpp"

using namespace std;

GenomicD::GenomicD(/*const string& _fastaFile, const string& _infoFile*/): Genomic_Segment()
{
  //this->ReadGenomicSegments(_fastaFile, _infoFile);
}

/*bool GenomicV::ReadGenomicSegments(const string& _fastFileName, const string& _infoFileName)
{
  //this is reading genomic V, so we need two files, one fasta, 
  if(_infoFileName.compare("")==0)//we really need the second file input
    {
      cout<<"ERROR: in reading genomic V segments, two files are necessary. Only one file supplied.\n";
      return false;
    }
  //now we need to check the existence of the file
  

  return true;
  }*/

GenomicD::~GenomicD()
{
  //empty;
}
