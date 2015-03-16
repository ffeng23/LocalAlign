#ifndef DATAIO_HPP
#define DATAIO_HPP

#include "GenomicJ.hpp"
#include "GenomicD.hpp"
#include "GenomicV.hpp"
#include "genomicSegments.hpp"
#include "LoadData.hpp"
#include "AlignmentSettings.hpp"
#include "Alignment.hpp"
#include "Alignment_V.hpp"
#include "Alignment_D.hpp"


//this module define the functions to read/write data files for alignment
//also take care of serialization deserialization of alignments

//====================>
//calling the ReadGenomicVDJ functions to do the reading, here we just 
//write code to put them together for future usage

//NOTE:: the caller don't initialize the _genVDJ array, since they don't
//know how many segments is going to read in.
//BUT THEY NEED TO clear up the memory
//finally return a boolean to indicating a successful or failed action
bool DoGenomicTemplateReading
(const string& _genomicVFileName,
 const string& _genomicDFileName,
 const string& _genomicJFileName,
 /*output*/
 GenomicV** _genV, GenomicD** _genD, GenomicJ** _genJ,
 unsigned& _totalNumV, unsigned& _totalNumD, unsigned& _totalJ
 );

//===========>
//for data sequence loading, no wrapper here, but simply 
//calling LoadData() function in genomicSegments module.
//bool DoSequenceDataLoading(

//=========>serialization of alignment datas
bool DoSerialization
(bool const* const* _vdj_align_ok_arrays, Alignment_Object const* const* _v_align_arrays,
 Alignment_D const* const* _d_align_arrays, Alignment_Object const* const* _j_align_arrays,
 vector<SequenceString>::const_iterator* _it_arrays , 
 const unsigned* _numOfAlignedSequenceForThread, 
 const unsigned* _numOfInputSequencesForThread, const unsigned& _numOfThread, 
 const string& _fileName
 );
//NOTE about the structure of the serialization
//1)first we store the version "ver 0.01"
//2) the total number of aligment sets, an unsigned, 
//3) each set contains seq, v_align, 
//           d align and v align.
//4) an unsigned used to indicating the number of set so far.
//       this number could be used to indicating the integraty of the serialization file
//       when doing deserialization.
//


//=========>deserialization of alignment datas
//reading the alignment from the disk, 
//the caller DON'T initialize the output arrays, but only declare it.
//also the caller need to delete/clean up the memory.
bool DoDeserialization
(const string& _fileName, unsigned& _numOfAligned,
 Alignment_Object** _v_align, Alignment_D** _d_align, Alignment_Object** _j_align,
 SequenceString** _seq
 );
#endif
