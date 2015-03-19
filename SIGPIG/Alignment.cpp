
#include <cstring>
#include <sstream>
#include "../SequenceString.hpp"
#include "genomicSegments.hpp"
#include "Alignment.hpp"
#include "MatrixFunctions.hpp"
#include "../string_ext.hpp"
#include "AlignmentSettings.hpp"

using namespace std;


Alignment_Object::Alignment_Object(const unsigned& numOfGenTemplates):
  numOfAligned(numOfGenTemplates), maximum_deletion(0), align_length(NULL),
  align_position(NULL), min_deletions(NULL), n_errors(NULL), error_positions(NULL),
  p_region_max_length(NULL), excess_error_positions(NULL), alleles_all(NULL),
  alleles_from_distinct_genes(NULL)
{
  //empty constructor
}

//copy constructor, we need a deep copy
Alignment_Object::Alignment_Object(const Alignment_Object& _ao) :
  numOfAligned(_ao.numOfAligned), maximum_deletion(_ao.maximum_deletion)
{
  //first determine the size of unsigned in term of char
  unsigned sizeOfUnsigned=sizeof(unsigned)/sizeof(char);
  
  align_length=new unsigned [numOfAligned];
  memcpy(align_length, _ao.align_length, numOfAligned*sizeOfUnsigned);

  align_position=new unsigned* [numOfAligned];
  for(unsigned i=0;i<numOfAligned;i++)
    {
      align_position[i]=new unsigned [2];
      memcpy(align_position[i], _ao.align_position[i], 2*sizeOfUnsigned);
    }

  min_deletions=new unsigned [numOfAligned];
  memcpy(min_deletions, _ao.min_deletions, numOfAligned*sizeOfUnsigned);

  n_errors=new unsigned [numOfAligned];
  memcpy(n_errors, _ao.n_errors, numOfAligned*sizeOfUnsigned);

  p_region_max_length=new unsigned*[numOfAligned];
  for(unsigned i=0;i<numOfAligned;i++)
    {
      p_region_max_length[i]=new unsigned[maximum_deletion+1];
      memcpy(p_region_max_length[i], _ao.p_region_max_length[i], (maximum_deletion+1)*sizeOfUnsigned);
    }

  excess_error_positions=new unsigned*[numOfAligned];
  for(unsigned i=0;i<numOfAligned;i++)
    {
      excess_error_positions[i]=new unsigned[AlignmentSettings::negative_excess_deletions_max];
      memcpy(excess_error_positions[i], _ao.excess_error_positions[i],AlignmentSettings::negative_excess_deletions_max * sizeOfUnsigned);
    }

  alleles_all=new unsigned[numOfAligned];
  memcpy(alleles_all, _ao.alleles_all, numOfAligned*sizeOfUnsigned);

  alleles_from_distinct_genes=new unsigned[numOfAligned];
  memcpy(alleles_from_distinct_genes, _ao.alleles_from_distinct_genes, numOfAligned*sizeOfUnsigned);
  //done!!!
}

//assignment operator
Alignment_Object& Alignment_Object::operator = (const Alignment_Object& _ao)
{
  if(this==&_ao)
    return *this;

  unsigned sizeOfUnsigned=sizeof(unsigned)/sizeof(char);
  
  numOfAligned=_ao.numOfAligned;
  maximum_deletion=_ao.maximum_deletion;
  
  align_length=new unsigned [numOfAligned];
  memcpy(align_length, _ao.align_length, numOfAligned*sizeOfUnsigned);

  align_position=new unsigned* [numOfAligned];
  for(unsigned i=0;i<numOfAligned;i++)
    {
      align_position[i]=new unsigned [2];
      memcpy(align_position[i], _ao.align_position[i], 2*sizeOfUnsigned);
    }

  min_deletions=new unsigned [numOfAligned];
  memcpy(min_deletions, _ao.min_deletions, numOfAligned*sizeOfUnsigned);

  n_errors=new unsigned [numOfAligned];
  memcpy(n_errors, _ao.n_errors, numOfAligned*sizeOfUnsigned);

  p_region_max_length=new unsigned*[numOfAligned];
  for(unsigned i=0;i<numOfAligned;i++)
    {
      p_region_max_length[i]=new unsigned[maximum_deletion+1];
      memcpy(p_region_max_length[i], _ao.p_region_max_length[i], (maximum_deletion+1)*sizeOfUnsigned);
    }

  excess_error_positions=new unsigned*[numOfAligned];
  for(unsigned i=0;i<numOfAligned;i++)
    {
      excess_error_positions[i]=new unsigned[AlignmentSettings::negative_excess_deletions_max];
      memcpy(excess_error_positions[i], _ao.excess_error_positions[i],AlignmentSettings::negative_excess_deletions_max * sizeOfUnsigned);
    }

  alleles_all=new unsigned[numOfAligned];
  memcpy(alleles_all, _ao.alleles_all, numOfAligned*sizeOfUnsigned);

  alleles_from_distinct_genes=new unsigned[numOfAligned];
  memcpy(alleles_from_distinct_genes, _ao.alleles_from_distinct_genes, numOfAligned*sizeOfUnsigned);
  //done!!!
  
  return *this;
}

void Alignment_Object::ResetData()
{
  //now we need to delete the data if there are 
  //we better do this in a reverse order than it is constructed
  
  //cout<<"reset data:"<<endl;
  //
  //cout<<"alleles_from_distinc_gene"<<endl;
  if(alleles_from_distinct_genes!=NULL)
    {
      //cout<<"alleles_from_distinc_gene"<<endl;
      delete[] alleles_from_distinct_genes;
      alleles_from_distinct_genes=NULL;
    }
  //
  //cout<<"alleles_all"<<endl;
  if(alleles_all!=NULL)
    {
      delete[] alleles_all;
      alleles_all=NULL;
    }
  //
  //cout<<"excess error"<<endl;
  if(excess_error_positions!=NULL)
    {
      for(unsigned i=0;i<numOfAligned;i++)
	{
	  if( excess_error_positions[i]!=NULL)
	    {
	      delete [] excess_error_positions[i];
	    }
	}
      delete [] excess_error_positions;
      excess_error_positions=NULL;
    }
  //
  //cout<<"p_region"<<endl;
  if(p_region_max_length!=NULL)
    {
      for(unsigned i=0;i<numOfAligned;i++)
	{
	  if(p_region_max_length[i]!=NULL)
	    {
	      delete [] p_region_max_length[i];
	    }
	}
      delete [] p_region_max_length;
      p_region_max_length=NULL;
    }
  //
  //cout<<"error position"<<endl;
  if(error_positions!=NULL)
    {
      for(unsigned i=0;i<numOfAligned;i++)
	{
	  if(error_positions[i]!=NULL)
	    {
	      delete [] error_positions[i];
	    }
	}
      delete [] error_positions;
      error_positions=NULL;
    }
  //
  //cout<<"n_errors"<<endl;
  if(n_errors!=NULL)
    {
      delete[] n_errors;
      n_errors=NULL;
    }
  //
  //cout<<"min_deletions"<<endl;
  if(min_deletions)
    {
      delete [] min_deletions;
      min_deletions=NULL;
    }
  //
  //cout<<"align_position"<<endl;
  if(align_position!=NULL)
    {
      for(unsigned i=0;i<numOfAligned;i++)
	{
	  if(align_position[i]!=NULL)
	    {
	      delete [] align_position[i];
	    }
	}
      delete [] align_position;
      align_position=NULL;
    }
  //
  //cout<<"align length"<<endl;
  if(align_length!=NULL)
    {
      //cout<<"***align_length"<<endl;
      delete[] align_length;
      align_length=NULL;
    }
  //cout<<"endl"<<endl;
  numOfAligned=0;
  maximum_deletion=0;
  
}  
  
Alignment_Object::~Alignment_Object()
{
  //calling to 
  //cout<<"***calling to destruct alignment"<<endl;
  ResetData();
}
/*define the text output of the alignment object
 */
string Alignment_Object::toString()
{
  stringstream ss;
  ss<<">>Alignment object:\n";
  //starting putting everything in order
  ss<<"align_length:";
  for(unsigned i=0;i<numOfAligned;i++)
    {
      ss<<align_length[i];
      if(i!=numOfAligned-1)
	ss<<"\t";
      else
	ss<<"\n";
    }
  ss<<"align_position:\n";
  for(unsigned i=0;i<numOfAligned;i++)
    {
      ss<<"\t";
      for(unsigned j=0;j<2;j++)
	{
	  ss<<align_position[i][j];
	  if(j!=1)
	    ss<<"\t";
	  else
	    ss<<"\n";
	}
    }
  ss<<"min_deletions:";
  for(unsigned i=0;i<numOfAligned;i++)
    {
      ss<<min_deletions[i];
      if(i!=numOfAligned-1)
	ss<<"\t";
      else
	ss<<"\n";
    }

  ss<<"n_errors:";
  for(unsigned i=0;i<numOfAligned;i++)
    {
      ss<<n_errors[i];
      if(i!=numOfAligned-1)
	ss<<"\t";
      else
	ss<<"\n";
    }

  ss<<"error_positions:\n";
  for(unsigned i=0;i<numOfAligned;i++)
    {
      ss<<"\t";
      for(unsigned j=0;j<n_errors[i];j++)
	{
	  ss<<error_positions[i][j];
	  if(j!=n_errors[i]-1)
	    ss<<"\t";
	  else
	    ss<<"\n";
	}
    }

  ss<<"p_region_max_length:\n";
  for(unsigned i=0;i<numOfAligned;i++)
    {
      ss<<"\t";
      for(unsigned j=0;j<maximum_deletion+1;j++)
	{
	  ss<<p_region_max_length[i][j];
	  if(j!=maximum_deletion+1-1)
	    ss<<"\t";
	  else
	    ss<<"\n";
	}
    }
  ss<<"excess_error_positions:\n";
  for(unsigned i=0;i<numOfAligned;i++)
    {
      ss<<"\t";
      for(unsigned j=0;j<AlignmentSettings::negative_excess_deletions_max;j++)
	{
	  ss<<excess_error_positions[i][j];
	  if(j!=AlignmentSettings::negative_excess_deletions_max-1)
	    ss<<"\t";
	  else
	    ss<<"\n";
	}
    }
  
  ss<<"alleles_all:";
  for(unsigned i=0;i<numOfAligned;i++)
    {
      ss<<alleles_all[i];
      if(i!=numOfAligned-1)
	ss<<"\t";
      else
	ss<<"\n";
    }
  
  ss<<"alleles_from_distinct_genes:";
  for(unsigned i=0;i<numOfAligned;i++)
    {
      ss<<alleles_from_distinct_genes[i];
      if(i!=numOfAligned-1)
	ss<<"\t";
      else
	ss<<"\n";
    }
  
  return ss.str();  
}

//this is the code to finds highest scoring alignement of seq1 and seq2 that forces their
//left ends to match.
//this is the one in the original code named "align_with_constraints_fixed_left_remove_right_errors". 
//we named it so here in order to make the distinction. no further changes made.

//input:
// _seq1 and _seq2, the sequence strings to align, with forcing to begin from the left
// _maximum_errors and error_cost
//output: Note: caller need to initialize the memory for the output
// _n_errors, pointer to the variable for the output,
// _error_positions, array of size maximum_errors, but only the first _n_errors are used.
//             error positions are relative to the beginning of the
//             alignment starts
//return: the align_length.
unsigned align_with_constraints_fixed_left_remove_right_errors
       (const string& _seq1, const string& _seq2, const unsigned& _maximum_errors, 
	const double& _error_cost,
	/*output*/ unsigned* _n_errors, unsigned* _error_positions)
{
  unsigned l_seq1=_seq1.size();
  unsigned l_seq2=_seq2.size();
  unsigned align_length=0;

  //cout<<"\t\tl_seq1:"<<l_seq1<<endl;
  //cout<<"\t\tl_seq2:"<<l_seq2<<endl;

  unsigned max_match_length=l_seq1;
  if(l_seq2<l_seq1)
    {
      max_match_length=l_seq2;
    }
  //now we need to bitwise matching between two strings.
  //there are a few things to do for this
  //1) go through the string to find match and mismatch for each bit
  //    to make it right, we need to go one beyond the maximum allowed errors
  //    till we reach the one next error if possible
  //2) we need to figure out i)whether this is empty array ii)error positions
  //    iii)align length;iv)num of errors v)best score so far
  //3) find the best score along the longest match with maximum errors. Just
  //    go through each errors position and calculate the best score for it and 
  //    compare to return the best one.
  
  bool matchFlag;
  
  *_n_errors=0;
  unsigned i;
 
  //cout<<"*****calling ...cost:"<<_error_cost<<endl;
  //do the alignment
  for(i=0;i<max_match_length;i++)
    {
      
      //check for match or mismatch
      if(_seq1.at(i)==_seq2.at(i))
	{
	  matchFlag=true;
	}
      else //mismatch, need to update the n_errors and error_position
	{
	  matchFlag=false;
	
	  *_n_errors=*_n_errors+1;
	}

      //check to see whether we have too many errors
      if(*_n_errors>_maximum_errors)
	{
	  //we are reaching the maximum errors, now it is _maximum_errors+1, done!
	  *_n_errors=*_n_errors-1;//here don't count the last one
	  break;
	}
      //this one is a mismatch, but not go beyond the maximum allowed ones, we remeber the position
      if(!matchFlag)
	{
	  _error_positions[*_n_errors-1]=i;
	}
      else
	{
	  align_length=i+1;
	}
    }
  //cout<<_seq1<<endl;
  //cout<<_seq2<<endl;
  //cout<<"\tafter match the align_length:"<<align_length<<endl;
  //cout<<"after the first round, alength:"<<align_length<<";n_error:"<<*_n_errors<<endl;
  
  //so far we know the longest alignment based on either it is not long <maximum_length
  //or because we have too many errors
  //now we need to recalculate the n_errors
  /*
  if(i>max_match_length)
    align_length=i; //because i starts at zero, so it is i-1+1
  else //due to maximum error
    align_length=i+1;//because it starts at zero, so it i+1;
  
  */
  //check to make sure the n_errors is not larger than maximum_errors
  //and we did not end at errors
  if(*_n_errors>0)
    {
      for(unsigned j=*_n_errors-1;j>0;j--)
	{
	  //from back to front
	  if(_error_positions[j]<align_length)
	    {
	      *_n_errors=j+1;
	      break;//we are done. the one error right before the last aligned one was found
	    }
	}
    }
  //cout<<"second round"<<endl;
  //cout<<"after the second round, alength:"<<align_length<<";n_error:"<<*_n_errors<<endl;
  //so far the align_length and n_error were correct.
  //now we need to figure out best score along the longest match
  double best_score=align_length-_error_cost* (*_n_errors);
  double running_score=best_score;
  unsigned best_n_errors=*_n_errors;

  //cout<<"after the 3 round, best_score:"<<best_score<<endl;//<<";n_error:"<<*_n_errors<<endl;
  for(unsigned j=0;j< *_n_errors ; j++)
    {
      //cout<<"\t@@@sub loop(#errors):"<<j<<endl;
      running_score=(_error_positions[j]-1+1)-_error_cost*j;
      //cout<<"\t@@@running score:"<<running_score<<endl;
      if(running_score>=best_score)
	{
	  align_length=_error_positions[j]-1+1;
	  best_score=running_score;
	  best_n_errors=j;
	  //cout<<"\t\t@@@@got one"<<endl;
	}
      //cout<<"why??j:"<<j<<endl;
    }
  *_n_errors=best_n_errors;
  //cout<<"after the 3 round, best_score:"<<best_score<<";align_length:"<<align_length<<";n_error:"<<*_n_errors<<endl;
  //prepare the error_position array based on the 
  //this is not necessary, we just need to according to the _n_errors to 
  //check the error_positions array elements
  for(unsigned j=*_n_errors;j<_maximum_errors;j++)
    {
      _error_positions[j]=0;
    }

  //clean up
  //delete [] c;
  //cout<<"done"<<endl;
  return align_length;
}

//input:
// _seq, the input sequence to align
// _target, the genomice template to align against.
// _maximium_errors, the maximum number of errors allow in this alignment
// _error_cost, the error cost, double
//
//output: Note, the caller need to initialize the memory
// _align_position, the array of size 2 to contain the alignment starting position 
//         for the best alignment. [0] for _seq, [1] for _target
// _n_errors, the pointer to the unsigned variable for holding the number of errors in the best alignment
// _error_positions, the array of size _maximum_errors. but it only have n_errors elements. Again, need to
//         initialized by the caller. 

//return: 
//   alignment_length, unsigned
unsigned align_with_constraints_fast_left
   (const string& _seq, const string& _target, 
    const unsigned& _maximum_errors,  const double& _error_cost, 
    /*output*/unsigned* _align_position, unsigned* _n_errors, 
					 unsigned* _error_positions)
{
  unsigned align_length=0;
  
  unsigned l_seq=_seq.size();
  unsigned l_target=_target.size();
  unsigned best_match_index=0;
  double score=-1E200;

  unsigned current_align_length, current_n_errors;
  unsigned* current_error_positions=new unsigned [_maximum_errors];
  double current_score;
  unsigned max_length;

  //cout<<"param:maximum_errors:"<<_maximum_errors<<";cost:"<<_error_cost<<endl;
  //looping through the positions of sequence to force the left side of target to align with
  //subsequence of seq starting at the current index of _seq
  for(unsigned i=0;i<l_seq;i++)
    {
      //cout<<"*************doing loop "<<i<<endl;
      max_length=l_seq-i;
      if(l_target<max_length)
	{
	  max_length=l_target;
	}
      //can not be better
      if(max_length<score)
	{
	  break;
	}
      string current_seq=_seq.substr(i, max_length);

      //string current_target=_target.substr(0,max_length);
      string current_target=_target;
      //cout<<"\tseq:"<<current_seq<<endl;
      //cout<<"\ttarget:"<<current_target<<endl;
      //call to do alignment forcing fixed left ends
      current_align_length=
	align_with_constraints_fixed_left_remove_right_errors
	(current_seq, current_target, _maximum_errors,_error_cost,
	 &current_n_errors, current_error_positions);
      
      current_score=current_align_length-_error_cost*current_n_errors;
      //cout<<"\tcurrent_score:"<<current_score<<";current_align_length:"<<current_align_length<<endl;

      //we got a better one or got an identical score, but long, we are good
      if(current_score>score||
	 ((current_score-score<1E-6||current_score-score>-1E-6)&&(current_align_length>align_length))
	 )
	{
	  align_length=current_align_length;
	  *_n_errors=current_n_errors;
	  //copy over the corrent error positions
	  for(unsigned j=0;j<current_n_errors;j++)
	    {
	      _error_positions[j]=current_error_positions[j];
	    }
	  score=current_score;
	  best_match_index=i;
	}//if loop
      
    }//end of looping over the seq positions
  
  _align_position[0]=best_match_index;
  _align_position[1]=0;
  
  //clean up
  delete [] current_error_positions;
  return align_length;
}



//see above for the definition (match_Vs)
//we should keep in mind that we are doing multithread functions
//remember Alignment_Obj* J has been declared and initialized, 
//but the vector memeber of it is not initialized. Do we have to 
//initialized it first?? or not necessarily? we could do the 
//initialization when we need to. 
bool match_J(const SequenceString& _seq,
	      const GenomicJ* _genJs, const unsigned& _numOfJSegs, 
	      const unsigned& _J_minimum_alignment_length, const unsigned& _J_maximum_deletion, 
	      const unsigned& _negative_excess_deletion_max, const unsigned& _J_allowed_errors, 
	     const unsigned& _error_cost,
	     /*output*/ Alignment_Object& _J
 )
{
  //***from Matlab code
  //#codegen
  //  %This function goes through each genomic J, finds highest scoring alignment
  //  %that begins where the J sequencing primer ends, sorts the alleles
  //  %by alignment score and then returns the ones that pass certain filters.
  //***above is from Matlab code
  //here, we might need to change a little bit. In this case, we want the genJ is fixed
  //on the right, so the alignment has to start at the end of the genomic template.
  //but the input sequence start at variable position since it contain some portion of IgG/M
  //constant regions

  //now we assume the return Aligment_obj is initialized, but the data they point to are filled and initialized here in this function. so we can use
  //the pointers to point the data and use them
  //***1) unsigned* align_length, size numOfJSegs maximumly ; ---temp_align_length//% contains the length of the alignment

  //***2) align_position, 2D, numOfJSegs by 2
  //% align_position(1,j) contains the first matched nt in seq (redundant, same as length(seq) - align_length +1 ),
  //% align_position_(2,j) is the first matched nt in genJ.seq{j} - also redundant information because we know where the sequencing primer ends.
  // temp_align_position = J->align_position;

  //***3) min_deletions, 1D, maximumly numOfJSegs long, temp_min_deletions;

  //***4)n_errors,  temp_n_errors=J->n_errors;//zeros(length(genJ),1); % number of mismatches in alignment

  //***5)error_positions = J->error_positions ;//zeros(length(genJ),J_allowed_errors); % positions of mismatches
  //     2D, maximumly numOfJSegs by J_allowed_errors; (numOfAligned x n_errors[i])

  //from here, all the array are of size numOfAligned, and declared/initialized after the alignment all done.
  //***6)p_region_max_length, temp_p_region_max_length 2D, numOfAligned by (max_deletion+1). 
  //    2ndDimension is fixed to be max_deletion+1. 1st dimension could be shorter	
  
  //***7)excess_error_positions, 2D, maximumly numOfJSegs by negative_deletion_max (3), usually, could be shorter, but keep it fixed

  //***8)alleles_all, 1D numOfAligned
  //***9)alleles_from_distinct_gene, 1D numOfAligned
  /*
  cout<<"***first****:sequence:"<<_seq.GetSequence()<<endl;
  cout<<"\t_numOfJSegs:"<<_numOfJSegs
      <<"\n\t_J_minmum_alignment_length:"<<_J_minimum_alignment_length
      <<"\n\t_J_maximum_deletion:"<<_J_maximum_deletion
      <<"\n\t_negative_excess_deletion_max:"<<_negative_excess_deletion_max
      <<"\n\t_J_allowed_errors"<<_J_allowed_errors
      <<"\n\t_error_cost"<<_error_cost<<endl;*/
  bool* j_large_deletion_flag=new bool [_numOfJSegs]; //=zeros(length(genJ),1); % Flag if deletions is too large
  unsigned l_seq = _seq.GetLength(); //length of the input sequence
  //cout<<"l_seq:"<<l_seq<<endl;
  unsigned* temp_align_length=new unsigned[_numOfJSegs];
  unsigned** temp_align_position=new unsigned*[_numOfJSegs];
  for(unsigned i=0;i<_numOfJSegs;i++)
    {
      temp_align_position[i]=new unsigned[2];
    }
  unsigned* temp_min_deletions=new unsigned[_numOfJSegs];

  unsigned* temp_n_errors=new unsigned[_numOfJSegs];
  std::memset(temp_n_errors, -1, _numOfJSegs);//fill the default value with -1
  unsigned** temp_error_positions=new unsigned* [_numOfJSegs];
  for(unsigned i=0;i<_numOfJSegs;i++)
    {
      temp_error_positions[i]=new unsigned[_J_allowed_errors];
    }
  //unsigned** temp_p_region_max_length=NULL;
  //unsigned** temp_excess_error_positions=NULL;
  //unsigned* temp_alleles_all; unsigned* alleles_from_distinct_gene;
  //temp_alleles_all=NULL; alleles_from_distinct_gene=NULL;

  //cout<<"***f2irst****"<<endl;
  
  //the following are only for setting up the output for the 
  //unsigned align_length_func;
  unsigned align_position_func[2];
  unsigned* error_position_func=new unsigned[_J_allowed_errors];
  SequenceString rev_seq=FlipSequenceString(_seq);
  //unsigned n_errors_func;
  //% Loop through template alleles for alignment
  //cout<<"***first****3"<<endl;
  for(unsigned int i=0;i<_numOfJSegs;i++) //for j=1:length(genJ)
    {
      //cout<<"\t***loop****"<<i<<endl;
     //%j=1
     //%disp(['loop: ' num2str(j)])
     //% Get highest scoring alignment for this allele, with acceptable number of errors
     //% NOTE: I use an alignment function that fixes position on the left (as in match_Vs) but I pass in the
     //% reversed sequences, because I want to fix the position on the right.

     //unsigned l_target=_genJs[i].Get_Seq().GetLength();//length of the genomic template allele for current one

     //before doing the alignment, we need to flip/reverse the sequence, since the function do it assume doing
     //fixed on the left. 
     
     SequenceString rev_target=FlipSequenceString(_genJs[i].Get_Seq());
     unsigned l_target =rev_target.GetLength();
     //cout<<"l_target:"<<l_target<<endl;
     //cout<<"\t&&&doing alignment :"<<endl;
     //cout<<"\trev_seq:"<<rev_seq.toString()<<endl;
     //cout<<"\trev_target:"<<rev_target.toString()<<endl;
     //now calling to do the alignment
     temp_align_length[i]=
       align_with_constraints_fast_left(rev_seq.GetSequence(), rev_target.GetSequence(), _J_allowed_errors, _error_cost,
					align_position_func, temp_n_errors+i, error_position_func);
					  
     //cout<<"\ttemp_align_length["<<i<<"]:"<<temp_align_length[i]<<endl;
     //cout<<"\ttemp_n_errors[i]"<<i<<"]:"<<temp_n_errors[i]<<endl;
     //cout<<"\talign_position_func"<<align_position_func[0]<<","<<align_position_func[1]<<endl;
     
     //_J->align_length.push_back(temp_align_length);
     
     //[align_length(j), rev_align_position, n_errors(j), rev_error_positions]=align_with_constraints_fast_no_fix(seq(end:-1:1) , genJ(j).seq(end:-1:1) , J_allowed_errors, error_cost);
     
     //_J->n_errors.push_back(temp_n_errors);

     // % Reverse the error positions to be relative to left-right orientation.
     //********be careful here, we decide to keep the j error position relative to the END OF J CHAIN
     //  %%%j  error positions is relative to the end of J chain
     //vector<unsigned> temp_error_position_vec;
     if( temp_n_errors[i]>0&&temp_n_errors[i]!=((unsigned)-1))  //the second case will not possible??
       {
	 //unsigned rev_j=temp_n_errors-1;
	 for(unsigned j=0;j<temp_n_errors[i];j++)
	   {
	     //unsigned tempValue=temp_error_positions[j]
	     temp_error_positions[i][j]= error_position_func[temp_n_errors[i]-j-1] ;
	     //temp_error_positions[rev_j]=tempValue;
	     //rev_j--;
	   }
       }//end
     //_J->error_positions.push_back(temp_error_position_vec);

     //vector<unsigned> temp_align_position_vec;
	//% Calculate (redundant) align_position values. These are used in model scripts.
     temp_align_position[i][0]
       /*temp_align_position[0]*/ = (l_seq-(align_position_func[0]+temp_align_length[i]-1)-1 ) ; //for _seq first position
     temp_align_position[i][1]
       /*temp_align_position[1]*/= l_target-(align_position_func[1]+temp_align_length[i]-1)-1 ;//for target second position
     //_J->align_position.push_back(temp_align_position_vec);

     // % Calculate deletions implied by alignment
     temp_min_deletions[i]=l_target-(align_position_func[1]+temp_align_length[i]-1) - 1;
     //cout<<"\tmin_deletion:"<<temp_min_deletions[i]<<endl;
     //    % Flag if number of deletions is too many
     j_large_deletion_flag[i]=false;
     if( temp_min_deletions[i] > _J_maximum_deletion)
	{
	  j_large_deletion_flag[i] = true;
	  //%continue;
        }//  end
     //cout<<"\tdone loop"<<endl;
   }//end of for loop to go through the J segs alleles for alignment.
    //% Check!  No need to check in c++ code.
 //	    assert(all(min_deletions>=0));

 //	  % sort the genomic Js in descending order of 'score' and ascending order of min_deletions

 //calculate the score first
  //cout<<"***first****3:_numOfJSegs:"<<_numOfJSegs<<endl;
  double* scores=new double[_numOfJSegs];
 //prepare the sorted index of the array.
  unsigned* sorted_index=new unsigned[_numOfJSegs];
  
  //cout<<"=========>before sorting:"<<endl;
 for(unsigned k=0;k<_numOfJSegs;k++)
   {
     sorted_index[k]=k;
     scores[k]=temp_align_length[k]-_error_cost*temp_n_errors[k];
     //cout<<scores[k]<<"-"<<sorted_index[k]<<"-"<<temp_min_deletions[k]<<",";
     
   }
 //cout<<"\n***first****3aa"<<endl;
 //sorted index also holding the gene allele index
 QuickSort<double>(scores, 0, _numOfJSegs-1, sorted_index, temp_min_deletions);
 //	      scores  = align_length - error_cost*n_errors;
 //S = [-scores, min_deletions];
 //	  [~,order]=sortrows(S);
 //cout<<"=========>after sorting:";
 /*for(unsigned k=0;k<_numOfJSegs;k++)
   {
     //sorted_index[k]=k;
     //scores[k]=temp_align_length[k]-_error_cost*temp_n_errors[k];
     cout<<scores[k]<<"-"<<sorted_index[k]<<",";
   }
 */
 //cout<<"\n***first****3bb"<<endl;
 //now we need to reverse the order, since the QuickSort is ascending, but for our purpose we need to descending.
 Reverse(sorted_index, _numOfJSegs);
 Reverse(scores, _numOfJSegs);
 Reverse(temp_min_deletions, _numOfJSegs);

 //cout<<"=========>after sorting:";
 /*for(unsigned k=0;k<_numOfJSegs;k++)
   {
     //sorted_index[k]=k;
     //scores[k]=temp_align_length[k]-_error_cost*temp_n_errors[k];
     cout<<scores[k]<<"-"<<sorted_index[k]<<",";
     }*/
 //cout<<"\n***first****3ccc"<<endl;
 //  % Set a score threshold for alleles. well, this is kind of arbitrary
 //we want to get the best ones, but limited numbers 
 double min_score=max_mf(scores,_numOfJSegs)-3*_J_minimum_alignment_length;
 //	      min_score = max(scores) - 3*J_minimum_alignment_length;

 //% Subset of alleles that have high enough score, long enough alignment, and
 //	      % not too many deletions.
 unsigned* ok_order = new unsigned[_numOfJSegs]; //this will directly used by J.alleles_all, so do NOT delete/clean later
 unsigned ok_count=0;
 //cout<<"***first****4"<<endl;
 for(unsigned i=0;i<_numOfJSegs;i++)
   { 
     //find( scores(order) >= min_score & align_length(order) >= J_minimum_alignment_length & j_large_deletion_flag(order)==0);
     if(scores[i]>=min_score&&temp_align_length[sorted_index[i]]>=_J_minimum_alignment_length&&!j_large_deletion_flag[sorted_index[i]])
       {
	 ok_order[ok_count]=sorted_index[i];
	 temp_min_deletions[ok_count]=temp_min_deletions[i];
	 ok_count++;
       }
   }
 /*cout<<"---------->showing the ok_order array:";
 for(unsigned i=0;i<ok_count;i++)
   {
     cout<<ok_order[i]<<",";
   }
 cout<<endl;
 */ 
 //bool seq_j_ok=true;
 //cout<<"*******first 4a check for status"<<endl;
 if(ok_count<=0)//empty
   {
     //cout<<"*****inside false condition"<<endl;
     _J.numOfAligned=0;
     //do we want to clean up the memory, not necessary???
     //we have to clean up the memory, by now some of the arrays have been allocated
     //<--------
     //clean up
     //cout<<"\t555delete 1"<<endl;
     delete [] j_large_deletion_flag;
     //cout<<"\t555delete 2"<<endl;
     delete [] temp_align_length;
     //cout<<"\t555delete 3"<<endl;
     CleanUpMemory(temp_align_position, _numOfJSegs);
     //cout<<"\t555delete 4"<<endl;
     delete [] temp_min_deletions;
     delete [] temp_n_errors;
     CleanUpMemory(temp_error_positions, _numOfJSegs);
     //CleanUpMemory(temp_p_region_max_length, _numOfJSegs);
     //CleanUpMemory(temp_excess_error_positions, _numOfJsegs);
     //delete [] temp_alleles_all;
     //delete [] alleles_from_distinc_gene;
     //cout<<"\t555delete 5"<<endl;
     delete [] error_position_func;
     //cout<<"\t555delete 6"<<endl;
     delete [] scores;
     //cout<<"\t555delete 7"<<endl;
     delete [] sorted_index;
     //cout<<"\t555delete 8"<<endl;
     delete [] ok_order; //ok order has not been referenced directly by J.alleles_all
     //and in this case, it will not. so we delete it.
 
     //delete [] genJ_ok_index;
     //cout<<">>>>>>>>>>>>failed from here1"<<endl;
     return false;
   }
 //cout<<"************first 4 bb keep going....."<<endl;
 //if we are here we are good, we still need to do more next
 //   % Store all the alignment information in J for acceptable alleles
 //NOTE: the following code to copy over the elements from one array to another
 //it returns false only because the array sizes are not set up correctly
 //if the size are right then we have no trouble. Here we are sure 
 //that the size are all right. and we will not have errors. so at the 
 //time of false, we will not clean up memeory, since it will never 
 //happen
 _J.numOfAligned=ok_count;
 _J.maximum_deletion=_J_maximum_deletion;
 //cout<<"***first****4aaa"<<endl;
 _J.align_length = new unsigned [ok_count];
 //cout<<"****ok_count***:"<<ok_count<<endl;
 if(!CopyElements(temp_align_length, _numOfJSegs,  _J.align_length, ok_count, ok_order, ok_count))
   {
     //cout<<">>>>>>>>>>>>failed from here2"<<endl;
     return false;
   }
 //cout<<"***first****4bbbb"<<endl;
 _J.align_position =new unsigned* [ok_count];
 for(unsigned m=0;m<ok_count;m++)
   {
     _J.align_position[m]=new unsigned [2];
   }
 if(!CopyElements(temp_align_position, _numOfJSegs, 2, _J.align_position, ok_count, 2, ok_order, ok_count))//(:,ok_order);
   {
     return false;
   }
 //cout<<"***first****4cccc"<<endl;
 _J.min_deletions =new unsigned [ok_count];
 /*cout<<"\t***&&&&&&showing min deletions"<<endl;
 cout<<"\t\t";
 for(unsigned _m=0;_m<_numOfJSegs;_m++)
   {
     cout<<temp_min_deletions[_m]<<","<<endl;
   }
   cout<<endl;*/
 /*if(!CopyElements(temp_min_deletions, _numOfJSegs, _J.min_deletions, ok_count, ok_order, ok_count))
   return false;
 */
 //****NOTE:important, temp_min_deletions, is different than others, it is sorted together with socres, so we only need to copy over
 //what is in the front (best values) according to ok_count
 for(unsigned _m=0;_m<ok_count;_m++)
   {
     _J.min_deletions[_m]=temp_min_deletions[_m];//<<","<<endl;
   }
 //cout<<"***first****4dddd"<<endl;
 _J.n_errors = new unsigned [ok_count];
 if(!CopyElements(temp_n_errors, _numOfJSegs, _J.n_errors, ok_count, ok_order, ok_count))
   return false;
 //cout<<"***first****4eee"<<endl;
 _J.error_positions = new unsigned* [ok_count];
 for(unsigned m=0;m<ok_count;m++)
   {
     _J.error_positions[m]=new unsigned [_J_allowed_errors];
   }
 //cout<<"***first****4fff"<<endl;
 if(!CopyElements(temp_error_positions, _numOfJSegs, _J_allowed_errors, _J.error_positions, ok_count, _J_allowed_errors, ok_order, ok_count))
   return false;//error_positions(ok_order,:);

 //now we are done so, and need to run initialization for the later assignments
 _J.p_region_max_length = new unsigned* [ok_count];

 //cout<<"***first5****"<<endl;
 for(unsigned i=0;i<ok_count;i++)
   {
     _J.p_region_max_length[i]=new unsigned[_J_maximum_deletion+1];
     std::memset(_J.p_region_max_length[i], 0, (_J_maximum_deletion+1)*sizeof(unsigned)/sizeof(char));//fill the default value with 0
   }
 /*cout<<"right after setting the memory:maxlength:"<<_J_maximum_deletion<<endl;
for(unsigned i=0;i<ok_count;i++)
   {
     cout<<"\t"<<endl;
     for(unsigned j=0;j<_J_maximum_deletion+1;j++)
       {
	 cout<<_J.p_region_max_length[i][j]<<",";
       }
     cout<<endl;
     }*/
 //zeros(numel(ok_order),J_maximum_deletion+1);
 _J.excess_error_positions =new unsigned* [ok_count];
 //zeros(numel(ok_order),negative_excess_deletions_max);
 for(unsigned i=0;i<ok_count;i++)
   {
     _J.excess_error_positions[i]=new unsigned [_negative_excess_deletion_max];
   }
 _J.alleles_all = ok_order;
 
 // % Find gene indices of the acceptable alleles.
 //		      %%%codegen compatible version of
 //		  %%%genJ_ok_gene_inds = [genJ(ok_order).gene_index];
 unsigned* genJ_ok_index = new unsigned[ok_count];
 //go through the genJ to figure out the distinct genes for the each allele, and then reture a array of indices to them
 
 //cout<<"\t^^^^^showing the gene index of the alleles:";
 for(unsigned i=0;i<ok_count;i++)
   {
     genJ_ok_index[i]=_genJs[_J.alleles_all[i]].Get_GeneIndex();
     // cout<<genJ_ok_index[i]<<",";
   }
 //cout<<endl;
 //genJ_ok_gene_inds = zeros(1,numel(ok_order));
//first we need to sort the J.alleles_all array 
 unsigned* sorted_genJ_ok_index_temp=new unsigned [ok_count];
 unsigned* sorted_genJ_ok_index_temp_index=new unsigned[ok_count];
 unsigned numOfUnique;
 Unique(genJ_ok_index, ok_count, sorted_genJ_ok_index_temp, sorted_genJ_ok_index_temp_index,numOfUnique);
 //cout<<"after unique:numOfUnique:"<<numOfUnique<<endl;
 //now we got the distinct gene index index, just need to sort it, to make it in order
 //because when we do the unique, we first sort it and the gene index might not in order
 //so the index of the gene index is not in order either
 QuickSort(sorted_genJ_ok_index_temp_index,0, numOfUnique-1);
 //cout<<"after sorting...."<<endl;
 //cout<<"the index index:"<< sorted_genJ_ok_index_temp_index[0]<<endl;

 //cout<<"***first****6"<<endl;
 //compare and then get it copied to the alignment object
 _J.alleles_from_distinct_genes=new unsigned[ok_count];

 //unsigned runningValue=genJ_ok_index[0];
 //unsigned runningIndexOfAlleleArray=1;
 //_J.alleles_from_distinct_genes[0]=_J.alleles_all[0];
 for(unsigned i=0;i<numOfUnique;i++)
   {
     //if(runningValue!=genJ_ok_index[i])
     //  {
	 _J.alleles_from_distinct_genes[i] = ok_order[sorted_genJ_ok_index_temp_index[i]];
	 //runningValue=genJ_ok_index[i];
	 //runningIndexOfAlleleArray++;
	 //}
   }
 //need to fill the rest of the distinc gene array with -1 as an indicator
 for(unsigned i=numOfUnique;i<ok_count;i++)
   {
     _J.alleles_from_distinct_genes[i]=-1;
   }
 //
 //% Get list of distinct genes from acceptable alleles.
 //[~,Jg] = unique(genJ_ok_gene_inds,"first" ); % Pick the first allele listed for each gene.
 //Jg = sort(Jg); % sort it because unique returns in ascending order of ok_order. This puts it back in order of best match.
 //% Only one (best) allele from each gene
 //J.alleles_from_distinct_genes = ok_order(Jg'); '
 //done so far
 //cout<<"******second ****6aaaa"<<endl; 
 //now call to generate panlindrol and negative access errors
 DeterminePalindromAndExcessError_J
   ( _seq, _genJs, ok_order,_J.min_deletions, 
     _negative_excess_deletion_max, _J_maximum_deletion,
     _J.align_length, _J.numOfAligned, _J.align_position,
     _J.p_region_max_length, _J.excess_error_positions  
     );
 //cout<<"***first****7"<<endl;
 //clean up
 delete [] j_large_deletion_flag;
 
 delete [] temp_align_length;
 CleanUpMemory(temp_align_position, _numOfJSegs);
 delete [] temp_min_deletions;
 delete [] temp_n_errors;
 CleanUpMemory(temp_error_positions, _numOfJSegs);
 //CleanUpMemory(temp_p_region_max_length, _numOfJSegs);
 //CleanUpMemory(temp_excess_error_positions, _numOfJsegs);
 //delete [] temp_alleles_all;
 //delete [] alleles_from_distinc_gene;
 
 delete [] error_position_func;

 delete [] scores;
 delete [] sorted_index;
 //delete [] ok_order; ok order is referenced directly by J.alleles_all
 
 delete [] genJ_ok_index;
 delete [] sorted_genJ_ok_index_temp;
 delete [] sorted_genJ_ok_index_temp_index;
 return true;
}

void CleanUpMemory(unsigned** p_mem, unsigned size_first_dim)
{
  for(unsigned i=0;i<size_first_dim;i++)
    {
      delete [] p_mem[i];
    }
  delete [] p_mem;
}

void DeterminePalindromAndExcessError_J
( const SequenceString& _seq, const GenomicJ* _genJs, const unsigned* _ok_order,
  const unsigned* _min_deletions, 
  const unsigned& _negative_excess_deletions_max, const unsigned& _J_maximum_deletion,
  const unsigned* _align_length, const unsigned& _numOfAligned, unsigned** _align_positions,
  /*output*/ unsigned** _p_region_max_length, unsigned** _excess_error_position  
  )
{
  unsigned p=0;
  /*  //cout<<"in the beginning p_Region_max_length:"<<endl;
for(unsigned i=0;i<_numOfAligned;i++)
   {
     //cout<<"\t"<<endl;
     for(unsigned j=0;j<_J_maximum_deletion+1;j++)
       {
	 cout<<_p_region_max_length[i][j]<<",";
       }
     cout<<endl;
   }
  */
  bool still_palindrome=true;
  //go through and find palindromic nucleotides for various number of deletions,
  //    %as well as error positions for 'negative' deletions.
  //cout<<"\tinside palindro.....before loop"<<endl;
  for(unsigned  j=0;j<_numOfAligned;j++)
    {
      //cout<<"\t*^^^^^^loop "<<j<<endl;
      string target=_genJs[_ok_order[j]].Get_Sequence();
      //% Loop over number of deletions (nd is actual number of genomic deletions).
      //  % Lower bound is maximum of 0 and min_deletions - negative_excess_deletions_max (usually 3).
      //  % Upper bound is minimum of J_maximum_deletion (usually 12) and the whole alignment being random insertions
      int nd=_min_deletions[j]-_negative_excess_deletions_max;
      if(nd<0)
	nd=0;
      int max_nd=_J_maximum_deletion;
      if(max_nd>(signed)(_align_length[j]+_min_deletions[j]))
	max_nd=_align_length[j]+_min_deletions[j];
      //cout<<"\t\t^^^^^second loop before, nd:"<<nd<<endl;
      for(;nd<=max_nd;nd++)
	{
	  //% For each value of deletions, find longest half-palindrome
          //  % from the implied end of the gene sequence
	  p=0;
	  still_palindrome=true;
	  //% Upper bound on p is the length of the J match (accounting for deletions) OR the never actually possible case of sequence length to the left of J match (accounting for deletions)
	  int max_p_length=(signed)(_align_length[j]-(nd-_min_deletions[j]));
	  if(max_p_length>(signed)(_align_positions[j][0]-1+nd-_min_deletions[j]))
	    max_p_length=(signed)(_align_positions[j][0]-1+nd-_min_deletions[j]);
	  if(max_p_length<0)
	    max_p_length=0;
	  while( still_palindrome && (signed)p < max_p_length)
	    {
	      
	      still_palindrome = target.at(_align_positions[j][1] + p + (nd - _min_deletions[j] )) == DnaComplement(_seq.GetSequence().at(_align_positions[j][0] - p + nd - _min_deletions[j] - 1) );
	      if (still_palindrome)
		{
		  p ++;
		}
	    }
	  
	  //% Store length of longest half-palindrome for each value of
	  //% deletions.
	  _p_region_max_length[j][ nd] = p;
	}
      //cout<<"\t\tend of second loop"<<endl;
      //        % Now for the 'negative' deletions, store where the mismatches are.
      // % This is so we can count number of errors for these possibilities.
      
      //        % Number of 'negative' deletions to consider:
      unsigned tempArry[]={_negative_excess_deletions_max, _min_deletions[j], _align_positions[j][0]-1 };
      unsigned n_excess = min_mf(tempArry,3);
      //cout<<"n_excess :"<<n_excess<<endl;
      unsigned k;    
      unsigned runningIndex_excessError=0;
      //cout<<"\tready to do the excess error"<<endl;
      // % Find mismatches between sequence and genomic J at those positions.
      for( k=0;k<n_excess;k++)
	{
	  //cout<<"\t\tthird for loop:"<<k<<endl;
	  //cout<<"\t\tposition:"<<_align_positions[j][0]<<endl;
	  if(target.at(_align_positions[j][1]-n_excess+k)!=_seq.GetSequence().at(_align_positions[j][0]-n_excess+k))
	    {
	      //cout<<"\t\tfound one!!pos:"<<_align_positions[j][0]-n_excess+k<<endl;
	      _excess_error_position[j][runningIndex_excessError]=_align_positions[j][0]-n_excess+k;
	      runningIndex_excessError++;
	    }
	  //cout<<"\t\t$$$$end of third for loop"<<endl;
	}
      //cout<<"runningIndex_excessError:"<<runningIndex_excessError<<endl;
      for(;runningIndex_excessError<_negative_excess_deletions_max;runningIndex_excessError++)
	{
	  //cout<<"\t\tfor loop"<<endl;
	  _excess_error_position[j][runningIndex_excessError]=0;
	}
      //cout <<"\t end of first for loop"<<endl;     
      //excess_err_pos = J.align_position(1,j) - ( n_excess + 1 - find((genJ(ok_order(j)).seq((J.align_position(2,j)-n_excess): (J.align_position(2,j)-1) )~=seq((J.align_position(1,j) -n_excess): (J.align_position(1,j)-1) ))) );
      //  J.excess_error_positions(j, 1:length(excess_err_pos)) = excess_err_pos;
    }//end of outer for loop for all the aligned strings
  //cout<<"Done for the function"<<endl;
  //return true;
}

void Alignment_Object::Serialize(ofstream& _ofs) const
{
  //do necessary checking
  if(!_ofs.is_open())
    {
      cout<<"**ERROR**: closed file buffer. quit..."<<endl;
      exit(-1);
    }
  //serializing......
  char* p_char; //pointer used to direct the writing to the file stream
  //save numOfAligned 
  p_char=(char*)&numOfAligned;
  _ofs.write(p_char, sizeof(unsigned));

  if(numOfAligned==0) //nothing to save, only save one this field
    {
      return ;
    }
  
  //save numOfAligned 
  p_char=(char*)&maximum_deletion;
  _ofs.write(p_char, sizeof(unsigned));

  //align_length, array
  p_char=(char*)align_length;
  _ofs.write(p_char, sizeof(unsigned)*numOfAligned);

  //align_position, 2D array
  for(unsigned i=0;i<numOfAligned;i++)
    {
      //write each elements
      p_char=(char*)align_position[i];
      _ofs.write(p_char, sizeof(unsigned)*2);
    }

  //min_deletion, array
  p_char=(char*)min_deletions;
  _ofs.write(p_char, sizeof(unsigned)*numOfAligned);
  
  //n_errors
  p_char=(char*)n_errors;
  _ofs.write(p_char, sizeof(unsigned)*numOfAligned);

  //error_positions
  for(unsigned i=0;i<numOfAligned;i++)
    {
      //write each elements
      p_char=(char*)error_positions[i];
      _ofs.write(p_char, sizeof(unsigned)*n_errors[i]);
    }
  
  //p_region_max_length
  for(unsigned i=0;i<numOfAligned;i++)
    {
      //write each elements
      p_char=(char*)p_region_max_length[i];
      _ofs.write(p_char, sizeof(unsigned)*(1+maximum_deletion));
    }

  //excess_error_positions
  for(unsigned i=0;i<numOfAligned;i++)
    {
      //write each elements
      p_char=(char*)excess_error_positions[i];
      _ofs.write(p_char, sizeof(unsigned)*AlignmentSettings::negative_excess_deletions_max);
    }

  //alleles_all
  p_char=(char*)alleles_all;
  _ofs.write(p_char, sizeof(unsigned)*numOfAligned);
  
  //n_errors
  p_char=(char*)alleles_from_distinct_genes;
  _ofs.write(p_char, sizeof(unsigned)*numOfAligned);

  //done!!
}

void Alignment_Object::Deserialize(ifstream& _ifs)
{//the order of reading is important!!!
  unsigned original_numOfAligned=numOfAligned;
  //first read numOfAligned
  char* p_char=(char*)&numOfAligned;
  if(_ifs.eof())
    {
      cout<<"****ERROR: the file end reached before reading the first field of the object"<<endl;
      exit(-1);
    }
  _ifs.read(p_char, sizeof(unsigned));
  if(numOfAligned==0)//nothing to read for this current object,
    {
      return ;
    }

  //we are here, means we need to read more
  p_char=(char*)&maximum_deletion;
  _ifs.read(p_char, sizeof(unsigned));

  //align_length
  if(align_length!=NULL)
    {
      delete[] align_length;
    }
  align_length=new unsigned[numOfAligned];
  p_char=(char*)align_length;
  _ifs.read(p_char, sizeof(unsigned)*numOfAligned);
  
  //align_position
  if(align_position!=NULL)
    {
      for(unsigned i=0;i<original_numOfAligned;i++)
	{
	  if(align_position[i]!=NULL)
	    {
	      delete [] align_position[i];
	    }
	}
      delete [] align_position;
    }
  
  align_position=new unsigned*[numOfAligned];

  for(unsigned i=0;i<numOfAligned;i++)
    {
      align_position[i]=new unsigned[2];
      p_char=(char*)align_position[i];
      _ifs.read(p_char, sizeof(unsigned)*2);
    }
  
  //min_deletions
  if(min_deletions)
    delete [] min_deletions;
  min_deletions=new unsigned[numOfAligned];
  p_char=(char*)min_deletions;
  _ifs.read(p_char, sizeof(unsigned)*numOfAligned);
  
  //
  if(n_errors!=NULL)
    delete[] n_errors;
  n_errors=new unsigned[numOfAligned];
  p_char=(char*)n_errors;
  _ifs.read(p_char, sizeof(unsigned)*numOfAligned);

  //
  if(error_positions!=NULL)
    {
      for(unsigned i=0;i<original_numOfAligned;i++)
	{
	  if(error_positions[i]!=NULL)
	    {
	      delete [] error_positions[i];
	    }
	}
      delete [] error_positions;
    }
  
  error_positions=new unsigned*[numOfAligned];

  for(unsigned i=0;i<numOfAligned;i++)
    {
      error_positions[i]=new unsigned[n_errors[i]];
      p_char=(char*)error_positions[i];
      _ifs.read(p_char, sizeof(unsigned)*n_errors[i]);
    }
  
  //
  if(p_region_max_length!=NULL)
    {
      for(unsigned i=0;i<original_numOfAligned;i++)
	{
	  if(p_region_max_length[i]!=NULL)
	    {
	      delete [] p_region_max_length[i];
	    }
	}
      delete [] p_region_max_length;
    }
  
  p_region_max_length=new unsigned*[numOfAligned];

  for(unsigned i=0;i<numOfAligned;i++)
    {
      p_region_max_length[i]=new unsigned[1+maximum_deletion];
      p_char=(char*)p_region_max_length[i];
      _ifs.read(p_char, sizeof(unsigned)*(1+maximum_deletion));
    }

  //excess error
  if(excess_error_positions!=NULL)
    {
      for(unsigned i=0;i<original_numOfAligned;i++)
	{
	  if( excess_error_positions[i]!=NULL)
	    {
	      delete [] excess_error_positions[i];
	    }
	}
      delete [] excess_error_positions;
    }
  excess_error_positions=new unsigned*[numOfAligned];
  for(unsigned i=0;i<numOfAligned;i++)
    {
      excess_error_positions[i]=new unsigned[AlignmentSettings::negative_excess_deletions_max];
      p_char=(char*)excess_error_positions[i];
      _ifs.read(p_char, sizeof(unsigned)*AlignmentSettings::negative_excess_deletions_max);
    }

  //alleles_all
  if(alleles_all!=NULL)
    {
      delete[] alleles_all;
    }
  alleles_all=new unsigned[numOfAligned];
  p_char=(char*)alleles_all;
  _ifs.read(p_char, sizeof(unsigned)*numOfAligned);

  
  //
  if(alleles_from_distinct_genes!=NULL)
    {
      delete[] alleles_from_distinct_genes;
    }
  alleles_from_distinct_genes=new unsigned[numOfAligned];
  p_char=(char*)alleles_from_distinct_genes;
  _ifs.read(p_char, sizeof(unsigned)*numOfAligned);
  
  //done!!!
}

