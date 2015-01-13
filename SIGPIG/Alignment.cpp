#include "Alignment.hpp"

//enum MatchState{ MATCH, MISMATCH, NOT_ALIGNED};
Alignment_Object::Alignment_Object(): numberOfAligned(0), align_length(NULL),
  align_position(NULL), min_deletions(NULL), n_errors(NULL), error_positions(NULL),
  p_region_max_length(NULL), excess_error_positions(NULL), alleles_all(NULL),
  alleles_from_distinct_genes(NULL)
{
  //empty constructor
}

Alignment_Object::~Alignment_Object()
{
  //now we need to delete the data if there are 
  //we better do this in a reverse order than it is constructed
  
  //
  if(alleles_from_distinct_genes!=NULL)
    {
      delete[] alleles_from_distinct_genes;
    }
  //
  if(alleles_all!=NULL)
    {
      delete[] alleles_all;
    }
  //
  if(excess_error_positions!=NULL)
    {
      for(unsigned i=0;i<numberOfAligned;i++)
	{
	  if( excess_error_positions[i]!=NULL)
	    {
	      delete [] excess_error_positions[i];
	    }
	}
      delete [] excess_error_positions;
    }
  //
  if(p_region_max_length!=NULL)
    {
      for(unsigned i=0;i<numberOfAligned;i++)
	{
	  if(p_region_max_length[i]!=NULL)
	    {
	      delete [] p_region_max_length[i];
	    }
	}
      delete [] p_region_max_length;
    }
  //
  if(error_positions!=NULL)
    {
      for(unsigned i=0;i<numberOfAligned;i++)
	{
	  if(error_positions[i]!=NULL)
	    {
	      delete [] error_positions[i];
	    }
	}
      delete [] error_positions;
    }
  //
  if(n_errors!=NULL)
    delete[] n_errors;
  //
  if(min_deletions)
    delete [] min_deletions;
  //
  if(align_position!=NULL)
    {
      for(unsigned i=0;i<numberOfAligned;i++)
	{
	  if(align_positions[i]!=NULL)
	    {
	      delete [] align_positions[i];
	    }
	}
      delete [] align_positions;
    }
  // 
  if(align_length!=NULL)
    {
      delete[] align_length;
    }
  
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
//return: the align_length.
unsigned align_with_constraints_fixed_left_remove_right_errors
       (const string& _seq1, const string& _seq2, const unsigned& _maximum_errors, 
	const double& _error_cost,
	/*output*/ unsigned* _n_errors, unsinged* _error_positions)
{
  unsigned l_seq1=_seq1.size();
  unsigned l_seq2=_seq2.size();
  unsigned align_length=0;

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

      //check to see whether we have 
      if(&_n_errors>_maximum_errors)
	{
	  //we are reaching the maximum errors, now it is _maximum_errors+1, done!
	  break;
	}
      //this one is a mismatch, but not go beyond the maximum allowed ones, we remeber the position
      if(!matchFlag)
	{
	  _error_positions[*n_errors-1]=i;
	}
    }

  //now we are done, we need to know the align length
  //need to know whether it is because we go through the maximum_length or due to too many errors
  if(i>max_match_length)
    align_length=i; //because i starts at zero, so it is i-1+1
  else //due to maximum error
    align_length=i+1;//because it starts at zero, so it i+1;
  double best_score=align_length-_error_cost* (*_n_errors);
  
  
  //now we need to figure out best score along the longest match
  double running_score=best_score;
  
  
  for(unsigned j=0;j< *_n_errors ; j++)
    {
      running_score=(_error_positions[j]-1+1)-_error_cost*j;
      
      if(running_score>best_score)
	{
	  align_length=_error_positions[j]-1+1;
	  best_score=running_score;
	  *_n_errors=j;
	}
    }
  //prepare the error_position array based on the 
  //this is not necessary, we just need to according to the _n_errors to 
  //check the error_positions array elements
  for(unsigned j=*_n_errors;j<_maximum_errors;j++)
    {
      _error_positions[j]=-1;
    }

  //clean up
  //delete [] c;
  
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
unsigned align_with_constraints_fast_left(const string& _seq, const string& _target, 
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
  unsigned* current_error_positions=new unsigned [_max_errors];
  double current_score;
  unsigned max_length;
  //looping through the positions of sequence to force the left side of target to align with
  //subsequence of seq starting at the current index of _seq
  for(unsigned i=0;i<l_seq;i++)
    {
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
      string current_target=_target.substr(0,max_length);
      //call to do alignment forcing fixed left ends
      current_align_length=
	align_with_constraints_fixed_left_remove_right_errors(current_seq, current_target,
							      _maximum_errors,_error_cost,
							      &current_n_errors, current_error_positions);
      
      current_score=current=current_align_length-_error_cost*current_n_errors;
      
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
  
  align_positions[0]=best_match_index;
  align_positions[1]=0;
  
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
	      const GenomicV* _genJs, const unsigned& _numOfJSegs, 
	      const unsigned& _J_minimum_alignment_length, const unsigned& _J_maximum_deletion, 
	      const unsigned& _nagative_excess_deletion_max, const unsigned& _J_allowed_errors, 
	     const unsigned& _error_cost,
	     /*output*/ Alignment_Obj* _J
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

  //now we assume the return Aligment_obj has been correctly initialized. so we can use it now
  //set up pointers to the data and use them
  //***1) vector<unsigned>* align_length = (J->align_length); //% contains the length of the alignment

  //***2) align_position, 
  //% align_position(1,j) contains the first matched nt in seq (redundant, same as length(seq) - align_length +1 ),
  //% align_position_(2,j) is the first matched nt in genJ.seq{j} - also redundant information because we know where the sequencing primer ends.
  //vector<unsigned[2]> align_position = J->align_position;

  //***3)n_errors,  vector<unsigned> n_errors=J->n_errors;//zeros(length(genJ),1); % number of mismatches in alignment

  //***4)vector<vector<unsigned> >error_positions = J->error_positions ;//zeros(length(genJ),J_allowed_errors); % positions of mismatches
  //***5)vecotr<unsigned>  min_deletions = J->min_deletions; //zeros(length(genJ),1); % number of deletions implied by alignment	
  bool j_large_deletion_flag[numOfJSegs]; //=zeros(length(genJ),1); % Flag if deletions is too large
 unsigned l_seq = _seq.GetLength(); //length of the input sequence

 
 unsigned temp_align_position[2];
 unsigned* temp_n_errors=new unsigned[_numOfJSegs];
 unsigned* temp_error_positions=new unsigned[_J_allowed_errors];
 unsigned temp_align_length
 //% Loop through template alleles for alignment
 for(unsigned int i=0;i<_numOfJSegs;i++) //for j=1:length(genJ)
   {
     //%j=1
     //%disp(['loop: ' num2str(j)])
     //% Get highest scoring alignment for this allele, with acceptable number of errors
     //% NOTE: I use an alignment function that fixes position on the left (as in match_Vs) but I pass in the
     //% reversed sequences, because I want to fix the position on the right.
     unsigned l_target=_genJs[i].Get_Seq().GetLength();//length of the genomic template allele for current one

     //before doing the alignment, we need to flip/reverse the sequence, since the function do it assume doing
     //fixed on the left. 
     SequenceString rev_seq=FlipSequenceString(_seq);
     SequenceString rev_target=FlipSequenceString(_genJs[j]);

     //now calling to do the alignment
     temp_align_length=
       aligned_with_constraints_fast_left(rev_seq.GetSeq(), rev.target.GetSeq(), _J_allowed_errors, _error_cost,
					  temp_align_position, &temp_n_errors, temp_error_positions);
					  

     _J->align_length.push_back(temp_align_length);
     
     //[align_length(j), rev_align_position, n_errors(j), rev_error_positions]=align_with_constraints_fast_no_fix(seq(end:-1:1) , genJ(j).seq(end:-1:1) , J_allowed_errors, error_cost);
     
     _J->n_errors.push_back(temp_n_errors);
     // % Reverse the error positions to be relative to left-right orientation.
     //  %%%j  error positions is relative to the end of J chain
     vector<unsigned> temp_error_position_vec;
     if( temp_n_errors>0)
       {
	 //unsigned rev_j=temp_n_errors-1;
	 for(unsigned j=0;j<temp_n_errors;j++)
	   {
	     //unsigned tempValue=temp_error_positions[j]
	     temp_error_position_vec.push_back( temp_error_positions[temp_n_errors-j-1]) ;
	     //temp_error_positions[rev_j]=tempValue;
	     //rev_j--;
	   }
       }//end
     _J->error_positions.push_back(temp_error_position_vec);

     vector<unsigned> temp_align_position_vec;
	//% Calculate (redundant) align_position values. These are used in model scripts.
     temp_align_position_vec.push_back(
				  /*temp_align_position[0]*/ = (l_seq-(temp_align_position[0]+temp_align_length-1)-1 ) ); //for _seq first position
     temp_align_position_vec.push_back(
				       /*temp_align_position[1]*/= l_target-(temp_align_position[1]+temp_align_length-1)-1 );//for target second position
     _J->align_position.push_back(temp_align_position_vec);

      // % Calculate deletions implied by alignment
     _J->min_deletions.push_back( l_target-(temp_align_position[1]+temp_align_length-1) - 1);
      //    % Flag if number of deletions is too many
     j_large_deletion_flag[j]=false;
     if( _J->min_deletions.at(j) > _J_maximum_deletion)
	{
	  j_large_deletion_flag[j] = true;
	  //%continue;
        }//  end
      //    end

   }//end of for loop to go through the J segs alleles for alignment.
    //% Check!  No need to check in c++ code.
 //	    assert(all(min_deletions>=0));

 //	  % sort the genomic Js in descending order of 'score' and ascending order of min_deletions

 //calculate the score first
 double* scores=new double[numOfJSegs];
 //prepare the sorted index of the array.
 unsigned* sorted_index=new unsigned[numOfJSegs];
 for(unsigned k=0;k<numOfSegs;k++)
   {
     sorted_index[k]=k;
     scores[k]=_J->align_length[k]-error_cost*_J->n_errors[k];
   }
 
 QuickSort(scores, 0, numOfJSegs-1, sorted_index, _J->min_deletions);
 //	      scores  = align_length - error_cost*n_errors;
 //S = [-scores, min_deletions];
 //	  [~,order]=sortrows(S);
 //now we need to reverse the order, since the QuickSort is ascending, but for our purpose we need to descending.
 

	  % Set a score threshold for alleles.
	      min_score = max(scores) - 3*J_minimum_alignment_length;

	  % Subset of alleles that have high enough score, long enough alignment, and
	      % not too many deletions.
	      j_list = find( scores(order) >= min_score & align_length(order) >= J_minimum_alignment_length & j_large_deletion_flag(order)==0);

	  if isempty(j_list)
		      //		      %disp(['No reasonable J match for sequence' ])
		      ok_order = [];
	  seq_j_ok = false; % Sequence not ok!

	  else
				ok_order = order(j_list);
	  seq_j_ok = true; % Sequence ok!
			       end

			       if seq_j_ok

				     % Store all the alignment information in J for acceptable alleles
												 J.align_length = align_length(ok_order);
	  J.align_position =  align_position(:,ok_order);
	  J.min_deletions = min_deletions(ok_order);
	  J.n_errors = n_errors(ok_order);
	  J.error_positions = error_positions(ok_order,:);
	  J.p_region_max_length = zeros(numel(ok_order),J_maximum_deletion+1);
	  J.excess_error_positions = zeros(numel(ok_order),negative_excess_deletions_max);
	  J.alleles_all = ok_order;

	      % Find gene indices of the acceptable alleles.
		      %%%codegen compatible version of
		  %%%genJ_ok_gene_inds = [genJ(ok_order).gene_index];
	  genJ_ok = genJ(ok_order);
	  genJ_ok_gene_inds = zeros(1,numel(ok_order));
	  for g=1:numel(ok_order)
		  genJ_ok_gene_inds(g) = genJ_ok(g).gene_index;
	      end

		    % Get list of distinct genes from acceptable alleles.
		[~,Jg] = unique(genJ_ok_gene_inds,"first" ); % Pick the first allele listed for each gene.
												       Jg = sort(Jg); % sort it because unique returns in ascending order of ok_order. This puts it back in order of best match.
															  % Only one (best) allele from each gene
															  J.alleles_from_distinct_genes = ok_order(Jg'); '

//    %Now go through and find palindromic nucleotides for various number of deletions,
//    %as well as error positions for 'negative' deletions.
    for j=1:numel(ok_order)
        % Loop over number of deletions (nd is actual number of genomic deletions).
        % Lower bound is maximum of 0 and min_deletions - negative_excess_deletions_max (usually 3).
        % Upper bound is minimum of J_maximum_deletion (usually 12) and the whole alignment being random insertions
        for nd=(max(0,J.min_deletions(j) - negative_excess_deletions_max)):min(J_maximum_deletion, J.align_length(j) + J.min_deletions(j))

            % For each value of deletions, find longest half-palindrome
            % from the implied end of the gene sequence
            p=0;
            still_palindrome=true;
            % Upper bound on p is the length of the J match (accounting for deletions) OR the never actually possible case of sequence length to the left of J match (accounting for deletions)
            while still_palindrome && p < min( J.align_length(j) - (nd-J.min_deletions(j)), J.align_position(1,j)-1 + nd - J.min_deletions(j))
                still_palindrome = genJ(ok_order(j)).seq(J.align_position(2,j) + p + (nd - J.min_deletions(j) )) == dnacomplement(seq(J.align_position(1,j) - p + nd - J.min_deletions(j) - 1) );
                if still_palindrome
                    p = p + 1;
                end
            end

            % Store length of longest half-palindrome for each value of
            % deletions.
            J.p_region_max_length(j, 1+nd) = p;
        end
   //        % Now for the 'negative' deletions, store where the mismatches are.
        % This is so we can count number of errors for these possibilities.

  //        % Number of 'negative' deletions to consider:
        n_excess = min([negative_excess_deletions_max, J.min_deletions(j), J.align_position(1,j)-1 ]);

        % Find mismatches between sequence and genomic J at those positions.
        excess_err_pos = J.align_position(1,j) - ( n_excess + 1 - find((genJ(ok_order(j)).seq((J.align_position(2,j)-n_excess): (J.align_position(2,j)-1) )~=seq((J.align_position(1,j) -n_excess): (J.align_position(1,j)-1) ))) );
        J.excess_error_positions(j, 1:length(excess_err_pos)) = excess_err_pos;
    end

else

    % Since sequence is not ok, return empty values.
    J.align_length = [];
    J.align_position =  [];
    J.min_deletions = [];
    J.n_errors = [];
    J.error_positions = [];
    J.p_region_max_length = [];
    J.excess_error_positions = [];
    J.alleles_all = [];
    J.alleles_from_distinct_genes = [];
end

end
  

//clean up
delete [] temp_error_positions;

  return true;
}




