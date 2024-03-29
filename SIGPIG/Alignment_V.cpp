
#include "Alignment.hpp"
#include "Alignment_V.hpp"
#include "../string_ext.hpp"



void DeterminePalindromAndExcessError_V
( const SequenceString& _seq, const GenomicV* _genVs, const unsigned* _ok_order,
  const unsigned* _min_deletions, 
  const unsigned& _negative_excess_deletions_max, const unsigned& _V_maximum_deletion,
  const unsigned* _align_length, const unsigned& _numOfAligned, unsigned** _align_positions,
  /*output*/ unsigned** _p_region_max_length, unsigned** _excess_error_position  
  )
{
  unsigned p=0;
  //cout<<"in the beginning p_Region_max_length:"<<endl;
  /*for(unsigned i=0;i<_numOfAligned;i++)
   {
     cout<<"\t"<<endl;
     for(unsigned j=0;j<_V_maximum_deletion+1;j++)
       {
	 cout<<_p_region_max_length[i][j]<<",";
       }
     cout<<endl;
     }*/
  bool still_palindrome=true;
  //go through and find palindromic nucleotides for various number of deletions,
  //    %as well as error positions for 'negative' deletions.
  //cout<<"\tinside palindro.....before loop. _numOfAligned:"<<_numOfAligned<<endl;
  for(unsigned  j=0;j<_numOfAligned;j++)
    {
      //cout<<"\t*^^^^^^loop "<<j<<endl;
      //cout<<"\t*^^^^^^seq:  "<<_seq.toString()<<endl;
      
      string target=_genVs[_ok_order[j]].Get_Sequence();
      //cout<<"\t*^^^^^^target:  "<<target<<endl;
      //cout<<"\talign_length:"<<_align_length[j]<<endl;
      //cout<<"\talign_position:"<<_align_positions[j][0]<<","<<_align_positions[j][1]<<endl;
      //cout<<"\t_min_deletion:"<<_min_deletions[j]<<endl;

      //% Loop over number of deletions (nd is actual number of genomic deletions).
      //  % Lower bound is maximum of 0 and min_deletions - negative_excess_deletions_max (usually 3).
      //  % Upper bound is minimum of J_maximum_deletion (usually 12) and the whole alignment being random insertions
      int nd=_min_deletions[j]-_negative_excess_deletions_max;
      if(nd<0)
	nd=0;
      int max_nd=_V_maximum_deletion;
      if(max_nd>(signed)(_align_length[j]+_min_deletions[j]))
	max_nd=_align_length[j]+_min_deletions[j];
      //cout<<"\t\t^^^^^second loop before, nd:"<<nd<<endl;
      for(;nd<=max_nd;nd++)
	{
	  //cout<<"\t\t\tsub loop:"<<nd<<endl;
	  //% For each value of deletions, find longest half-palindrome
          //  % from the implied end of the gene sequence
	  p=0;
	  still_palindrome=true;
	  //% Upper bound on p is the length of the J match (accounting for deletions) OR the never actually possible case of sequence length to the left of J match (accounting for deletions)
	  int max_p_length=_align_length[j]-(nd-_min_deletions[j]);
	  
	  if(max_p_length>(signed)(_seq.GetLength()-(_align_positions[j][0])-_align_length[j]+nd-_min_deletions[j]))
	    max_p_length=(signed)(_seq.GetLength()-(_align_positions[j][0])-_align_length[j]+nd-_min_deletions[j]);

	  if(max_p_length<0)
	    max_p_length=0;
	  //cout<<"\t\t\tchecking while loop:max_p_length:"<<max_p_length<<endl;
	  while( still_palindrome && (signed)p < max_p_length)
	    {	  
	      //cout<<"\t\t\tstarting while......."<<endl;
	      //cout<<"target at:"<<_align_positions[j][1] +_align_length[j] -1 - p - (nd - _min_deletions[j] )-0<<endl;
	      //cout<<"seq at:"<<_align_positions[j][0]+ _align_length[j]-1 -(nd - _min_deletions[j])+p +1<<endl;
	      still_palindrome = target.at(_align_positions[j][1] +_align_length[j] -1 - p - (nd - _min_deletions[j] )-0) == DnaComplement(_seq.GetSequence().at(_align_positions[j][0]+ _align_length[j]-1 -(nd - _min_deletions[j])+p +1) );
	      //cout<<"\t\t\tafter while.........."<<endl;
	      if (still_palindrome)
		{
		  p ++;
		}
	    }
	  //cout<<"\t\t\tdone!!!"<<endl;
	  //% Store length of longest half-palindrome for each value of
	  //% deletions.
	  _p_region_max_length[j][ nd] = p;
	}
      //cout<<"\t\tend of second loop"<<endl;
      //        % Now for the 'negative' deletions, store where the mismatches are.
      // % This is so we can count number of errors for these possibilities.
      
      //        % Number of 'negative' deletions to consider:
      unsigned tempArry[]={_negative_excess_deletions_max, _min_deletions[j], _seq.GetLength()-(_align_positions[j][0])-_align_length[j] };
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
	  if(target.at(_align_positions[j][1]+_align_length[j]+k)!=_seq.GetSequence().at(_align_positions[j][0]+_align_length[j]+k))
	    {
	      //cout<<"\t\tfound one!!pos:"<<_align_positions[j][0]-n_excess+k<<endl;
	      _excess_error_position[j][runningIndex_excessError]=/*_align_positions[j][0]+*/_align_length[j]+k;
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
    }//end of outer for loop for all the aligned strings
  //cout<<"Done for the function"<<endl;
  //return true;
}


//defining the functions doing alignment

//the function to do the V alignments
//input:
// _seq, the input seq that is being aligned against. constant reference 
//_genVs, the V gen segments arrays.
//_numOfVSegs, the number of v gene segments in the array
//_V_minimum_alignment_length, length of v mininum alignement
//_V_maximum_deletion, maximum deletion
//_negative_excess_deletion_max, maximum excess deletion, usually 3
//_v_allowed_errors, number maximum allowed errors in the alignment
//_error_cost, mismatch error cost, -4 or -5 (?)
//
//output:
// _V, the Alignment_obj pointer V holding the alignment details. The caller need to allocate the memory
//bool, indicating whether the alignment is successful.
bool match_V(const SequenceString& _seq,
	     const GenomicV* _genVs, const unsigned& _numOfVSegs, 
	     const unsigned& _V_minimum_alignment_length,
	     const unsigned& _V_maximum_deletion, 
	     const unsigned& _negative_excess_deletion_max,
	     const unsigned& _V_allowed_errors, 
	     const unsigned& _error_cost,
	     /*output*/ Alignment_Object& _V
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
  //***1) unsigned* align_length, size numOfVSegs maximumly ; ---temp_align_length//% contains the length of the alignment

  //***2) align_position, 2D, numOfVSegs by 2
  //% align_position(1,j) contains the first matched nt in seq (redundant, same as length(seq) - align_length +1 ),
  //% align_position_(2,j) is the first matched nt in genJ.seq{j} - also redundant information because we know where the sequencing primer ends.
  // temp_align_position = J->align_position;

  //***3) min_deletions, 1D, maximumly numOfVSegs long, temp_min_deletions;

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
  cout<<"\t_numOfVSegs:"<<_numOfVSegs
      <<"\n\t_V_minmum_alignment_length:"<<_V_minimum_alignment_length
      <<"\n\t_V_maximum_deletion:"<<_V_maximum_deletion
      <<"\n\t_negative_excess_deletion_max:"<<_negative_excess_deletion_max
      <<"\n\t_V_allowed_errors"<<_V_allowed_errors
      <<"\n\t_error_cost"<<_error_cost<<endl;*/

  bool* v_large_deletion_flag=new bool [_numOfVSegs]; //=zeros(length(genJ),1); % Flag if deletions is too large
  //unsigned l_seq = _seq.GetLength(); //length of the input sequence
  //cout<<"l_seq:"<<l_seq<<endl;
  unsigned* temp_align_length=new unsigned[_numOfVSegs];
  unsigned** temp_align_position=new unsigned*[_numOfVSegs];
  for(unsigned i=0;i<_numOfVSegs;i++)
    {
      temp_align_position[i]=new unsigned[2];
    }
  unsigned* temp_min_deletions=new unsigned[_numOfVSegs];

  unsigned* temp_n_errors=new unsigned[_numOfVSegs];
  //std::memset(temp_n_errors, -1, _numOfVSegs);//fill the default value with -1
  unsigned** temp_error_positions=new unsigned* [_numOfVSegs];
  for(unsigned i=0;i<_numOfVSegs;i++)
    {
      temp_error_positions[i]=new unsigned[_V_allowed_errors];
    }
  //unsigned** temp_p_region_max_length=NULL;
  //unsigned** temp_excess_error_positions=NULL;
  //unsigned* temp_alleles_all; unsigned* alleles_from_distinct_gene;
  //temp_alleles_all=NULL; alleles_from_distinct_gene=NULL;

  //cout<<"***f2irst****"<<endl;
  
  //the following are only for setting up the output for the 
  //unsigned align_length_func;
  //unsigned align_position_func[2];
  //unsigned* error_position_func=new unsigned[_V_allowed_errors];
  //SequenceString rev_seq=FlipSequanceString(_seq);
  //unsigned n_errors_func;
  //% Loop through template alleles for alignment
  //cout<<"***first****3"<<endl;
  for(unsigned int i=0;i<_numOfVSegs;i++) //for j=1:length(genJ)
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
     
     SequenceString target=_genVs[i].Get_Seq();
     unsigned l_target =target.GetLength();
     //cout<<"\t&&&doing alignment :"<<endl;
     //cout<<"\t_seq:"<<_seq.toString()<<endl;
     //cout<<"\ttarget:"<<target.toString()<<endl;
     //now calling to do the alignment
     temp_align_length[i]= align_with_constraints_fast_no_fix(_seq.GetSequence(), target.GetSequence(), _V_allowed_errors, _V_minimum_alignment_length,_error_cost,
					temp_align_position[i], temp_n_errors[i], temp_error_positions[i]);
					  
     //cout<<"\ttemp_align_length["<<i<<"]:"<<temp_align_length[i]<<endl;
     //cout<<"\ttemp_n_errors[i]"<<i<<"]:"<<temp_n_errors[i]<<endl;
     //cout<<"\talign_position_func"<<temp_align_position[i][0]<<","<<temp_align_position[i][1]<<endl;
     //cout<<"\tscore"<<CalculateScore(temp_align_length[i], temp_error_positions[i], temp_n_errors[i], _error_cost);
     
     //_J->align_length.push_back(temp_align_length);
     
     //[align_length(j), rev_align_position, n_errors(j), rev_error_positions]=align_with_constraints_fast_no_fix(seq(end:-1:1) , genJ(j).seq(end:-1:1) , J_allowed_errors, error_cost);
     
     //_J->n_errors.push_back(temp_n_errors);
     /*
     // % Reverse the error positions to be relative to left-right orientation.
     // ********be careful here, we decide to keep the j error position relative to the END OF J CHAIN
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
        = (l_seq-(align_position_func[0]+temp_align_length[i]-1)-1 ) ; //for _seq first position
     temp_align_position[i][1]
       = l_target-(align_position_func[1]+temp_align_length[i]-1)-1 ;//for target second position
     //_J->align_position.push_back(temp_align_position_vec);
     */
     
     // % Calculate deletions implied by alignment
     temp_min_deletions[i]=l_target-(temp_align_position[i][1]+temp_align_length[i]-1) - 1;
     //cout<<"\tmin_deletion:"<<temp_min_deletions[i]<<endl;
     //    % Flag if number of deletions is too many
     v_large_deletion_flag[i]=false;
     if( temp_min_deletions[i] > _V_maximum_deletion)
	{
	  v_large_deletion_flag[i] = true;
	  //%continue;
        }//  end
      //    end

   }//end of for loop to go through the J segs alleles for alignment.
    //% Check!  No need to check in c++ code.
 //	    assert(all(min_deletions>=0));

 //	  % sort the genomic Js in descending order of 'score' and ascending order of min_deletions

 //calculate the score first
 // cout<<"***first****3"<<endl;
 double* scores=new double[_numOfVSegs];
 //prepare the sorted index of the array.
 unsigned* sorted_index=new unsigned[_numOfVSegs];
 //cout<<"=========>before sorting:";
 for(unsigned k=0;k<_numOfVSegs;k++)
   {
     sorted_index[k]=k;
     scores[k]=CalculateScore(temp_align_length[k], temp_error_positions[k], temp_n_errors[k], _error_cost);//temp_align_length[k]-_error_cost*temp_n_errors[k];
     //cout<<scores[k]<<"-"<<sorted_index[k]<<"-"<<temp_min_deletions[k]<<"-"<<temp_align_length[k]<<",";
     
   }
 //cout<<endl;
 // cout<<"\n***first****3aa"<<endl;
 //sorted index also holding the gene allele index
 QuickSort<double>(scores, 0, _numOfVSegs-1, sorted_index, temp_min_deletions);
 //	      scores  = align_length - error_cost*n_errors;
 //S = [-scores, min_deletions];
 //	  [~,order]=sortrows(S);
 //cout<<"=========>after sorting:";
 /*for(unsigned k=0;k<_numOfVSegs;k++)
   {
     //sorted_index[k]=k;
     //scores[k]=temp_align_length[k]-_error_cost*temp_n_errors[k];
     cout<<k<<":"<<scores[k]<<"-"<<sorted_index[k]<<"-"<<temp_min_deletions[k]<<"-"<<temp_align_length[sorted_index[k]]<<",";
   }
   cout<<endl;*/
 //cout<<"\n***first****3bb"<<endl;
 //now we need to reverse the order, since the QuickSort is ascending, but for our purpose we need to descending.
 Reverse(sorted_index, _numOfVSegs);
 Reverse(scores, _numOfVSegs);
 Reverse(temp_min_deletions, _numOfVSegs);
 /*
 cout<<"=========>after sorting:";
 for(unsigned k=0;k<_numOfVSegs;k++)
   {
     //sorted_index[k]=k;
     //scores[k]=temp_align_length[k]-_error_cost*temp_n_errors[k];
     cout<<scores[k]<<"-"<<sorted_index[k]<<",";
     }*/
 //cout<<"\n***first****3ccc"<<endl;
 //  % Set a score threshold for alleles. well, this is kind of arbitrary
 //we want to get the best ones, but limited numbers 
 double min_score=max_mf(scores,_numOfVSegs)-3*_V_minimum_alignment_length;
 unsigned max_error=temp_n_errors[sorted_index[0]]+3;//another arbitrary threshold
 //cout<<"min_score:"<<min_score<<endl;
 //% Subset of alleles that have high enough score, long enough alignment, and
 //	      % not too many deletions.
 unsigned* ok_order = new unsigned[_numOfVSegs]; //this will directly used by J.alleles_all, so do NOT delete/clean later
 unsigned ok_count=0;
 //cout<<"***first****4"<<endl;
 for(unsigned i=0;i<_numOfVSegs;i++)
   { 

     //cout<<"index i:"<<i<<"----****"<<ok_count<<endl;
     
     //find( scores(order) >= min_score & align_length(order) >= J_minimum_alignment_length & v_large_deletion_flag(order)==0);
     if(scores[i]>=min_score&&temp_align_length[sorted_index[i]]>=_V_minimum_alignment_length&&!v_large_deletion_flag[sorted_index[i]]
	&&temp_n_errors[sorted_index[i]]<=max_error)
       {
	 ok_order[ok_count]=sorted_index[i];
	 temp_min_deletions[ok_count]=temp_min_deletions[i];
	 ok_count++;
       }

   }
 //cout<<"---------->showing the ok_order array:"<<ok_count<<endl;
 /*for(unsigned i=0;i<ok_count;i++)
   {
   cout<<i<<"-"<<ok_order[i]<<",";
   }
 cout<<endl;
 */
 //bool seq_j_ok=true;
 //cout<<"*******first 4a check for status"<<endl;
 if(ok_count<=0)//empty
   {
     //cout<<"*****inside false condition"<<endl;
     _V.numOfAligned=0;
     _V.maximum_deletion=0;
     //do we want to clean up the memory, not necessary???
     //we have to clean up the memory, by now some of the arrays have been allocated
     //<--------
     //clean up
     //cout<<"\t555delete 1"<<endl;
     delete [] v_large_deletion_flag;
     v_large_deletion_flag=NULL;
     //cout<<"\t555delete 2"<<endl;
     delete [] temp_align_length;
     temp_align_length=NULL;
     //cout<<"\t555delete 3"<<endl;
     CleanUpMemory(temp_align_position, _numOfVSegs);
     //cout<<"\t555delete 4"<<endl;
     delete [] temp_min_deletions;
     temp_min_deletions=NULL;
     delete [] temp_n_errors;
     temp_n_errors=NULL;
     CleanUpMemory(temp_error_positions, _numOfVSegs);
 //CleanUpMemory(temp_p_region_max_length, _numOfVSegs);
 //CleanUpMemory(temp_excess_error_positions, _numOfVsegs);
 //delete [] temp_alleles_all;
 //delete [] alleles_from_distinc_gene;
     //cout<<"\t555delete 5"<<endl;
     //delete [] error_position_func;
     //cout<<"\t555delete 6"<<endl;
     delete [] scores;
     scores=NULL;
     //cout<<"\t555delete 7"<<endl;
     delete [] sorted_index;
     sorted_index=NULL;
     //cout<<"\t555delete 8"<<endl;
     delete [] ok_order; //ok order has not been referenced directly by J.alleles_all
     ok_order=NULL;
     //and in this case, it will not. so we delete it.
 
     //delete [] genJ_ok_index;
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
 _V.numOfAligned=ok_count;
 _V.maximum_deletion=_V_maximum_deletion;
 //cout<<"***first****4aaa"<<endl;
 _V.align_length = new unsigned [ok_count];
 //cout<<"****ok_count***:"<<ok_count<<endl;
 if(!CopyElements(temp_align_length, _numOfVSegs,  _V.align_length, ok_count, ok_order, ok_count))
   return false;
 //cout<<"***first****4bbbb"<<endl;
 _V.align_position =new unsigned* [ok_count];
 for(unsigned m=0;m<ok_count;m++)
   {
     _V.align_position[m]=new unsigned [2];
   }
 if(!CopyElements(temp_align_position, _numOfVSegs, 2, _V.align_position, ok_count, 2, ok_order, ok_count))//(:,ok_order);
   {
     return false;
   }
 //cout<<"***first****4cccc"<<endl;
 _V.min_deletions =new unsigned [ok_count];
 //cout<<"\t***&&&&&&showing min deletions"<<endl;
 //cout<<"\t\t";
 /*for(unsigned _m=0;_m<_numOfVSegs;_m++)
   {
     cout<<temp_min_deletions[_m]<<","<<endl;
   }
   cout<<endl;*/
 //****NOTE:important, temp_min_deletions, is different than others, it is sorted together with socres, so we only need to copy over
 //what is in the front (best values) according to ok_count
 for(unsigned _m=0;_m<ok_count;_m++)
   {
     _V.min_deletions[_m]=temp_min_deletions[_m];//<<","<<endl;
   }
 //cout<<"***first****4dddd"<<endl;
 _V.n_errors = new unsigned [ok_count];
 if(!CopyElements(temp_n_errors, _numOfVSegs, _V.n_errors, ok_count, ok_order, ok_count))
   return false;
 //cout<<"***first****4eee"<<endl;
 _V.error_positions = new unsigned* [ok_count];
 for(unsigned m=0;m<ok_count;m++)
   {
     _V.error_positions[m]=new unsigned [_V_allowed_errors];
   }
 //cout<<"***first****4fff"<<endl;
 if(!CopyElements(temp_error_positions, _numOfVSegs, _V_allowed_errors, _V.error_positions, ok_count, _V_allowed_errors, ok_order, ok_count))
   return false;//error_positions(ok_order,:);

 //now we are done so, and need to run initialization for the later assignments
 _V.p_region_max_length = new unsigned* [ok_count];

 //cout<<"***first5****"<<endl;
 for(unsigned i=0;i<ok_count;i++)
   {
     _V.p_region_max_length[i]=new unsigned[_V_maximum_deletion+1];
     std::memset(_V.p_region_max_length[i], 0, (_V_maximum_deletion+1)*sizeof(unsigned)/sizeof(char));//fill the default value with 0
   }
 /*cout<<"right after setting the memory:maxlength:"<<_V_maximum_deletion<<endl;
for(unsigned i=0;i<ok_count;i++)
   {
     cout<<"\t"<<endl;
     for(unsigned j=0;j<_V_maximum_deletion+1;j++)
       {
	 cout<<_V.p_region_max_length[i][j]<<",";
       }
     cout<<endl;
     }*/
 //zeros(numel(ok_order),J_maximum_deletion+1);
 _V.excess_error_positions =new unsigned* [ok_count];
 
 for(unsigned i=0;i<ok_count;i++)
   {
     _V.excess_error_positions[i]=new unsigned [_negative_excess_deletion_max];
   }
 _V.alleles_all = ok_order;
 
 // % Find gene indices of the acceptable alleles.
 //		      %%%codegen compatible version of
 //		  %%%genJ_ok_gene_inds = [genJ(ok_order).gene_index];
 unsigned* genV_ok_index = new unsigned[ok_count];
 //go through the genJ to figure out the distinct genes for the each allele, and then reture a array of indices to them
 
 //cout<<"\t^^^^^showing the gene index of the alleles:";
 for(unsigned i=0;i<ok_count;i++)
   {
     genV_ok_index[i]=_genVs[_V.alleles_all[i]].Get_GeneIndex();
     //cout<<genV_ok_index[i]<<",";
   }
 //cout<<endl;
 //genJ_ok_gene_inds = zeros(1,numel(ok_order));
//first we need to sort the J.alleles_all array 
 unsigned* sorted_genV_ok_index_temp=new unsigned [ok_count];
 unsigned* sorted_genV_ok_index_temp_index=new unsigned[ok_count];
 unsigned numOfUnique;
 Unique(genV_ok_index, ok_count, sorted_genV_ok_index_temp, sorted_genV_ok_index_temp_index,numOfUnique);
 //cout<<"after unique:numOfUnique:"<<numOfUnique<<endl;
 //now we got the distinct gene index index, just need to sort it, to make it in order
 //because when we do the unique, we first sort it and the gene index might not in order
 //so the index of the gene index is not in order either
 QuickSort(sorted_genV_ok_index_temp_index,0, numOfUnique-1);
 //cout<<"after sorting...."<<endl;
 //cout<<"the index index:"<< sorted_genV_ok_index_temp_index[0]<<endl;

 //cout<<"***first****6"<<endl;
 //compare and then get it copied to the alignment object
 _V.alleles_from_distinct_genes=new unsigned[ok_count];

 //unsigned runningValue=genJ_ok_index[0];
 //unsigned runningIndexOfAlleleArray=1;
 //_J.alleles_from_distinct_genes[0]=_J.alleles_all[0];
 for(unsigned i=0;i<numOfUnique;i++)
   {
     //if(runningValue!=genJ_ok_index[i])
     //  {
	 _V.alleles_from_distinct_genes[i] = ok_order[sorted_genV_ok_index_temp_index[i]];
	 //runningValue=genJ_ok_index[i];
	 //runningIndexOfAlleleArray++;
	 //}
   }
 //need to fill the rest of the distinc gene array with -1 as an indicator
 for(unsigned i=numOfUnique;i<ok_count;i++)
   {
     _V.alleles_from_distinct_genes[i]=-1;
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
 DeterminePalindromAndExcessError_V
   ( _seq, _genVs, ok_order,_V.min_deletions, 
     _negative_excess_deletion_max, _V_maximum_deletion,
     _V.align_length, _V.numOfAligned, _V.align_position,
     _V.p_region_max_length, _V.excess_error_positions  
     );

  
 //cout<<"***first****7"<<endl;
 //clean up
 delete [] v_large_deletion_flag;
 
 delete [] temp_align_length;
 CleanUpMemory(temp_align_position, _numOfVSegs);
 delete [] temp_min_deletions;
 delete [] temp_n_errors;
 CleanUpMemory(temp_error_positions, _numOfVSegs);
 //CleanUpMemory(temp_p_region_max_length, _numOfJSegs);
 //CleanUpMemory(temp_excess_error_positions, _numOfJsegs);
 //delete [] temp_alleles_all;
 //delete [] alleles_from_distinc_gene;
 
 // delete [] error_position_func;

 delete [] scores;
 delete [] sorted_index;
 //delete [] ok_order; ok order is referenced directly by J.alleles_all
 
 delete [] genV_ok_index;
 delete [] sorted_genV_ok_index_temp;
 delete [] sorted_genV_ok_index_temp_index;
 return true;
}

unsigned align_with_constraints_fast_no_fix
     (const string& _seq, const string& _target, 
      const unsigned& _maximum_errors, const unsigned& _minimum_align_length,  
      const double& _cost,
      /*output*/unsigned* _align_position, unsigned& _n_errors, 
      unsigned* _error_positions)
{
  unsigned align_length=0;
    
  unsigned l_seq=_seq.size();
  unsigned l_target=_target.size();
  int best_match_index=0;
  double score=-1E200;

  unsigned current_align_length, current_n_errors;
  unsigned* current_error_positions=new unsigned [_maximum_errors];
  double current_score;
  unsigned max_length;
  unsigned target_start=0, seq_start=0;

  string current_seq, current_target;
  unsigned align_position_start=0;
  unsigned curr_l_seq, curr_l_target;
  cout<<"param:maximum_errors:"<<_maximum_errors<<";cost:"<<_cost
      <<"minimum_align_length:"<<_minimum_align_length
      <<endl;
  //looping through the positions of sequence not to force  the left side  to align 
  //here start at -l_target+1 since they have to be overlap otherwise no way to
  //align. Same here for l_seq-1, we need them to be overlap at least with one nucleotide
  //the total aligned position is l_seq+l_target-2.
  for(signed int i=-1*(l_target-1);i<=((signed)(l_seq-1));i++)  //this i is the length, not index; so the index =i-1?
    {
      cout<<"\t4*************doing loop "<<i<<endl;
      target_start=0;
      seq_start=0;

      curr_l_seq=l_seq;
      curr_l_target=l_target;

      if(i<0) //in this case, we are not aligning the whole length of target,
	//I mean, starting from the middle of target
	{
	  target_start=-1*i;
	  curr_l_target=l_target-target_start+1;
	}
      else //i>=0, we are starting from beginning of target, but in the middle of seq
	{
	  seq_start=i;
	  curr_l_seq=l_seq-seq_start+1;
	}
      
      max_length=curr_l_seq;
      if(curr_l_target<max_length)
	{
	  max_length=curr_l_target;
	}
      cout<<"\t4 max_length:"<<max_length<<": seq_start:"<<seq_start<<"; target start:"<<target_start<<endl;
      //can not be better, this is for i>0, in the middle of seq with whole target case
      if(i>0&&max_length<score)
	{
	  break;
	}
      if(i>0&&max_length< _minimum_align_length)
	{
	  cout<<"===>breaking........"<<endl;
	  break;
	}
      if(i<0&&max_length< _minimum_align_length)
	{
	  cout<<"--->continue"<<endl;
	  continue;//don't do this one.
	}
      current_seq=_seq.substr(seq_start, max_length);
      current_target=_target.substr(target_start, max_length);
      unsigned current_align_position_start;//this is the one relative to the begining of input sequences
      //cout<<"\tseq:"<<current_seq<<endl;
      //cout<<"\ttarget:"<<current_target<<endl;
      //call to do alignment forcing fixed left ends
      current_align_length= align_with_constraints_fixed_left_remove_both_errors
	(current_seq, current_target, _maximum_errors, _cost,
	 /*output*/ current_n_errors, current_error_positions, current_align_position_start);
      //***---> NOTE: the current_align_position_start and current_error_positions are all relative to the begining of the current sub seq and target 
      //cout<<"\t\t===>right after calling removeBoth error"<<endl;
      if(current_align_length==0)
	continue;
      
      current_score=CalculateScore(current_align_length, current_error_positions, current_n_errors, _cost);
      cout<<"\tcurrent_score:"<<current_score<<";current_align_length:"<<current_align_length<<endl;
      //cout<<"\tbest score so far:"<<score<<endl;
      //we got a better one or got an identical score, but long, we are good
      if(current_score>score||
	 ((current_score-score<1E-6&&current_score-score>-1E-6)&&(current_align_length>align_length))
	 )
	{
	  align_length=current_align_length;
	  _n_errors=current_n_errors;
	  //copy over the corrent error positions
	  for(unsigned j=0;j<current_n_errors;j++)
	    {
	      _error_positions[j]=current_error_positions[j];
	    }
	  score=current_score;
	  best_match_index=i;
	  align_position_start=current_align_position_start;
	}//if loop
      
    }//end of looping over the seq positions
  //cout<<"best_match_index:"<<best_match_index<<endl;
  //now to prepare the output
  if(best_match_index<0)  //starting in the middle of target
    {
      target_start=-1*best_match_index;
      seq_start=0;
    }
  else //alignment starts in the middle of seq
    {
      target_start=0;
      seq_start=best_match_index;
    }
  
  _align_position[0]=seq_start+align_position_start;
  _align_position[1]=target_start+align_position_start;

  //error number and error_positions keeps as the return from function "...remove_both_errrors()"
  
  //clean up
  delete [] current_error_positions;
  return align_length;
}

/*the function wrapper for calculate the score of 
 *the alignment.
 *so far the score is simple, but we might want to change
 *later. so make it easy for now.
 * return: the score
 */ 
double CalculateScore
  (const unsigned& _total_aligned_len, 
   const unsigned* _error_positions,
   const unsigned& _n_errors,
   const double& _cost)
{
  //make the score positive anyway.
  double temp_cost=_cost;
  if(_cost<0)
    temp_cost =-1*_cost;
  return _total_aligned_len-temp_cost*_n_errors;
}

/*the function used to figure out the sub_alignment with the best
 *score along the full alignment. Here we add one more parameter
 *allow the subalignment could move freely on the left side, 
 *"_curr_error_position_index" (see below for explanation)
 *
 *input:
 *   _total_align_len:the full length of the aligment we 
 *            are working to get the best sub alignment 
 *   _all_error_positions: the error positoins of the full 
 *            alignment. The positions are relative to the
 *            beginning of the alignment.
 *   _curr_error_position_index: we add this one to make the 
 *            best alignment could start not necessarily with 
 *            the fixed left end. The outer caller can set up 
 *            to make the search from the left, ie varying
 *            left to make it removing left errors.
 *            THIS VALUE COULD BE -1, which
 *            means we started from very beginning of the full
 *            alignment and include all the errors. If it is 0, 
 *            it means we don't include first error and start
 *            to consider the sequence beginning from 
 *            error_position[0]+1 positoins.!!!!!
 *output:
 *   _best_n_errors: the n_errors give the best score. it starts
 *            right after the curr_error_position_index (DON'T
 *            include this one).
 *   NOTE: we don't include the error positoins, since we can
 *         figure it out by other information. 
 * return: the best aligned length 
 */
unsigned  GetBestScoreAmongTheAlignment
   (const unsigned& _total_align_len, 
    const unsigned* _all_error_positions,
    const unsigned& _all_n_errors,  
    const unsigned& _curr_error_position_index,
    const double& _cost,
    /*output*/ 
    unsigned& _best_n_errors  )
{

  //_best_n_errors has been defined in the function definition
  unsigned best_align_len;
  double best_score;
  //double best_n_errors;

  //
  unsigned temp_align_len;
  double temp_score;
  unsigned temp_n_errors;
  //get the longest one as the best to start with
  //check to see whether are doing the special case, -1
  if(((signed int)_curr_error_position_index)==-1)
    {
      _best_n_errors=_all_n_errors;
      best_align_len=_total_align_len;
      best_score=CalculateScore(best_align_len, _all_error_positions,
				_best_n_errors, _cost);	

      
    }
  else
    {
      _best_n_errors=_all_n_errors-(_curr_error_position_index+1);
      best_align_len=
      (_total_align_len-1)-_all_error_positions[_curr_error_position_index];
      best_score=
	CalculateScore(best_align_len, 
	      _all_error_positions+_curr_error_position_index+1, 
		       _best_n_errors, _cost);
    }
  //cout<<"-->start sub looping"<<endl;
  //now we have the starting point, go through to find the alignment
  //that starts from position at the 
  //error_positions[current_error_position_index]+1 but removing
  //right errors to see wehether we have a better score.
  for(unsigned i=_curr_error_position_index+1;i<_all_n_errors;
      i++)
    {
      //cout<<"\t-->loop:"<<i<<endl;
      //check for special case
      if(((signed int)_curr_error_position_index)==-1)
	{
	  //cout<<"\t\t-->special case"<<endl;
	  temp_n_errors=i;
	  temp_align_len=_all_error_positions[i]-1+1;
	  temp_score=CalculateScore(temp_align_len, 
				    _all_error_positions,
				    temp_n_errors, _cost);	
	}
      else
	{
	  temp_n_errors=i-_curr_error_position_index-1;
	  temp_align_len=_all_error_positions[i]-1
	    -_all_error_positions[_curr_error_position_index];
	  temp_score=CalculateScore(temp_align_len,
				    _all_error_positions+_curr_error_position_index+1,
				    temp_n_errors, _cost);
	}
      //cout<<"\t\t--->temp score:"<<temp_score<<endl;
      //now check to see whether this is the best,
      if(temp_score>best_score||
	 (temp_score==best_score&&temp_align_len>best_align_len))//update everything
	{
	  _best_n_errors=temp_n_errors;
	  best_align_len=temp_align_len;
	  best_score=temp_score;
	}
      //cout<<"\t\t--------end"<<endl;
    }
  //again, we did not explicitly return best_error_positions, it
  //it is redundent.
  return best_align_len;
}

//this is exactly identical to GetBestScoreAmongTheAlignment(),
//see above for usage.
unsigned  align_with_constraints_after_left_remove_right_errors
   (const unsigned& _total_align_len, 
    const unsigned* _all_error_positions,
    const unsigned& _all_n_errors,  
    const unsigned& _curr_error_position_index,
    const double& _cost,
    /*output*/ 
    unsigned& _best_n_errors  )
{
  //cout<<"@@@@@Calling surrogate function"<<endl;
  return GetBestScoreAmongTheAlignment(_total_align_len,
				       _all_error_positions,
				       _all_n_errors,
				       _curr_error_position_index,
				       _cost, _best_n_errors);
  
}

/*align the two strings, but need to find the best scores for
 *the subalignments that could remove both left and right errors
 *
 *input:
 *     seq1 and seq2: two string to align fixed left (both start
 *           at position 0
 *     _cost, error cost
 *output:
 *     _n_errors: number of errors after return
 *     _error_positions: holding the position, need 
 *               to be initialized by the caller and 
 *               the size of _v_allowed_errors
 *     _align_position_start
 */
unsigned align_with_constraints_fixed_left_remove_both_errors
(const string& _seq1, const string& _seq2, 
 const unsigned& _v_allowed_errors, 
 const double& _cost,
 /*output*/ unsigned& _n_errors, unsigned* _error_positions, 
 unsigned& _align_position_start)
{
  //cout<<"=======inside remove both errors function"<<endl;
  unsigned l_seq1=_seq1.size();
  unsigned l_seq2=_seq2.size();
  unsigned align_length=0;
  unsigned max_match_length=l_seq1;
  if(l_seq2<l_seq1)
    {
      max_match_length=l_seq2;
    }
  unsigned* temp_error_positions=new unsigned[max_match_length];
  _align_position_start=0;//starting from the beginning
//now we need to bitwise matching between two strings.
  //there are a few things to do for this
  //1) go through the string to find match and mismatch for each bit
  //    /*to make it right, we need to go one beyond the maximum allowed errors
  //    till we reach the one next error if possible*/
  //   for this function, we need to go through the whole length and
  //   will take care of max_error later
  
  //bool matchFlag;
  _n_errors=0;
  //  unsigned i;
 
  //cout<<"*****calling ...cost:"<<_cost<<endl;
  //do the alignment
  for(unsigned i=0;i<max_match_length;i++)
    {
      
      //check for match or mismatch
      if(_seq1.at(i)!=_seq2.at(i))
	{
	  //matchFlag=false;
	 temp_error_positions[_n_errors]=i;
	  _n_errors=_n_errors+1;
	}

    }
  align_length=max_match_length;
  /*cout<<"\tafter matching:"<<endl;
  cout<<"\tn_errors:"<<_n_errors<<endl;
  
  cout<<"\talength:"<<align_length<<endl;
  cout<<"\tn_errors:"<<_n_errors<<endl;
  cout<<"\talign_position_start:"<<_align_position_start<<endl;
  cout<<"\terror_positions:";
  
  for(unsigned int i=0;i<_n_errors;i++)
    {
      cout<<temp_error_positions[i]<<",";
    }
    cout<<endl;
  */
  //now we need to go through the alignment to see whether we end
  //up at error. if so we remove the error and stop at match
  //this might not be necessary, but will get a short sequence
  //so make the calculation faster.
  for(unsigned j=_n_errors-1;j<0;j--)
    {
      //go from back to front to check the error positions
      if(temp_error_positions[j]<align_length-1)
	{
	  break;
	}
      else //at this point, the error is at the very 
	//end of alignment, so we need to remove it and update 
	//align_length and 
	{
	  align_length=temp_error_positions[j];
	  _n_errors--;
	}
    }
  
  /*  cout<<"After removing the trailing errors"<<endl;
  cout<<"\talength:"<<align_length<<endl;
  cout<<"\tn_errors:"<<_n_errors<<endl;
  cout<<"\talign_position_start:"<<_align_position_start<<endl;
  cout<<"\terror_positions:";
  for(unsigned int i=0;i<_n_errors;i++)
    {
      cout<<temp_error_positions[i]<<",";
    }
    cout<<endl;
  */
  //now we are ready to do the alignment
  double best_score=CalculateScore(align_length,
				    temp_error_positions,
				    _n_errors, _cost);
  _align_position_start=0;
  unsigned best_align_length=align_length;
  unsigned best_n_errors=_n_errors;
  unsigned best_align_position_start=_align_position_start;
  unsigned* best_error_position=temp_error_positions;

  unsigned curr_error_position;
  unsigned temp_align_length;
  unsigned temp_n_errors;  
  double temp_best_score;
  //cout<<"Best score so far:"<<best_score<<endl;
  //cout<<"\tstart the looping...."<<endl;
  
  //move along the errors to remove left errors to see which gives
  //the best score
  for(signed int j=-1;j<(signed)_n_errors;j++)
    {
      //cout<<"%%%%%loop:"<<j<<";best error_position:"<<best_error_position<<endl;
      if(j==-1)
	{
	  curr_error_position=-1;
	}
      else
	{
	  curr_error_position=temp_error_positions[j];
	}
      if((align_length-1-curr_error_position)<best_score)
	{//in this case, the length is too small and we will not
	  //get a better score than best_score.
	  //cout<<"break out!!!"<<endl;
	  break;
	}
      //call alignment
      temp_align_length=align_with_constraints_after_left_remove_right_errors
	(align_length, temp_error_positions, _n_errors, (unsigned)j, _cost,
	 /*output*/ temp_n_errors  );
      if(j==-1)
	temp_best_score=CalculateScore(temp_align_length, 
				     temp_error_positions, 
				     temp_n_errors,
				       _cost);
      else
	temp_best_score=CalculateScore(temp_align_length, 
				     temp_error_positions+j+1, 
				     temp_n_errors,
				       _cost);
      /*cout<<"after aligning:"<<endl;
      cout<<"\talength:"<<temp_align_length<<endl;
      cout<<"\tn_errors:"<<temp_n_errors<<endl;
      cout<<"\tbest score:"<<temp_best_score<<endl;
      cout<<"\talign_position_start:"<<_align_position_start<<endl;
      cout<<"\terror_positions:";
      for(unsigned int i=j+1;i<_n_errors;i++)
	{
	  cout<<temp_error_positions[i]<<",";
	}
      cout<<endl;
      */

      //check for best score
      if(temp_best_score>best_score)
	{
	  //cout<<"&&&&&&&&---------find a better one"<<endl;
	  best_score=temp_best_score;
	  best_align_length=temp_align_length;
	  best_n_errors=temp_n_errors;
	  
	  if(j==-1)
	    {
	      best_align_position_start=0;
	      best_error_position=temp_error_positions;
	    }
	  else
	    {
	      best_align_position_start=temp_error_positions[j]+1;
	      best_error_position=temp_error_positions+j+1;
	    }
	}
    }//end of for loop
  /*
  cout<<"after aligning final:"<<endl;
  cout<<"\talength:"<<best_align_length<<endl;
  cout<<"\tn_errors:"<<best_n_errors<<endl;
  cout<<"\tbest score:"<<best_score<<endl;
  cout<<"\talign_position_start:"<<best_align_position_start<<endl;
  cout<<"\terror_positions:";
  for(unsigned int i=0;i<best_n_errors;i++)
    {
      cout<<best_error_position[i]<<",";
    }
  cout<<endl;
  */


  //cout<<"\t\t preparing the outputs........."<<endl;
  if(best_n_errors>_v_allowed_errors)
    {
      //not good
      //cout<<">>>>>>a failed alignment"<<endl;
      align_length=0;
      _n_errors=0;

    }
  else
    {
      align_length=best_align_length;
      _n_errors=best_n_errors;
      _align_position_start=best_align_position_start;
      
      for(unsigned m=0;m<_n_errors;m++)
	{
	  _error_positions[m]=best_error_position[m]-
	    _align_position_start;
	} 
      
    }
  //cout<<"======>cleaning up mem"<<endl;
  //clean up memory
  delete[] temp_error_positions;
  //cout<<"==========>done with align both errors"<<endl;
  return align_length;
}//nd of align remove both errors
 
