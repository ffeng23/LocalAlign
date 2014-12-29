#include "Alignment.hpp"

//see above for the definition (match_Vs)
//we should keep in mind that we are doing multithread functions
bool match_J(const SequenceString& _seq,
	      const GenomicV* _genJs, unsigned _numOfJSegs, 
	      unsigned _J_minimum_alignment_length, unsigned _J_maximum_deletion, 
	      unsigned _nagative_excess_deletion_max, unsigned _J_allowed_errors, unsigned _error_cost,
	      Alignment_Obj* _J
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
  vector<unsigned>* align_length = (J->align_length); //% contains the length of the alignment

  //% align_position(1,j) contains the first matched nt in seq (redundant, same as length(seq) - align_length +1 ),
  //% align_position_(2,j) is the first matched nt in genJ.seq{j} - also redundant information because we know where the sequencing primer ends.
  vector<unsigned[2]> align_position = J->align_position;
  vector<unsigned> n_errors=J->n_errors;//zeros(length(genJ),1); % number of mismatches in alignment

  vector<unsigned> >error_positions = zeros(length(genJ),J_allowed_errors); % positions of mismatches
											      min_deletions = zeros(length(genJ),1); % number of deletions implied by alignment

																	     j_large_deletion_flag=zeros(length(genJ),1); % Flag if deletions is too large
																										l_seq = length(seq);


  % Loop through alleles
      for j=1:length(genJ)
	      %j=1
	      %disp(['loop: ' num2str(j)])
	      % Get highest scoring alignment for this allele, with acceptable number of errors
							 % NOTE: I use an alignment function that fixes position on the left (as in match_Vs) but I pass in the
							 % reversed sequences, because I want to fix the position on the right.
							 l_target=length(genJ(j).seq);
  %[align_length(j), rev_align_position, n_errors(j), rev_error_positions] = align_with_constraints_fast_left_on_target( seq(end:-1:1) , genJ(j).seq(end:-1:1) , J_allowed_errors, error_cost);
  [align_length(j), rev_align_position, n_errors(j), rev_error_positions]=align_with_constraints_fast_no_fix(seq(end:-1:1) , genJ(j).seq(end:-1:1) , J_allowed_errors, error_cost);
      % Reverse the error positions to be relative to left-right orientation.
	  %%%j  error positions is relative to the end of J chain
	  if n_errors(j)>0
		       error_positions(j,1:n_errors(j)) = (rev_error_positions(n_errors(j):-1:1)) ;
      end

	% Calculate (redundant) align_position values. These are used in model scripts.
	align_position(:,j) = [ (l_seq-(rev_align_position(1)+align_length(j)-1)+1) , l_target-(rev_align_position(2)+align_length(j)-1)+1 ];

          % Calculate deletions implied by alignment
	      min_deletions(j) = align_position(2,j) - 1;
          % Flag if number of deletions is too many
			     if( min_deletions(j) > J_maximum_deletion)
			       j_large_deletion_flag(j) = 1;
      %continue;
          end
	    end

	    % Check!
	    assert(all(min_deletions>=0));

	  % sort the genomic Js in descending order of 'score' and ascending order of min_deletions
	      scores  = align_length - error_cost*n_errors;
	  S = [-scores, min_deletions];
	  [~,order]=sortrows(S);

	  % Set a score threshold for alleles.
	      min_score = max(scores) - 3*J_minimum_alignment_length;

	  % Subset of alleles that have high enough score, long enough alignment, and
	      % not too many deletions.
	      j_list = find( scores(order) >= min_score & align_length(order) >= J_minimum_alignment_length & j_large_deletion_flag(order)==0);

	  if isempty(j_list)
		      %disp(['No reasonable J match for sequence' ])
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
		[~,Jg] = unique(genJ_ok_gene_inds,'first' ); % Pick the first allele listed for each gene.
												       Jg = sort(Jg); % sort it because unique returns in ascending order of ok_order. This puts it back in order of best match.
															  % Only one (best) allele from each gene
															  J.alleles_from_distinct_genes = ok_order(Jg');

    %Now go through and find palindromic nucleotides for various number of deletions,
    %as well as error positions for 'negative' deletions.
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

        % Now for the 'negative' deletions, store where the mismatches are.
        % This is so we can count number of errors for these possibilities.

        % Number of 'negative' deletions to consider:
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
  


  return true;
}




