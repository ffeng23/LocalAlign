note: about this mapping isotype

It has not been tested extensively. so use with caution. But it should work. 

see examples in ./IsotypeMapping
see examples in ./IsotypeMapping_2 updated 11/7/2018

calling
./ngsmapping_Isotype -s top_1_clone1073_seqs.picked.fa.noIsotype.fa -f isotype_demux_fullLength.fasta -d 2 -g -5

things to pay attention in this case use 3' mapping (the program will rev
complement the isotype sequence before mapping) -d 2 #for 3'prime mapping. -d
1 for 5' mapping. For 5' or 3', we only mean to reverse complement the isotype
sequence, but not the sequence input (pattern for local alignment). the input
sequence will be kept unchanged.

also need to pay attention to minimum overlap length and offset as well as
match rate.

This has been written to use cutadapt stile.


using minimum overlap to control the alignment not too short. match rate will
control not too many errors. Since we are doing local alignment, the above
both will not conflict with each other. We could have very few errors (or no
errors), but very short alignment, which is not good and could be resulted in
by chance. offset will control the alignment/match will not be too far from
the end.
 
-d 2  for 3' mapping
-g 10 to gap opening of 10 might be good
-e 5 to extension of score of 5
-n 0.7 70% matching is a good empirical rate.

the rest should be using the default one.
