note: about this mapping isotype

It has not been tested extensively. so use with caution. But it should work. 

see examples in ./IsotyeMapping

calling
./ngsmapping_Isotype -s top_1_clone1073_seqs.picked.fa.noIsotype.fa -f isotype_demux_fullLength.fasta -d 2 -g -5

things to pay attention in this case use 3' mapping (the program will rev
complement the isotype sequence before mapping)

also need to pay attention to minimum overlap length and offset as well as
match rate.
