

calling commands
./ngsmapping_demux -s ./Isotype/top_1_clone1073_seqs.picked.fa.noIsotype.fa -f ./Isotype/top_1_clone1073_seqs.picked.fa.noIsotype.fa -b ./Isotype/top_1_clone1073_seqs.picked.fa.noIsotype.fa -a ./Isotype/top_1_clone1073_seqs.picked.fa.noIsotype.fa -x

	# to demux so to reading -s sequences.
	
	
../ngsmapping_demux -b Barcode7b_R1.fasta -a Barcode7b_R2.fasta -t -s Undetermined_S0_L002_R1_001_1000.fasta -q Undetermined_S0_L002_R2_001_1000.fasta -d >temp.txt
		#-t read indexes from the name of the sequences.
		#-s and -q are the sequence data 
		#-d is to do demux
		#-b and -a, the expected barcodes 
		
./ngsmapping_demux, the program used to demux or runs stats of indexes. It would read in
barcodes and sequence indexes and then match(instead of align) to find them.
We assume that the index (read) and barcode (expected) are having identical length.
Well, actually, the index could be long, but the barcodes are having fixed length.
we always match/align barcodes again index(as template) and we simply match (compare)
instead of aligning them. We also allow degeneracy on barcode, but not on indexes.
Therefore, the match score, mmf, is asymmetric. We also allow to read sequence 
index from the name of sequences instead of sequence read file. When sequences
are specified, we can do demux. Otherwise, we will only output stats about
barcoded sequences. We allow single or dual indexes. It is also possible to 
do 5' or 3' aligned on the second read depending on how the input specified.
	Option string: "hva:b:f:e:ptdrxs:n:" 
Usage:
	./ngsmapping_demux [-s sequence file] -a adaptor file (R1) [-b barcode file (R2)]
		 [-f  sequence index R1 file] [-e sequence index R2 file]
		 [-p (pair end read,true or false)] [-t (read indexes from the name of sequences] 
		 [-x (demux or not)] [-r (reverse implement r2 sequence index)]
		 [-q sequence file r2, optional]
		 [-m maximal num of mismatch allowed] [-d (dual index)]
		 [-h] [-v]

	-a: the input file barcode name (r2 for dual index input), 
		assuming has identical name and order
	-b: the input barcode file name (r1 or for single index input
	 a and b list the expected barcode to demux the input
	-f: the input file name read 1 index
	-e: the input file name read 2 index
	f and e contain the read index from sequencing to be demux'ed
	-p: paired end read, boolean, false by default
	-t: read indexs from the names of each sequence and assuming 
		the illumina style:
			  xxxxxx:xxxxx:xxxx:actgkd+acdfad
				separated by ':' and the last field contains the barcode (pair)
				false by default
	-x: demux (true or false), meaning to write out the sequences
 		or simple output the stats. False by default	-d: dual index (true or false), true by default
 		or simple output the stats	-r: reverse complement the read2 barcode or not to do demux
		only work for Read2, to reverse complement the sequence index 2,
		we never reverse complement the barcode.False(FivePrime) by default
	-q: sequence data file (R2), optional, could be
		specified or not. if dumux, then needed to check for this.
	-s: sequece data file, contains the real sequence data, 
		having the same order and name.
	-m: maximal number of mismatches allowed, 
		allow degenerate barcode, meaning the barcode (expected)
		 could be degenerated, but not the other way around.
		 barcode (expected)-->row, index (sequence read) ---->column
		 also this number is the total number of mismatches allow 
		 for both barcodes if there are dual indexeds. 1 by default
	-v and -h: version and help
Note: 1)when do comparison, always use barcode(expected) to 
	align index (data), in terms of alignment, barcode is sujbect,
	 index is pattern
	             index (pattern)
	             |||||
	             barcode (subject)
	2)allow degeneracy on barcode, but less so on index
	3)it is possible to have index longer than barcode, but not vice versa
	4)barcode and index must start on 5' position 1 (position 0 in c)
	5)we don't do alignment in comparison, but string comparison, 
		see comparison matrix, match.matrix.feng (mmf)
	6)the match score could be fraction number, because of degeneracy
	7)all the files should have the same name or order in order to reference them
			 *************@4/2/2014 by Feng
			*********updated by Feng @ BU 2018
			
####1/8/2019, run this to parse out the indexes, no demux, but reading from the output of
				#the ngsmapping_getbarcode, index1 and index2, allowing maximal 2 mismatches.
				#we want to see how the barcodes working on the data\
				#for WL02nano first try
###this below calling is correct, but the input indexes are not including all. are unique ones from ngsmapping_getbarcode run. so run the next one below this one.
./ngsmapping_demux -f Undetermined_S0_L001_R1_001.fastq.gzbar1.fasta -e Undetermined_S0_L001_R1_001.fastq.gzbar2.fasta -b WL02Nano_barcode1.fasta -a WL02Nano_barcode2.fasta -d -m2

####calling on reading the 