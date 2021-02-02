

calling commands
./ngsmapping_demux -s ./Isotype/top_1_clone1073_seqs.picked.fa.noIsotype.fa -f ./Isotype/top_1_clone1073_seqs.picked.fa.noIsotype.fa -b ./Isotype/top_1_clone1073_seqs.picked.fa.noIsotype.fa -a ./Isotype/top_1_clone1073_seqs.picked.fa.noIsotype.fa -x

	# to demux (-x) so to reading -s sequences.
	
../ngsmapping_demux -b Barcode7b_R1.fasta -a Barcode7b_R2.fasta -t -s Undetermined_S0_L002_R1_001_1000.fasta -q Undetermined_S0_L002_R2_001_1000.fasta -d >temp.txt
		#-t read indexes from the name of the sequences.
		#-s and -q are the sequence data 
		#-d is to do dual index
		#-b and -a, the expected barcodes 
######====doing on fastq files and
    #read from sequence names the index and do pairend and demux and dual index.
../ngsmapping_demux -b Barcode7b_R1.fasta -a Barcode7b_R2.fasta -t -s Undetermined_S0_L001_R1_001_200.fastq -q Undetermined_S0_L001_R2_001_200.fastq -x -d -p 

###08/22/2019 run test on read data WL03_TS03 undetermined first Undetermined_S0_L001_R1_001_200
../ngsmapping_demux -b Barcode7b_R1.fasta -a Barcode7b_R2.fasta -t -s Undetermined_S0_L001_R1_001_200.fastq -q Undetermined_S0_L001_R2_001_200.fastq -x -d -p 
		
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
There are a few things to mention:
         1) pair end read vs dual indexes. they are two different things.
         it is possible to have pair end read, but no dual indexes (single index in this case).
         it is not possilbe to have dual indexes on single end read data (program will quit on this).
         2) index vs barcode. index is the sequence(s) on the read data to indicate the "identity" 
                of the sequences. Barcode is the expected index sequences as reference.
         3) degenaracy vs mismatch. we allow degeneracy only on barcode, but not on index.
         the mismatch is unmatched nt at the same posible between index and barcode.
        Option string: "hva:b:f:e:ptdrxs:n:" 
Usage:
        ./ngsmapping_demux [-s sequence file (R1 if this is pair end read)] [-q sequence file r2, optional]
                 -b barcode file (R1) [-a adaptor file (R2) ]
                 [-f  sequence index R1 file] [-e sequence index R2 file]
                 [-p (pair end read,true or false)]
                 [-t (read indexes from the name of sequences] 
                 [-x (demux or not)] 
                 [-r (reverse implement r2 sequence index)]
                 [-m maximal num of mismatch allowed]
                 [-d (dual index)]
                 [-l length of barcode/index]
                 [-h] [-v]

        -b: the input barcode file name (r1 or for single index input
        -a: the input file barcode name (r2 for dual index input), 
                assuming has identical name and order
         **a and b list the expected barcode to demux/stat the inputs
        -s: sequece data file, contains the real sequence data, 
                having the same order and name.
        -q: sequence data file (R2), optional, could be
                specified or not. if demux and pair end, then needed to check for this.
        -f: the input file name for read 1 index
        -e: the input file name for read 2 index
        **f and e contain the read index from sequencing
        -p: paired end read, boolean, FALSE by default
        -t: read indexs from the names of each sequence and assuming 
                the illumina style:
                          xxxxxx:xxxxx:xxxx:actgkd+acdfad
                                separated by ':' and the last field contains the barcode (pair).
                even in the case of dual index and pair end read, we only read the seq 1 to get the index pair.
                assuming seq r1 contains both indexes by xxx+xxx
                                FALSE by default
        -x: demux (true or false), meaning to write out the sequences by barcodes
                or simple output the stats about barcode groups (when it is false). False by default
        -d: dual index (true or false), TRUE by default
        -r: reverse complement the read2 indexes/sequences for matching.***do we need to rc th sequence output???need to add elaboration.
                only work for Read2, to reverse complement the sequence index 2,
                we never reverse complement the barcode.False(FivePrime) by default
        -m: maximal number of mismatches allowed, 
                allow degenerate barcode, meaning the barcode (expected)
                 could be degenerated, but not the other way around.
                 barcode (expected)-->row, index (sequence read) ---->column
                 also this number is the total number of mismatches allow 
                 for both barcodes if there are dual indexeds. 1 by default
        -l: length of barcodes. 8 by defult
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
example scenarios:
         1)reading from index file, single end reads and single index 
         2)reading from index file, pair end reads and single index 
         3)reading from index file, pair end reads and dual index 
         4)reading from index file, paire end reads and dual index with reverse complement on read 2 index
         5)-8) case 1)-4) with read index from sequence name 
         9)-16) case 1)-8) with demux 
         17-32) case 1)-16) with mismatch 0 or 2+ 
         33)-64) case 1)-32) with different length 10 and case where barcode length specified (-l)
                 is short than file input or shorter than index or vice versa. 
         65)-128) case 1)-64) with fasta fastq gzip
                         *************@4/2/2014 by Feng
                        *********updated by Feng @ BU 2018
                        *********updated by Feng @ BU Aug 2019

			
####1/8/2019, run this to parse out the indexes, no demux, but reading from the output of
				#the ngsmapping_getbarcode, index1 and index2, allowing maximal 2 mismatches.
				#we want to see how the barcodes working on the data\
				#for WL02nano first try
###this below calling is correct, but the input indexes are not including all. are unique ones from ngsmapping_getbarcode run. so run the next one below this one.
./ngsmapping_demux -f Undetermined_S0_L001_R1_001.fastq.gzbar1.fasta -e Undetermined_S0_L001_R1_001.fastq.gzbar2.fasta -b WL02Nano_barcode1.fasta -a WL02Nano_barcode2.fasta -d -m2

####calling on reading the 
