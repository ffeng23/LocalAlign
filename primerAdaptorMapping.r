###############this is R code to read fasta file and then do some preliminary diagnosis######################
########the job to do is to align 1)adaptor if they are not trimmed; 2)primer; 3)barcode.
##########the idea is to use "overlap" (means no end gap penalty to align), finnd the best score ones.
##########also, here we generate some statistics
##########last thing, we can do trimming too.

#source("http://bioconductor.org/biocLite.R")
# biocLite("Biostrings")

 #load
 library('Biostrings')

##########section to prepare the adaptor, primer, barcode##############
########this section might need to change according different platform, different primer sets, barcoded, etc##########

########for 454 Roche, with 9 barcodes, adaptor A,B and NestedUPM and IgG/M/D
MID1<-DNAString("ACGAGTGCGT"); MID2<-DNAString("ACGCTCGACA"); MID3<-DNAString("AGACGCACTC"); MID4<-DNAString("AGCACTGTAG");
MID5<-DNAString("ATCAGACACG"); MID6<-DNAString("ATATCGCGAG");
MID7<-DNAString("CGTGTCTCTA"); MID8<- DNAString("CTCGCGTGTC");
MID9<- DNAString("TAGTATCAGC");
Barcodes<-list(MID1, MID2, MID3, MID4, MID5, MID6, MID7, MID8, MID9);
Barcode_names<-list("MID1", "MID2", "MID3", "MID4", "MID5", "MID6", "MID7", "MID8", "MID9");

adaptorA<-DNAString("CGTATCGCCTCCCTCGCGCCATCAG");
adaptorB<-DNAString("CTATGCGCCTTGCCAGCCCGCTCAG");

##for the universal primer, it is kind of tricky, the long UPM+Magic code is the 
##whole string introduced on one end, but the nest UPM is the one within the
## LongUPM+MagicCode, so the actually introduced on the final reads is the
## nestedUPM+restUPM+Magic
NestedUPM<-DNAString("ACGACTCACTATAGGGCAAGCAG");
RestOfUPMWithOligo<-DNAString("TGGTAACAACGCAGAGTACGCGGG");

FullUPMSeqenceOnSequence<-c(NestedUPM, RestOfUPMWithOligo);


###constant side of primer
###only the Nested IgG/M/D
NestedIgG<-DNAString("gttcggggaagtagtccttgac");
NestedIgM<-DNAString("caggagacgagggggaa");
NestedIgD<-DNAString("acgttgggtggtacccagttat");

#####################################
####now we need to put together them for full combination in order to run them
### alignment against sequence reads for mapping
Adaptor<-list(adaptorA, adaptorB);
#Barcode<-list(MID1,MID2,MID3,MID4,MID5,MID6,MID7,MID8,MID9);
Constant<-list(NestedIgG, NestedIgM, NestedIgD);
Constant_names<-list("NestedIgG", "NestedIgM", "NestedIgD");
UPM<-list(FullUPMSeqenceOnSequence);
UPM_names<-list("Nested UPM + TS Oligo");

########create set holding both ends sequences
ForwardSet<-list();
ForwardNameSet<-list();
ReverseSet<-list();
ReverseNameSet<-list();


########structor always is adaptor+Barcode+primer
###adaptorA is mapping the foward 5'-3'side/beginning 

i<-1;j<-1;k<-1;
for( j in 1:length(Barcodes))
{
	for(k in 1:length(Constant))
	{
		temp<-c(adaptorA, Barcodes[[j]], Constant[[k]]);
		
		ForwardSet<-c(ForwardSet,temp);
		temp_name<-paste("Adaptor A +", Barcode_names[[j]],"+", Constant_names[[k]]);
		ForwardNameSet<-c(ForwardNameSet,temp_name );
	}
	#the below part can happen, altough, in our current setting they are not used, we keep it here for future
	for(m in 1:length(UPM))
	{
		temp<-c(adaptorA,Barcodes[[j]], UPM[[m]]);
		ForwardSet<-c(ForwardSet,temp);
		ForwardNameSet<-c(ForwardNameSet, paste("Adaptor A +", Barcode_names[[j]], "+", UPM_names[[m]]));
	}
}
###adaptorB is mapping the reverse 3'-5'side/ending
for( j in 1:length(Barcodes))
{
	
	for(m in 1:length(UPM))
	{
		temp<-c(adaptorB,Barcodes[[j]], UPM[[m]]);
		ReverseSet<-c(ReverseSet,temp);
		ReverseNameSet<-c(ReverseNameSet, paste("Adaptor B +", Barcode_names[[j]], "+", UPM_names[[m]]));
	}
	###this follow is not possible in the current setting, but we keep it for the future.
	for(k in 1:length(Constant))
	{
		temp<-c(adaptorB, Barcodes[[j]], Constant[[k]]);
		ReverseSet<-c(ReverseSet,temp);
		ReverseNameSet<-c(ReverseNameSet, paste("Adaptor B +", Barcode_names[[j]], "+", Constant_names[[k]]));
	}
}

####how far we allow the alignment to be away from the ends. can not be too far, since they are supposed to be aligned on the ends.

offsetForward<-10;###10 might too big????
offsetReverse<-10;

#########now we need to read in the data
####first the score matrix, we used NUC.4.4, NUC.4.4 is the file on the disk, 
######we need to find the good working directory
##########code example############
## Compare our BLOSUM62 with BLOSUM62 from ftp://ftp.ncbi.nih.gov/blast/matrices/
##  data(BLOSUM62)
##  BLOSUM62["Q", "Z"]
##  file <- "ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62"
##  b62 <- as.matrix(read.table(file, check.names=FALSE))
##  b62["Q", "Z"]

#setwd("E:/MSI_software/deepsequencing")
NUC44 <- as.matrix(read.table("NUC.4.4", check.names=FALSE))
#define gap
gapOpening<- -4;
gapExtension<- -4;

#####now read the input file with reads
ff <- readDNAStringSet("IgG_Region1_Reads2000.fasta", format="fasta")#, strip.descs=TRUE );#get rid of ">" symbol for desc
#desc <- sapply(ff, function(x) x$desc);  #get description/name line

####now we start doing the job
i<-46;j<-1;k<-1;
MismatchRateThreshold<-0.80; ###not too many mismatch.
MinimumOverlapLength<-10; ##with barcode and primer, if with adaptor this should be even longer, otherwise 10minimum

####the list for holding the output
mappedReadsBoth<-NULL;
mappedReadsForward<-NULL;
mappedReadsReverse<-NULL;
mappedReadsNone<-NULL;
numOfSeq<-10; ##this is maximum number of sequences each file can hold

#for statistics
mappedBothLen_arr<-c();
mappedForwardLen_arr<-c();
mappedReverseLen_arr<-c();
mappedNoneLen_arr<-c();

##
#fileCounter_mpBoth<-1;
#fileCounter_mpForward<-1;
#fileCounter_mpReverse<-1;
#fileCounter_mpNone<-1;


for(i in 1:length(ff))
{
	tempLen<-nchar(ff[[i]]);
	if(i%%50==0)
	{
		cat("..",i,"/",length(ff));
	}
	#align each forward and reverse sequence
	bestForwardScore<- -10000000;
	bestForwardAlign<-NULL;
	bestReverseScore<- -1000000;
	bestReverseAlign<-NULL;
	
	bestForwardIndex<-0;
	bestReverseIndex<-0;
	
	foundForwardFlag<-FALSE;
	foundReverseFlag<-FALSE;
	
	###forward has to be mapped to the beginning!!!
	for(j in 1:length(ForwardSet))
	{
		overlapAlign <-    pairwiseAlignment(ff[[i]], ForwardSet[[j]], type = "overlap", substitutionMatrix = NUC44, gapOpening = gapOpening, gapExtension = gapExtension)
		
		###need to get the mismatch rate
		strPattern<-as.character(pattern(overlapAlign));
		strSubject<-as.character(subject(overlapAlign));
		strPattern<-strsplit(strPattern,"");
		strSubject<-strsplit(strSubject,"");
		strPattern<-strPattern[[1]];
		strSubject<-strSubject[[1]];
		lenOfOverlap<-length(strPattern);##total length of overlap, including the mismatch
		mismatch_rate<-length(strPattern[strPattern==strSubject])/lenOfOverlap;
		#check the best score and mismatch rate
		if(score(overlapAlign)>bestForwardScore&&mismatch_rate>MismatchRateThreshold&&lenOfOverlap>MinimumOverlapLength)
		{
			##we need to see the offset too, only important to the pattern, here the 
			##the pattern is the long one, we align the primer sequence against
			###the primer should be in the beginning, not too far,
			if(start(pattern(overlapAlign))<offsetForward )
			{##we good
				bestForwardScore<-score(overlapAlign);
				bestForwardAlign<-overlapAlign;#with name
				bestForwardIndex=j;
				foundForwardFlag<-TRUE;
				
				##need to build a data structure to 
			}
		}
	}
	
	####reverse side should be mapped on the end of the reads, need to reverse complement the sequence too
	for(k in 1:length(ReverseSet))
	{
		reverseComplementReverse<-reverseComplement(ReverseSet[[k]])
		overlapAlign <-    pairwiseAlignment(ff[[i]], reverseComplementReverse, type = "overlap", substitutionMatrix = NUC44, gapOpening = gapOpening, gapExtension = gapExtension)
		###need to get the mismatch rate
		strPattern<-as.character(pattern(overlapAlign));
		strSubject<-as.character(subject(overlapAlign));
		strPattern<-strsplit(strPattern,"");
		strSubject<-strsplit(strSubject,"");
		strPattern<-strPattern[[1]];
		strSubject<-strSubject[[1]];
		lenOfOverlap<-length(strPattern);##total length of overlap, including the mismatch
		mismatch_rate<-length(strPattern[strPattern==strSubject])/lenOfOverlap;
		
		#check the best score
		if(score(overlapAlign)>bestReverseScore&&mismatch_rate>MismatchRateThreshold&&lenOfOverlap>MinimumOverlapLength)
		{
			##we need to see the offset too, only important to the pattern, here the 
			##the pattern is the long one, we align the primer sequence against
			###the primer should on the both ends, not too far,
			if( 
				(nchar(ff[[i]])-end(pattern(overlapAlign)))<offsetReverse )
			{##we good
				bestReverseScore<-score(overlapAlign);
				bestReverseAlign<-overlapAlign;#with name
				bestReverseIndex<-k;
				foundReverseFlag<-TRUE;
			}
		}
	}
	
	###done with mapping, now we need to make the output ready
		
	#the forward 
	leadingSpaceOriginal<-NULL;
	leadingSpaceForward<-NULL;
	leadingSpaceReverse<-NULL;
	if(foundForwardFlag)
	{
		###we need to figure out the leading spaces in front of the original sequences
		startOriginal=start(pattern(bestForwardAlign));
		startSubject=start(subject(bestForwardAlign));
		if(startSubject>startOriginal)
		{
			#leadingSpaceOriginal=as.character();
			leadingSpaceOriginal<-paste(rep("-",startSubject-startOriginal), collapse = '')
		}else
		{
			#leadingSpaceForward=as.character();
			leadingSpaceForward<-paste(rep("-",startOriginal-startSubject), collapse = '')
		}
		
		#now we need to add the aligned sequence to replace the original one
		#we assume the original one is longer than the aligned one, it has to be
		 #this is the intial part
		replaceOne<-substr(ForwardSet[[bestForwardIndex]],1, start(subject(bestForwardAlign))-1);
		 #aligned part
		replaceOne<-paste(replaceOne, subject(bestForwardAlign),sep="");
		 #last part unaligned
		replaceOne<-paste(replaceOne, substr(ForwardSet[[bestForwardIndex]],end(subject(bestForwardAlign)), nchar(ForwardSet[[bestForwardIndex]])), sep="");
		tempLstF<-list(desc=ForwardNameSet[[bestForwardIndex]], seq=replaceOne)
	}else
	{
		#no need to add leading space
		tempLstF<-list(desc="NoMatch", seq="");
	}
	############here to do!!!!!!!!!
	
	tempLstF$seq<-paste(leadingSpaceForward, tempLstF$seq, sep="");
	
	#the reverse
	if(foundReverseFlag)
	{
		#now we need to add the aligned sequence to replace the original one
		#we assume the original one is longer than the aligned one, it has to be
		 #this is the intial part
		rcReverseSeq=reverseComplement(ReverseSet[[bestReverseIndex]])
		replaceOneR<-substr(rcReverseSeq,1, start(subject(bestReverseAlign))-1);
		 #aligned part
		replaceOneR<-paste(replaceOneR, subject(bestReverseAlign),sep="");
		 #last part unaligned
		replaceOneR<-paste(replaceOneR, substr(rcReverseSeq,end(subject(bestReverseAlign))+1, nchar(rcReverseSeq)), sep="");
		tempLstR<-list(desc=ReverseNameSet[[bestReverseIndex]], seq=replaceOneR)
		
		#now we need to figure out how the leading space to put in front of reverse one
		if(start(pattern(bestReverseAlign))>=start(subject(bestReverseAlign)))
		{
			leadingSpaceReverse<-paste( rep("-",start(pattern(bestReverseAlign))-start(subject(bestReverseAlign))), collapse='');
			tempLstR$seq<-paste(leadingSpaceReverse, tempLstR$seq, sep="");
		}else
		{#here, in this case, the adaptor+primer is longer than the seqs just by alignment, then we need to simply remove some leading part of the adaptor primer
			tempLstR$seq<-substr(tempLstR$seq, start(subject(bestReverseAlign))-start(pattern(bestReverseAlign))+1, nchar(as.character(tempLstR$seq)));
		}
	}else
	{
		tempLstR<-list(desc="NoMatch", seq="");
	}
	
	#now we need to take care of the read sequence alignment string
	#on the reverse part first
	spaceCarryOverFromForwardToReverse=0;####this is the leading space for reverse adaptor primer, because the insertion in the forward alignment
	if(foundReverseFlag)
	{
		replaceOneSeq<-substr(ff[[i]],1, start(pattern(bestReverseAlign))-1);
		 #aligned part
		replaceOneSeq<-paste(replaceOneSeq, pattern(bestReverseAlign),sep="");
		 #last part unaligned
		replaceOneSeq<-paste(replaceOneSeq, substr(ff[[i]],end(pattern(bestReverseAlign))+1, nchar(ff[[i]])), sep="");
		ff[[i]]<-DNAString(replaceOneSeq);
	}
	if(foundForwardFlag)
	{
		replaceOneSeq<-substr(ff[[i]],1, start(pattern(bestForwardAlign))-1);
		 #aligned part
		replaceOneSeq<-paste(replaceOneSeq, pattern(bestForwardAlign),sep="");
		 #last part unaligned
		replaceOneSeq<-paste(replaceOneSeq, substr(ff[[i]],end(pattern(bestForwardAlign))+1, nchar(ff[[i]])), sep="");
		ff[[i]]<-DNAString(replaceOneSeq);
		spaceCarryOverFromForwardToReverse=nchar(pattern(bestForwardAlign))- (end(pattern(bestForwardAlign))-start(pattern(bestForwardAlign))+1);	
	}
	tempStr<-paste(rep("-",spaceCarryOverFromForwardToReverse), collapse='');
	tempStr<-paste(tempStr, as.character(tempLstR$seq), sep="");
	tempLstR$seq<-paste(leadingSpaceOriginal, tempStr, sep="");
	ff[[i]]<-DNAString(paste(leadingSpaceOriginal, ff[[i]], sep=""))
	
	#now put the sequences to the correct files
	if(foundForwardFlag)
	{
		if(foundReverseFlag)
		{
			tempStrSet<-DNAStringSet(ff[[i]]);
			names(tempStrSet)[1]=names(ff)[i];
			mappedReadsBoth<-append(mappedReadsBoth, tempStrSet);
			tempStrSet<-DNAStringSet(tempLstF$seq);
			names(tempStrSet)[1]=tempLstF$desc;
			mappedReadsBoth<-append(mappedReadsBoth, tempStrSet);
			tempStrSet<-DNAStringSet( tempLstR$seq);
			names(tempStrSet)[1]=tempLstR$desc;
			mappedReadsBoth<-append(mappedReadsBoth, tempStrSet);
			
			mappedBothLen_arr<-c(mappedBothLen_arr,tempLen);
		}else
		{
			tempStrSet<-DNAStringSet(ff[[i]]);
			names(tempStrSet)[1]=names(ff)[i];
			mappedReadsForward<-append(mappedReadsForward, tempStrSet);
			
			tempStrSet<-DNAStringSet(tempLstF$seq);
			names(tempStrSet)[1]=tempLstF$desc;
			mappedReadsForward<-append(mappedReadsForward, tempStrSet);
			
			tempStrSet<-DNAStringSet( tempLstR$seq);
			names(tempStrSet)[1]=tempLstR$desc;
			mappedReadsForward<-append(mappedReadsForward, tempStrSet);
			
			#mappedReadsForward<-c(mappedReadsForward, list(ff[[i]], tempLstF, tempLstR));
			mappedForwardLen_arr<-c(mappedForwardLen_arr,tempLen);
		}
	}else
	{
		if(foundReverseFlag)
		{
			tempStrSet<-DNAStringSet(ff[[i]]);
			names(tempStrSet)[1]=names(ff)[i];
			mappedReadsReverse<-append(mappedReadsReverse, tempStrSet);
			
			tempStrSet<-DNAStringSet(tempLstF$seq);
			names(tempStrSet)[1]=tempLstF$desc;
			mappedReadsReverse<-append(mappedReadsReverse, tempStrSet);
			
			tempStrSet<-DNAStringSet( tempLstR$seq);
			names(tempStrSet)[1]=tempLstR$desc;
			mappedReadsReverse<-append(mappedReadsReverse, tempStrSet);
		
			#mappedReadsReverse<-c(mappedReadsReverse, list(ff[[i]], tempLstF, tempLstR));
			mappedReverseLen_arr<-c(mappedReverseLen_arr,tempLen);
		}else
		{
			tempStrSet<-DNAStringSet(ff[[i]]);
			names(tempStrSet)[1]=names(ff)[i];
			mappedReadsNone<-append(mappedReadsNone, tempStrSet);
			
			tempStrSet<-DNAStringSet(tempLstF$seq);
			names(tempStrSet)[1]=tempLstF$desc;
			mappedReadsNone<-append(mappedReadsNone, tempStrSet);
			
			tempStrSet<-DNAStringSet( tempLstR$seq);
			names(tempStrSet)[1]=tempLstR$desc;
			mappedReadsNone<-append(mappedReadsNone, tempStrSet);
			
			#mappedReadsNone<-c(mappedReadsNone, list(ff[[i]], tempLstF, tempLstR));
			mappedNoneLen_arr<-c(mappedNoneLen_arr,tempLen);
		}
	}
	
	if(i%%numOfSeq==0) #write once every 1000 sequences
	{
		if(! is.null(mappedReadsBoth))
		{
			writeXStringSet(mappedReadsBoth, file="mappedBoth.fas", format="fasta", append=TRUE,width=100);
			#fileCounter_mpBoth<-fileCounter_mpBoth+1;
			mappedReadsBoth<-NULL;
		}
		if(! is.null(mappedReadsForward))
		{
			writeXStringSet(mappedReadsForward, file="mappedForward.fas", format="fasta",append=TRUE,width=100);
			#fileCounter_mpForward<-fileCounter_mpForward+1;
			mappedReadsForward<-NULL;
		}
		if(! is.null(mappedReadsReverse))
		{
			writeXStringSet(mappedReadsReverse, file="mappedReverse.fas", format="fasta",append=TRUE,width=100);
			#fileCounter_mpReverse<-fileCounter_mpReverse+1;
			mappedReadsReverse<-NULL;
		}
		if(! is.null(mappedReadsNone))
		{
			writeXStringSet(mappedReadsNone, file="mappedNone.fas",format="fasta",append=TRUE,width=100);
			#fileCounter_mpNone<-fileCounter_mpNone+1;
			mappedReadsNone<-NULL;
		}

		#also write out the stats
		write(mappedBothLen_arr, "mappedBoth_stats.txt", append=TRUE);
		mappedBothLen_arr<-c();
		write(mappedForwardLen_arr, "mappedForward_stats.txt", append=TRUE);
		mappedForwardLen_arr<-c();
		write(mappedReverseLen_arr, "mappedReverse_stats.txt",append=TRUE);
		mappedReveseLen_arr<-c();
		write(mappedNoneLen_arr, "mappedNone_stats.txt",append=TRUE);
		mappedNoneLen_arr<-c();

	}#end of each 1000 seqs read write
}
#now write
	if(! is.null(mappedReadsBoth))
	{
		writeXStringSet(mappedReadsBoth, file="mappedBoth.fas", format="fasta",append=TRUE, width=100);
		#fileCounter_mpBoth<-fileCounter_mpBoth+1;
		#mappedReadsBoth<-NULL;
	}
	if(! is.null(mappedReadsForward))
	{
		writeXStringSet(mappedReadsForward, file="mappedForward.fas", format="fasta",append=TRUE,width=100);
		#fileCounter_mpForward<-fileCounter_mpForward+1;
		#mappedReadsForward<-NULL;
	}
	if(!is.null(mappedReadsReverse))
	{
		writeXStringSet(mappedReadsReverse, file="mappedReverse.fas", format="fasta",append=TRUE,width=100);
		#fileCounter_mpReverse<-fileCounter_mpReverse+1;
		#mappedReadsReverse<-NULL;
	}
	if(!is.null(mappedReadsNone))
	{
		writeXStringSet(mappedReadsNone, file="mappedNone.fas",format="fasta",append=TRUE,width=100);
		#fileCounter_mpNone<-fileCounter_mpNone+1;
		#mappedReadsNone<-NULL;
	}
	
	write(mappedBothLen_arr, "mappedBoth_stats.txt",append=TRUE);
	write(mappedForwardLen_arr, "mappedForward_stats.txt",append=TRUE);
	write(mappedReverseLen_arr, "mappedReverse_stats.txt",append=TRUE);
	write(mappedNoneLen_arr, "mappedNone_stats.txt",append=TRUE);
	
	
	
