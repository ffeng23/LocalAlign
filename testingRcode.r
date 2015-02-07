
 #load
 library('Biostrings')


NUC44 <- as.matrix(read.table("NUC.4.4", check.names=FALSE))
BLOSUM50<-as.matrix(read.table("BLOSUM50", check.names=FALSE))
tm<-as.matrix(read.table("TestMatrix", check.names=FALSE))
tm2<-as.matrix(read.table("TestMatrix2", check.names=FALSE))
#define gap
gapOpening<- -8;
gapExtension<- -8;

seq1<-AAString("VSPAGMASGYDPGKA");
seq2 <-AAString("IPGKATREYDVSPAG");

LocalAlign <-    pairwiseAlignment(seq1, seq2, type = "global", substitutionMatrix = BLOSUM50, gapOpening = gapOpening, gapExtension = gapExtension)


seq1<-AAString("ATAT");
seq2 <-AAString("ACHKAT");

seq1<-AAString("ATCGA");
seq2 <-AAString("GATTGA");
seq1<-AAString("AGCTAGAGACCAGTCTGAGGTAGA");
seq2 <-AAString("AGCTAGAGACCAGCTATCTAGAGGTAGA");

seq1<-DNAString("CCAATCTACTACTGCTTGCAGTAC");

seq2<-DNAString("AGTCCGAGGGCTACTCTACTGAAC");

seq1<-DNAString("ATGACGGTGATGATGGTAGTGGCCATTTCTTTGCCTACCCACTGTGGCCTTGTTGTGAAGAACTTTGATGCCCTGTTGTATGTGATCGTTATGACAATCCTGTGAAGTTTTTCAATGATGAAGGAAACTGAGGCTTAGAGAGTTGGAGTAATCTGTGAAGGCTCAGGATGGGCAAGAGGTCCCGCCCAGGTTTGAGCCCAGATGCGAGGTTACCACGCTTCCTGGTGAGGTGTTTTACAACTAAGGCCAAGCCAGGCAAAACCCATTGTTCTGCAGCTTCTGGCTTGGATTGGGTGTCTTGTTGAGTATGTGGGCAGTGGATCTGATGTTTTCCACTTCCACCAAGGGCCCATCGGTCTTCCCCCTGT");

seq2<-DNAString("GTATTATGATTACGTTTGGGGGAGTTATCGTTATACC");



LocalAlign <-    pairwiseAlignment(seq1, seq2, type = "local", substitutionMatrix = NUC44, gapOpening = gapOpening, gapExtension = gapExtension)
