###Makefile for LocalAlign
###By Feng 03/08/2014 @bu
#######v1.0

###Define variables######

GCC = gcc
GXX = g++

CFLAG= -Wall -g -Werror -O -std=c++11
LOADFLAG=-s -lm -lz

CXXFLAG=${CFLAG}

ACCESSDIR=Accessory/

SRCS_0=main.cpp ${ACCESSDIR}string_ext.cpp score.cpp ${ACCESSDIR}SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp LocalAlignment.cpp GlobalAlignment.cpp OverlapAlignment.cpp ${ACCESSDIR}FastaHandler.cpp SequenceHandler.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp TracebackTable.cpp SequenceHandlerCommon.cpp ${ACCESSDIR}GzTools.cpp ${ACCESSDIR}FileHandler.cpp ${ACCESSDIR}FASTQ.cpp ${ACCESSDIR}FastqHandler.cpp

SRCS_1=	NGSMapping_Adaptor_main.cpp ${ACCESSDIR}string_ext.cpp score.cpp ${ACCESSDIR}SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp ${ACCESSDIR}FastaHandler.cpp SequenceHandler.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp TracebackTable.cpp SequenceHandlerCommon.cpp LocalAlignment.cpp ${ACCESSDIR}GzTools.cpp ${ACCESSDIR}FileHandler.cpp ${ACCESSDIR}FASTQ.cpp ${ACCESSDIR}FastqHandler.cpp

SRCS_2=	NGSMapping_PrimerDimer_main.cpp ${ACCESSDIR}string_ext.cpp score.cpp ${ACCESSDIR}SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp ${ACCESSDIR}FastaHandler.cpp SequenceHandler.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp TracebackTable.cpp SequenceHandlerCommon.cpp LocalAlignment.cpp ${ACCESSDIR}GzTools.cpp ${ACCESSDIR}FileHandler.cpp ${ACCESSDIR}FASTQ.cpp ${ACCESSDIR}FastqHandler.cpp

SRCS_3=	NGSMapping_Constant_main.cpp ${ACCESSDIR}string_ext.cpp score.cpp ${ACCESSDIR}SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp ${ACCESSDIR}FastaHandler.cpp SequenceHandlerConstant.cpp TracebackTable.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp LocalAlignment.cpp SequenceHandlerCommon.cpp ${ACCESSDIR}GzTools.cpp ${ACCESSDIR}FileHandler.cpp ${ACCESSDIR}FASTQ.cpp ${ACCESSDIR}FastqHandler.cpp

SRCS_4=	NGSMapping_Isotype_main.cpp ${ACCESSDIR}string_ext.cpp score.cpp ${ACCESSDIR}SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp ${ACCESSDIR}FastaHandler.cpp SequenceHandlerIsotype.cpp TracebackTable.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp LocalAlignment.cpp SequenceHandlerCommon.cpp ${ACCESSDIR}GzTools.cpp ${ACCESSDIR}FastqHandler.cpp ${ACCESSDIR}FileHandler.cpp ${ACCESSDIR}FASTQ.cpp

SRCS_5=NGSMapping_Demux_main.cpp ${ACCESSDIR}string_ext.cpp score.cpp ${ACCESSDIR}SequenceString.cpp AlignmentString.cpp ${ACCESSDIR}FastaHandler.cpp SequenceHandlerCommon.cpp SequenceHandlerBarcode.cpp ${ACCESSDIR}GzTools.cpp ${ACCESSDIR}FileHandler.cpp ${ACCESSDIR}FastqHandler.cpp ${ACCESSDIR}FASTQ.cpp

SRCS_6=NGSMapping_getBarcode_main.cpp ${ACCESSDIR}FASTQ.cpp ${ACCESSDIR}SequenceString.cpp ${ACCESSDIR}GzTools.cpp ${ACCESSDIR}FastqHandler.cpp ${ACCESSDIR}string_ext.cpp SequenceHandlerBarcode.cpp SequenceHandlerCommon.cpp ${ACCESSDIR}FastaHandler.cpp score.cpp ${ACCESSDIR}FileHandler.cpp

SRCS_7=localAlign_main.cpp ${ACCESSDIR}string_ext.cpp score.cpp ${ACCESSDIR}SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp LocalAlignment.cpp GlobalAlignment.cpp OverlapAlignment.cpp ${ACCESSDIR}FastaHandler.cpp SequenceHandler.cpp GapModel.cpp AffineGapModel.cpp TracebackTable.cpp SequenceHandlerCommon.cpp MarkovChainGapModel_454.cpp ${ACCESSDIR}GzTools.cpp ${ACCESSDIR}FileHandler.cpp ${ACCESSDIR}FASTQ.cpp ${ACCESSDIR}FastqHandler.cpp

SRCS_8=NGS_Concate_main.cpp ${ACCESSDIR}string_ext.cpp ${ACCESSDIR}SequenceString.cpp ${ACCESSDIR}FastaHandler.cpp  ${ACCESSDIR}GzTools.cpp ${ACCESSDIR}FileHandler.cpp ${ACCESSDIR}FASTQ.cpp ${ACCESSDIR}FastqHandler.cpp SequenceHandlerCommon.cpp score.cpp

#SRCS_2=remove_replicates.cpp string_ext.cpp 

OBJS_0=${SRCS_0:.cpp=.o}
OBJS_1=${SRCS_1:.cpp=.o}
OBJS_2=${SRCS_2:.cpp=.o}
OBJS_3=${SRCS_3:.cpp=.o}
OBJS_4=${SRCS_4:.cpp=.o}
OBJS_5=${SRCS_5:.cpp=.o}
OBJS_6=${SRCS_6:.cpp=.o}
OBJS_7=${SRCS_7:.cpp=.o}
OBJS_8=${SRCS_8:.cpp=.o}

#OBJS_1=${SRCS_1:.cpp=.o}
#OBJS_2=${SRCS_2:.cpp=.o}

PROG_0=align
PROG_1=ngsmapping_adaptor
PROG_2=ngsmapping_primer_dimer
PROG_3=ngsmapping_constant
PROG_4=ngsmapping_Isotype
PROG_5=ngsmapping_demux
PROG_6=ngsmapping_getBarcode
PROG_7=localAlign
PROG_8=ngs_concate
#PROG_1=rmxs_at
#PROG_2=remove_replicate
DEPEND=$(GXX) $(CFLAG) -MM

######Rules######

all: $(PROG_0) $(PROG_1) $(PROG_2) $(PROG_3) $(PROG_4) $(PROG_5) $(PROG_6) $(PROG_7) $(PROG_8)

.PHONY: clean all depend

clean:
	rm -fr *.o *~ core $(PROG_0) $(PROG_1) $(PROG_2) $(PROG_3) $(PROG_4) $(PROG_5) $(PROG_6) $(PROG_7) $(PROG_8)
	cd $(ACCESSDIR) && rm -fr *.o *~ core; #ls ;   #pwd ; # note makefile is a script and each line (of commands) is doing its own sub-precess (&& is trying to stop if the "cd" command fails.

.cpp.o:   #old fasion suffix rule, double suffix rule
	$(GXX) $(CXXFLAG) -c $< -o $(addsuffix .o, $(basename $<))

$(PROG_0): $(OBJS_0)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+		
	@echo "******Make complete";echo "";

$(PROG_1): $(OBJS_1)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo "******Make complete"; echo "";
	

$(PROG_2): $(OBJS_2)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo "" 	; echo "******Make complete"

$(PROG_3): $(OBJS_3)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo "";	
	@echo "******Make complete"

$(PROG_4): $(OBJS_4)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo "";	
	@echo "******Make complete"

$(PROG_5): $(OBJS_5)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo "";	
	@echo "******Make complete"

$(PROG_6): $(OBJS_6)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo "";	
	@echo "******Make complete"

$(PROG_7): $(OBJS_7)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo "";	
	@echo "******Make complete"

$(PROG_8): $(OBJS_8)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo "";	
	@echo "******Make complete"
	
#$(PROG_1): $(OBJS_1)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"


#$(PROG_2): $(OBJS_2)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"


depend: .depend

.depend: Makefile $(SRCS_0) $(SRCS_1) $(SRCS_2) $(SRCS_3) $(SRCS_4) $(SRCS_5) $(SRCS_6) $(SRCS_7) $(SRCS_8)
	$(GXX) -MM *.cpp >.depend
	@echo " "
	@echo "****Dependencies generated successfully."


include .depend 

