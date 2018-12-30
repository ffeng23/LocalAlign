
###Makefile for LocalAlign
###By Feng 03/08/2014 @bu
#######v1.0

###Define variables######

GCC = gcc
GXX = g++

CFLAG= -Wall -g -Werror -O 
LOADFLAG=-s -lm -lz

CXXFLAG=${CFLAG}
ALIGNMENTDIR=../
SIGPIGDIR=../SIGPIG/
SRCS_0=main.cpp Poisson.cpp 

SRCS_1=	main_fastq.cpp ${ALIGNMENTDIR}FASTQ.cpp ${ALIGNMENTDIR}SequenceString.cpp GzTools.cpp ${ALIGNMENTDIR}FastqHandler.cpp string_ext.cpp ${ALIGNMENTDIR}SequenceHandlerBarcode.cpp ${ALIGNMENTDIR}SequenceHandlerCommon.cpp ${ALIGNMENTDIR}FastaHandler.cpp ${ALIGNMENTDIR}score.cpp

#SRCS_1=	main_reformat.cpp ${ALIGNMENTDIR}string_ext.cpp ${ALIGNMENTDIR}SequenceString.cpp ${ALIGNMENTDIR}FastaHandler.cpp ${SIGPIGDIR}LoadData.cpp

#SRCS_2=	main_showAlignments.cpp ${ALIGNMENTDIR}string_ext.cpp ${ALIGNMENTDIR}score.cpp ${ALIGNMENTDIR}SequenceString.cpp ${ALIGNMENTDIR}AlignmentString.cpp ${ALIGNMENTDIR}pairwiseAlignment.cpp ${ALIGNMENTDIR}GapModel.cpp ${ALIGNMENTDIR}AffineGapModel.cpp ${ALIGNMENTDIR}MarkovChainGapModel_454.cpp ${ALIGNMENTDIR}TracebackTable.cpp ${ALIGNMENTDIR}LocalAlignment.cpp ${SIGPIGDIR}GenomicJ.cpp ${SIGPIGDIR}GenomicD.cpp ${SIGPIGDIR}GenomicV.cpp ${SIGPIGDIR}genomicSegments.cpp ${SIGPIGDIR}genomicSegment.cpp ${SIGPIGDIR}AlignmentSettings.cpp ${SIGPIGDIR}Alignment.cpp ${SIGPIGDIR}Alignment_V.cpp ${SIGPIGDIR}Alignment_D.cpp ${SIGPIGDIR}DataIO.cpp ${SIGPIGDIR}MatrixFunctions.cpp ${ALIGNMENTDIR}FastaHandler.cpp

#SRCS_3=	NGSMapping_Constant_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp SequenceHandlerConstant.cpp TracebackTable.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp LocalAlignment.cpp SequenceHandlerCommon.cpp

OBJS_0=${SRCS_0:.cpp=.o}
OBJS_1=${SRCS_1:.cpp=.o}
#OBJS_2=${SRCS_2:.cpp=.o}
#OBJS_3=${SRCS_3:.cpp=.o}

PROG_0=test
PROG_1=test_fastq
#PROG_2=show_alignment
#PROG_3=ngsmapping_constant


DEPEND=$(GXX) $(CFLAG) -MM


######Rules######

all: $(PROG_0) $(PROG_1) #$(PROG_2) #$(PROG_3) 

.PHONY: clean all depend

clean:
	rm -fr *.o *~ core $(PROG_0) $(PROG_1) #$(PROG_2) $(PROG_3) ###$(PROG_2) 

.cpp.o:
	$(GXX) $(CXXFLAG) -c $< -o $(addsuffix .o, $(basename $<))

$(PROG_0): $(OBJS_0)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo ""
	@echo "******Make complete:"

$(PROG_1): $(OBJS_1)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo ""
	@echo "******Make complete"

#$(PROG_2): $(OBJS_2)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"

#$(PROG_3): $(OBJS_3)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"

depend: .depend

.depend: Makefile $(SRCS_0) $(SRCS_1) #$(SRCS_2) #$(SRCS_3)
	$(GXX) -MM *.cpp >.depend
	@echo " "
	@echo "****Dependencies generated successfully."


include .depend 

