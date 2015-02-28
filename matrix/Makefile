###Makefile for LocalAlign
###By Feng 03/08/2014 @bu
#######v1.0

###Define variables######

GCC = gcc
GXX = g++

CFLAG= -Wall -g -Werror -O 
LOADFLAG=-s -lm

CXXFLAG=${CFLAG}

SRCS_0=main.cpp Matrix.cpp
#SRCS_1=	NGSMapping_Adaptor_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp SequenceHandler.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp TracebackTable.cpp SequenceHandlerCommon.cpp LocalAlignment.cpp

#SRCS_2=	NGSMapping_PrimerDimer_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp SequenceHandler.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp TracebackTable.cpp SequenceHandlerCommon.cpp LocalAlignment.cpp

#SRCS_3=	NGSMapping_Constant_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp SequenceHandlerConstant.cpp TracebackTable.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp LocalAlignment.cpp SequenceHandlerCommon.cpp

#SRCS_1=	rm_x_s_probes.cpp string_ext.cpp
#SRCS_2=remove_replicates.cpp string_ext.cpp 

OBJS_0=${SRCS_0:.cpp=.o}
#OBJS_1=${SRCS_1:.cpp=.o}
#OBJS_2=${SRCS_2:.cpp=.o}
#OBJS_3=${SRCS_3:.cpp=.o}
#OBJS_1=${SRCS_1:.cpp=.o}
#OBJS_2=${SRCS_2:.cpp=.o}

PROG_0=test_main
#PROG_1=ngsmapping_adaptor
#PROG_2=ngsmapping_primer_dimer
#PROG_3=ngsmapping_constant
#PROG_1=rmxs_at
#PROG_2=remove_replicate
DEPEND=$(GXX) $(CFLAG) -MM


######Rules######

all: $(PROG_0) #$(PROG_1) $(PROG_2) $(PROG_3) #(PROG_3)

.PHONY: clean all depend

clean:
	rm -fr *.o *~ core $(PROG_0) #$(PROG_1) $(PROG_2) $(PROG_3) ###$(PROG_2) 

.cpp.o:
	$(GXX) $(CXXFLAG) -c $< -o $(addsuffix .o, $(basename $<))

$(PROG_0): $(OBJS_0)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo ""
	@echo "******Make complete:"

#$(PROG_1): $(OBJS_1)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"

#$(PROG_2): $(OBJS_2)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"

#$(PROG_3): $(OBJS_3)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"


#$(PROG_1): $(OBJS_1)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"


#$(PROG_2): $(OBJS_2)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"


depend: .depend

.depend: Makefile $(SRCS_0) #$(SRCS_1) $(SRCS_2) $(SRCS_3)
	$(GXX) -MM *.cpp >.depend
	@echo " "
	@echo "****Dependencies generated successfully."


include .depend 

