###Makefile for sigpig (statistical inference of generation probability of IG)
###By Feng 09/11/2014 @bu
#######v1.0

###Define variables######

GCC = gcc
GXX = g++
MAKE= make

CFLAG= -Wall -g -Werror -O -lpthread
LOADFLAG=-s -lm
AlignmentLibDir=../

CXXFLAG=${CFLAG}


SRCS_0=main.cpp ${AlignmentLibDir}string_ext.cpp ${AlignmentLibDir}score.cpp ${AlignmentLibDir}SequenceString.cpp ${AlignmentLibDir}AlignmentString.cpp ${AlignmentLibDir}AffineGapModel.cpp ${AlignmentLibDir}MarkovChainGapModel_454.cpp ${AlignmentLibDir}TracebackTable.cpp ${AlignmentLibDir}pairwiseAlignment.cpp ${AlignmentLibDir}OverlapAlignment.cpp ${AlignmentLibDir}GapModel.cpp ${AlignmentLibDir}FastaHandler.cpp ${AlignmentLibDir}LocalAlignment.cpp ${AlignmentLibDir}SequenceHandlerCommon.cpp ${AlignmentLibDir}GlobalAlignment.cpp GenomicV.cpp genomicSegments.cpp genomicSegment.cpp GenomicJ.cpp GenomicD.cpp LoadData.cpp

#SRCS_1=	NGSMapping_Adaptor_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp SequenceHandler.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp TracebackTable.cpp SequenceHandlerCommon.cpp

#SRCS_2=	NGSMapping_PrimerDimer_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp SequenceHandler.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp TracebackTable.cpp SequenceHandlerCommon.cpp

#SRCS_3=	NGSMapping_Constant_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp SequenceHandlerConstant.cpp TracebackTable.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp LocalAlignment.cpp SequenceHandlerCommon.cpp

#SRCS_4=	NGSChecking_UPMGCG_ATG_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp TracebackTable.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp LocalAlignment.cpp SequenceHandlerCommon.cpp SequenceHandlerCheckUPM_GCGvsATG.cpp

#SRCS_1=	rm_x_s_probes.cpp string_ext.cpp
#SRCS_2=remove_replicates.cpp string_ext.cpp 

OBJS_0=${SRCS_0:.cpp=.o}
#OBJS_1=${SRCS_1:.cpp=.o}
#OBJS_2=${SRCS_2:.cpp=.o}
#OBJS_3=${SRCS_3:.cpp=.o}
#OBJS_4=${SRCS_4:.cpp=.o}
#OBJS_1=${SRCS_1:.cpp=.o}
#OBJS_2=${SRCS_2:.cpp=.o}

PROG_0=sigpig
#PROG_1=ngsmapping_adaptor
#PROG_2=ngsmapping_primer_dimer
#PROG_3=ngsmapping_constant
#PROG_4=ngschecking_gcg_atg
#PROG_1=rmxs_at
#PROG_2=remove_replicate
DEPEND=$(GXX) $(CFLAG) -MM


######Rules######

all: buildLib $(PROG_0) #$(PROG_1) $(PROG_2) $(PROG_3) $(PROG_4) #(PROG_3)

.PHONY: clean all depend buildLib

buildLib:
	@echo "calling to rebuild the alignment library"
	@cd ../ && $(MAKE) ##this line here is recursively calling another make to compile the alignment lib. The reason to call this is just want to make sure the library is in a good shape on its own and then for this project
	@echo "Done!"


clean:
	rm -fr *.o *~ core $(PROG_0) ###$(PROG_1) $(PROG_2) $(PROG_3) $(PROG_4)###$(PROG_2) 

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


#$(PROG_1): $(OBJS_1)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"


#$(PROG_2): $(OBJS_2)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"


depend: .depend

.depend: Makefile $(SRCS_0) #$(SRCS_1) $(SRCS_2) $(SRCS_3) $(SRCS_4)
	$(GXX) -MM *.cpp >.depend
	@echo " "
	@echo "****Dependencies generated successfully."


include .depend 

