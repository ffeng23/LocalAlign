###Makefile for sigpig (statistical inference of generation probability of IG)
###By Feng 09/11/2014 @bu
#######v1.0

###Define variables######
##	@echo "doing making here"

GCC = gcc
GXX = g++
MAKE= make

CFLAG=-g -Werror -O -I ../Lib/include
LOADFLAG=-s -lm  -lpthread
#LOADFLAG=-s -lm -L../Lib/lib/ -lpthread -lboost_serialization
#  the -L ../Lib/lib is the searching dir for shared libraries. -lboost_serialization is the boost serilization library, but this is not what we do now. so remove it for now.
AlignmentLibDir=../

CXXFLAG=${CFLAG}

SRCS_0=main.cpp ${AlignmentLibDir}string_ext.cpp ${AlignmentLibDir}score.cpp ${AlignmentLibDir}SequenceString.cpp ${AlignmentLibDir}AlignmentString.cpp ${AlignmentLibDir}AffineGapModel.cpp ${AlignmentLibDir}MarkovChainGapModel_454.cpp ${AlignmentLibDir}TracebackTable.cpp ${AlignmentLibDir}pairwiseAlignment.cpp ${AlignmentLibDir}OverlapAlignment.cpp ${AlignmentLibDir}GapModel.cpp ${AlignmentLibDir}FastaHandler.cpp ${AlignmentLibDir}LocalAlignment.cpp ${AlignmentLibDir}SequenceHandlerCommon.cpp ${AlignmentLibDir}GlobalAlignment.cpp GenomicV.cpp genomicSegments.cpp genomicSegment.cpp GenomicJ.cpp GenomicD.cpp LoadData.cpp MatrixFunctions.cpp Alignment.cpp AlignmentSettings.cpp Alignment_V.cpp Alignment_D.cpp do_VDJ_alignment.cpp DataIO.cpp

SRCS_1=tester_main.cpp MatrixFunctions.cpp Alignment.cpp genomicSegments.cpp ${AlignmentLibDir}string_ext.cpp ${AlignmentLibDir}SequenceString.cpp genomicSegment.cpp GenomicV.cpp GenomicD.cpp GenomicJ.cpp ${AlignmentLibDir}FastaHandler.cpp AlignmentSettings.cpp Alignment_D.cpp ${AlignmentLibDir}score.cpp ${AlignmentLibDir}AlignmentString.cpp ${AlignmentLibDir}LocalAlignment.cpp ${AlignmentLibDir}pairwiseAlignment.cpp ${AlignmentLibDir}TracebackTable.cpp ${AlignmentLibDir}AffineGapModel.cpp ${AlignmentLibDir}MarkovChainGapModel_454.cpp ${AlignmentLibDir}GapModel.cpp


SRCS_2=tester_main_matchV.cpp MatrixFunctions.cpp Alignment.cpp genomicSegments.cpp ${AlignmentLibDir}string_ext.cpp ${AlignmentLibDir}SequenceString.cpp genomicSegment.cpp GenomicV.cpp GenomicD.cpp GenomicJ.cpp ${AlignmentLibDir}FastaHandler.cpp AlignmentSettings.cpp Alignment_V.cpp Alignment_D.cpp ${AlignmentLibDir}score.cpp ${AlignmentLibDir}AlignmentString.cpp ${AlignmentLibDir}LocalAlignment.cpp ${AlignmentLibDir}pairwiseAlignment.cpp ${AlignmentLibDir}TracebackTable.cpp ${AlignmentLibDir}AffineGapModel.cpp ${AlignmentLibDir}MarkovChainGapModel_454.cpp ${AlignmentLibDir}GapModel.cpp

SRCS_3=tester_main_matchD.cpp MatrixFunctions.cpp Alignment.cpp genomicSegments.cpp ${AlignmentLibDir}string_ext.cpp ${AlignmentLibDir}SequenceString.cpp genomicSegment.cpp GenomicV.cpp GenomicD.cpp GenomicJ.cpp ${AlignmentLibDir}FastaHandler.cpp AlignmentSettings.cpp Alignment_V.cpp Alignment_D.cpp ${AlignmentLibDir}score.cpp ${AlignmentLibDir}AlignmentString.cpp ${AlignmentLibDir}LocalAlignment.cpp ${AlignmentLibDir}pairwiseAlignment.cpp ${AlignmentLibDir}TracebackTable.cpp ${AlignmentLibDir}AffineGapModel.cpp ${AlignmentLibDir}MarkovChainGapModel_454.cpp ${AlignmentLibDir}GapModel.cpp


#SRCS_3=tester_serialization.cpp


#SRCS_2=	NGSMapping_PrimerDimer_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp SequenceHandler.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp TracebackTable.cpp SequenceHandlerCommon.cpp

#SRCS_3=	NGSMapping_Constant_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp SequenceHandlerConstant.cpp TracebackTable.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp LocalAlignment.cpp SequenceHandlerCommon.cpp

#SRCS_4=	NGSChecking_UPMGCG_ATG_main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp TracebackTable.cpp GapModel.cpp AffineGapModel.cpp MarkovChainGapModel_454.cpp LocalAlignment.cpp SequenceHandlerCommon.cpp SequenceHandlerCheckUPM_GCGvsATG.cpp

#SRCS_1=	rm_x_s_probes.cpp string_ext.cpp
#SRCS_2=remove_replicates.cpp string_ext.cpp 

OBJS_0=${SRCS_0:.cpp=.o}
OBJS_1=${SRCS_1:.cpp=.o}
OBJS_2=${SRCS_2:.cpp=.o}
OBJS_3=${SRCS_3:.cpp=.o}

#OBJS_3=${SRCS_3:.cpp=.o}
#OBJS_4=${SRCS_4:.cpp=.o}
#OBJS_1=${SRCS_1:.cpp=.o}
#OBJS_2=${SRCS_2:.cpp=.o}

PROG_0=sigpig
PROG_1=tester
PROG_2=tester_matchV
PROG_3=tester_matchD
#PROG_3=tester_serialization
#PROG_3=ngsmapping_constant
#PROG_4=ngschecking_gcg_atg
#PROG_1=rmxs_at
#PROG_2=remove_replicate
DEPEND=$(GXX) $(CFLAG) -MM


######Rules######

all: buildLib $(PROG_0) $(PROG_1) $(PROG_2) $(PROG_3) ####$(PROG_4) #$(PROG_3)

.PHONY: clean all depend buildLib

buildLib:
	@echo "calling to rebuild the alignment library"
	@cd ../ && $(MAKE) ##this line here is recursively calling another make to compile the alignment lib. The reason to call this is just want to make sure the library is in a good shape on its own and then for this project
	@echo "Done!"


clean:
	rm -fr *.o *~ core $(PROG_0) $(PROG_1) $(PROG_2) ##$(PROG_3) ###$(PROG_4)###$(PROG_2) 

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

$(PROG_2): $(OBJS_2)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo ""
	@echo "******Make complete"
	@echo ""

$(PROG_3): $(OBJS_3)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo ""
	@echo "******Make complete"
	@echo ""



#$(PROG_1): $(OBJS_1)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"


#$(PROG_2): $(OBJS_2)
#	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
#	@echo ""
#	@echo "******Make complete"


depend: .depend

.depend: Makefile $(SRCS_0) $(SRCS_1) $(SRCS_2) $(SRCS_3) #$(SRCS_4)
	$(GXX) -MM *.cpp >.depend
	@echo " "
	@echo "****Dependencies generated successfully."


include .depend 

