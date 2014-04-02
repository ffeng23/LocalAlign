###Makefile for LocalAlign
###By Feng 03/08/2014 @bu
#######v1.0

###Define variables######

GCC = gcc
GXX = g++

CFLAG= -Wall -g -Werror -O 
LOADFLAG=-s -lm

CXXFLAG=${CFLAG}

SRCS_0=main.cpp string_ext.cpp score.cpp SequenceString.cpp AlignmentString.cpp pairwiseAlignment.cpp LocalAlignment.cpp GlobalAlignment.cpp OverlapAlignment.cpp FastaHandler.cpp
#SRCS_1=	rm_x_s_probes.cpp string_ext.cpp
#SRCS_2=remove_replicates.cpp string_ext.cpp 

OBJS_0=${SRCS_0:.cpp=.o}
#OBJS_1=${SRCS_1:.cpp=.o}
#OBJS_2=${SRCS_2:.cpp=.o}

PROG_0=testalign
#PROG_1=rmxs_at
#PROG_2=remove_replicate
DEPEND=$(GXX) $(CFLAG) -MM


######Rules######

all: $(PROG_0) #$(PROG_1) #$(PROG_2)

.PHONY: clean all depend

clean:
	rm -fr *.o *~ core $(PROG_0) $(PROG_1) ###$(PROG_2) 

.cpp.o:
	$(GXX) $(CXXFLAG) -c $< -o $(addsuffix .o, $(basename $<))

$(PROG_0): $(OBJS_0)
	$(GXX) -o $@ $(CXXFLAG) $(LOADFLAG) $+
	@echo ""
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

.depend: Makefile $(SRCS_0) #$(SRCS_1) $(SRCS_2)
	$(GXX) -MM *.cpp >.depend
	@echo " "
	@echo "****Dependencies generated successfully."


include .depend 
