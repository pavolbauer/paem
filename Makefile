CC     = gcc
###########################################################################
# this file is normally executed by the MATLAB interface
###########################################################################

ifeq ($(URDME_ROOT),)
    $(error URDME_ROOT is not set)
endif
URDME=$(URDME_ROOT)
LIB = $(URDME)/lib
BIN = $(URDME)/bin
INCLUDE = $(URDME)/include
SRC = $(URDME)/src

###########################################################################
MATLAB = $(shell $(URDME_ROOT)/bin/urdme_init -m)
MATLAB_ARCH = $(shell $(URDME_ROOT)/bin/urdme_init -a)
###########################################################################
MATLAB_EXT = $(MATLAB)/extern
MATLAB_INC=$(MATLAB_EXT)/include
###########################################################################
#ARCHS=$(shell $(URDME_ROOT)/bin/urdme_init -s)
ifneq ($(ARCHS),)
SET_ARCH=-arch $(ARCHS)
endif
###########################################################################

URDME_LIBMAT=YES
ifeq ($(URDME_LIBMAT),)
LMAT = -lmx -lmat -lmex -lmwservices 
MATLAB_LIB =-L$(MATLAB)/bin/$(MATLAB_ARCH)
else
DMAT = -DURDME_LIBMAT
endif

# Sparse output
#DSPARSE = -DURDME_OUTPUT_SPARSE
ifdef DSPARSE
ifdef URDME_LIBMAT
$(error URDME_LIBMAT does not support sparse output, please use Matlab libraries.)
endif
endif

GLIB = /usr/include/glib-2.0/

CFLAGS = -c $(SET_ARCH) -O0 -g -std=gnu99 $(MEMOPTS) -I$(MATLAB_LIB) -I$(MATLAB_INC) -I$(INCLUDE) $(DMAT) $(DSPARSE) `pkg-config --cflags --libs glib-2.0 gsl` -pthread -D_GNU_SOURCE -lm

LFLAGS = `pkg-config --libs glib-2.0 gsl` -pthread -lm
ifneq ($(MATLAB_LIB),)
LFLAGS += $(LMAT) $(MATLAB_LIB) 
endif

###########################################################################
ifeq ($(URDME_MODEL),)
     $(error URDME_MODEL not set)
endif
###########################################################################
# destination directory
OUT=./
###########################################################################
###########################################################################
###########################################################################
###########################################################################

all:  aem

aem: $(OUT)$(URDME_MODEL).aem

.PHONY: all clean 

###########################################################################

DEPENDS = -MT $@ -MD -MP -MF $(subst .o,.d,$@)

SOURCES= $(URDME_MODEL).c $(wildcard $(SRC)/*.c) $(wildcard $(SRC)/aem/*.c) $(SRC)/nsm/binheap.c
OBJS=$(SOURCES:.c=.o)

#ifeq ($(URDME_LIBMAT),)
#OBJS=nsm.o report.o matmodel.o $(URDME_MODEL).o
#else
#OBJS=nsm.o nsmcore.o report.o matmodel.o read_matfile.o $(URDME_MODEL).o
#endif


%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

testmakefile: 
	-echo $(URDME_MODEL).aem
	-echo $(URDME)


$(OUT)$(URDME_MODEL).aem: $(OBJS)
	$(CC) $(DEPENDS) $(OBJS) -o $(OUT)$(URDME_MODEL).aem $(LFLAGS)

$(OUT):
	-mkdir -p $(OUT)


clean:
	rm -f $(OBJS) $(OUT)$(URDME_MODEL).aem


