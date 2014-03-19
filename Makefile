# C Compiler
CC = gcc

# Various machine and compiler specific flags
#    Select one set.

# LINUX Optimizing:
CFLAGS = -W -Wall -O3 -msse -msse2 -malign-double
# LINUX Debugging:
#CFLAGS = -Wall -W -ansi -pedantic -ggdb -pg -I. -O2 -mpentiumpro -malign-double

# Use the following if we are using DMALLOC
#DMALLOCLINK = -ldmalloc
#CFLAGS += -DUSEDMALLOC

# Note:  If you want to use a machine-optimized BLAS, the following
#        lines muct be uncommented and corrected.
#        You should really test the timings!  This may actually
#        be slower than the non BLAS version.
# Path to the BLAS libraries
#BLASLIBDIR = /usr/local/lib
# How to link with the PGPLOT libs
#BLASLINK = -L$(BLASLIBDIR) -lblas
# Adjustments to CFLAGS to use BLAS
#CFLAGS += -DUSE_BLAS

# Objects 
OBJS = eventdft.o eventdft_cmd.o period.c utils.o stats.o \
	dcdflib.o ipmpar.o

# Set the Date
DATE = $(shell date +%d%b%y)

all:  eventdft

eventdft: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(BLASLINK) $(DMALLOCLINK) -lm
clig:
	clig -d eventdft_cmd.cli

psman:
	nroff -man eventdft.1 | a2x -man > eventdft_man.ps

clean:
	rm -f *.o *~ *#

tar: squeaky package

squeaky:
	rm -f eventdft *.out *.o *~ *.tgz *#

package:
	cd ..; tar --exclude CVS -czf eventdft$(DATE).tgz event_DFT
	mv ../eventdft$(DATE).tgz .










