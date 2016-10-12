THISDIR = src
LIBDIR = ../../../lib

ALLOC_SPACE = $(LIBDIR)/alloc_space
COMPLEX = $(LIBDIR)/complex

INCL	= /usr/include
NETCDF = $(INCL)/netcdf


AOFILES = main.o \
	lsm.o \
	$(ALLOC_SPACE).o $(COMPLEX).o
ACFILES = main.c  \
	lsm.c \
	$(ALLOC_SPACE).c $(COMPLEX).c

BASEH = $(COMPLEX).h $(ALLOC_SPACE).h 

MAINH = $(BASEH) lsm.h

LSMH	= $(BASEH) $(NETCDF).h

CC = cc
CFLAGS = -c
MFLAGS = -lm
LIBFLAGS = -I$(LIBDIR) -I$(INCL) -lnetcdff -lnetcdf
OPTFLAGS = -O -Wall
CPPFLAGS =
LDFLAGS =
LINTFLAGS =

.KEEP_STATE:

all: eis3.1

clean: rm -f *.o

eis3.1: $(AOFILES)
	$(CC) $(OPTFLAGS) $(AOFILES) -o eis3.1 $(LIBFLAGS) $(MFLAGS)

main.o: main.c $(MAINH)
	$(CC) $(OPTFLAGS) $(CFLAGS) main.c $(LIBFLAGS)

lsm.o: lsm.c $(LSMH)
	$(CC) $(OPTFLAGS) $(CFLAGS) lsm.c $(LIBFLAGS)



