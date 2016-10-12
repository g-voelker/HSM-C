THISDIR = .

LOFILES = alloc_space.o ceigs.o clinalg.o complex.o fft2D.o fft3D.o \
		specialfn.o differentiate.o integrate.o linalg.o \
		myfft.o interpolate.o myerror.o mylaguer.o myzroots.o \
		mywavelet.o polyroots.o transforms.o \
		myzbrent.o scalefld.o matrixOps.o stagger.o statistics.o \
		convert_xy_rt.o get_dmatrix_average.o \
		smoothfunction.o acceptvalues.o \
		readPlot.o writePlot.o readImg.o writeImg.o \
		eigenvalues.o
LCFILES = alloc_space.c ceigs.c clinalg.c complex.c fft2D.c fft3D.c \
		specialfn.c differentiate.c integrate.c linalg.c \
		myfft.c interpolate.c myerror.c mylaguer.c myzroots.c \
		mywavelet.c polyroots.c transforms.c \
		myzbrent.c scalefld.c matrixOps.c stagger.c statistics.c \
		convert_xy_rt.c get_dmatrix_average.c \
		smoothfunction.c acceptvalues.c \
		readPlot.c writePlot.c readImg.c writeImg.c \
		eigenvalues.c

BASEH = complex.h alloc_space.h macros.h constants.h

SCALDFLDH = $(BASEH) interpolate.h differentiate.h

TRANSFORMSH = $(BASEH) specialfn.h

CC = cc 
CFLAGS = -c -Wall
MFLAGS = -lm -limage 
LIBFLAGS = -I/usr/include/malloc -I/opt/local/include
JPEGLIB = /opt/local/lib/libjpeg.a
OPTFLAGS = -O2
CPPFLAGS =
LDFLAGS =
LINTFLAGS = -wlint

.KEEP_STATE:

all: lib

lib: $(LOFILES)

lintlib: $(LOFILES)
	$(CC) $(LINTFLAGS) -c $(LOFILES) $(MFLAGS)  


alloc_space.o: alloc_space.c complex.h
	$(CC) $(OPTFLAGS) $(CFLAGS) alloc_space.c $(LIBFLAGS)

ceigs.o: ceigs.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) ceigs.c $(LIBFLAGS)

clinalg.o: clinalg.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) clinalg.c $(LIBFLAGS)

complex.o: complex.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) complex.c $(LIBFLAGS)

fft2D.o: fft2D.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) fft2D.c $(LIBFLAGS)

fft3D.o: fft3D.c
	$(CC) $(OPTFLAGS) $(CFLAGS) fft3D.c $(LIBFLAGS)
		
specialfn.o: specialfn.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) specialfn.c $(LIBFLAGS)

differentiate.o: differentiate.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) differentiate.c $(LIBFLAGS)

integrate.o: integrate.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) integrate.c $(LIBFLAGS)

linalg.o: linalg.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) linalg.c $(LIBFLAGS)
		
myfft.o: myfft.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) myfft.c $(LIBFLAGS)

interpolate.o: interpolate.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) interpolate.c $(LIBFLAGS)

myerror.o: myerror.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) myerror.c $(LIBFLAGS)

mylaguer.o: mylaguer.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) mylaguer.c $(LIBFLAGS)

myzroots.o: myzroots.c
	$(CC) $(OPTFLAGS) $(CFLAGS) myzroots.c $(LIBFLAGS)

mywavelet.o: mywavelet.c
	$(CC) $(OPTFLAGS) $(CFLAGS) mywavelet.c $(LIBFLAGS)
		
myzbrent.o: myzbrent.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) myzbrent.c $(LIBFLAGS)

polyroots.o: polyroots.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) polyroots.c $(LIBFLAGS)

scalefld.o: scalefld.c $(SCALEFLDH)
	$(CC) $(OPTFLAGS) $(CFLAGS) scalefld.c $(LIBFLAGS)

matrixOps.o: matrixOps.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) matrixOps.c $(LIBFLAGS)

stagger.o: stagger.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) stagger.c $(LIBFLAGS)

statistics.o: statistics.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) statistics.c $(LIBFLAGS)

convert_xy_rt.o: convert_xy_rt.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) convert_xy_rt.c $(LIBFLAGS)

get_dmatrix_average.o: get_dmatrix_average.c
	$(CC) $(OPTFLAGS) $(CFLAGS) get_dmatrix_average.c $(LIBFLAGS)
		
smoothfunction.o: smoothfunction.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) smoothfunction.c $(LIBFLAGS)

acceptvalues.o: acceptvalues.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) acceptvalues.c $(LIBFLAGS)

readPlot.o: readPlot.c 
	$(CC) $(OPTFLAGS) $(CFLAGS) readPlot.c $(LIBFLAGS)

writePlot.o: writePlot.c
	$(CC) $(OPTFLAGS) $(CFLAGS) writePlot.c $(LIBFLAGS)

writeImg.o: writeImg.c
	$(CC) $(OPTFLAGS) $(CFLAGS) writeImg.c $(LIBFLAGS)

readImg.o: readImg.c
	$(CC) $(OPTFLAGS) $(CFLAGS) readImg.c $(LIBFLAGS)

transforms.o: transforms.c $(TRANFORMSH)
	$(CC) $(OPTFLAGS) $(CFLAGS) transforms.c $(LIBFLAGS)

eigenvalues.o: eigenvalues.c $(EIGENVALUESH)
	$(CC) $(OPTFLAGS) $(CFLAGS) eigenvalues.c $(LIBFLAGS)

