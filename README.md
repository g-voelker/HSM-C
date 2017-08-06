# HSM-C - The Hybrid Slab Model in C

## Introduction

   HSM-C is a scientific model that simulates energy fluxes associated with the excitation of near-inertial internal gravity waves by wind stress following Voelker et. al. (2017). The Software is distributed under the GPL 3.0 license and my be copied, redistributed or modified accordingly. A copy of the [license](GPL3-LICENSE) can be found as part of the software repository.


## Obtaining and Building from Source

   The actual stable source can be obtained from a git repository found [here](https://g-voelker@bitbucket.org/g-voelker/hsm-c.git)

   All of the following dependencies are freely available and should be part of the toolchain

   - cmake (ver. 3.6)
   - the FFTW3 library (ver. 3.3.6)
   - the NETCDF-C libary (ver. 4.4.1.1)
   - the libm (math) c library

   The versions indicated in brackets are the stable version at the development time of the code - other versions might work but are not tested to date.
   To compile the code you need to clone the code (or download and unpack), navigate to the root directory and run

```cmake CMakeLists.txt

make Makefile```

   Alternatively, an IDE like Visual Studio, Eclipse or CLion can be used. On linux machines the compilation will be straight forward as c compilers are commonly part of the distribution or easily accessible. On windows / Mac a compiler and the libraries will have to be installed manually and communicate accordingly. I recommend to have a look at cygwin (Win) / macports (Mac) to collect the requirements.


## Execution of the Simulation

   To execute the simulation the parameter file "input.txt" that has to be in the same directory as the binary and tailored to the desired setup. See input.txt for documentation. The model assumes [NCEP-CFSR](http://cfs.ncep.noaa.gov/cfsr/) wind stress data and accordingly regridded [MIMOC](https://www.pmel.noaa.gov/mimoc/) data for execution. I prepared my data using the grib_to_netcdf routine deliered with the GRIB2 library and the nccopy tool from the NETCDF Software package. For more information on data conversion please contact voelker@uni-bremen.de.

   Run the code from the command line with (UNIX)

```./HSM-C```

   or (Windows)

```./HSM-C.exe```

## Modification of the Model

   Feel free to modify the source according to the above mentioned license. All code is documented / commented and checked with valgrind against memory leaks and access to unallocated space. If you find any bug please report at the repository.
   To track the recent changes made to the stable source checkout the recent commits [here](https://bitbucket.org/g-voelker/hsm-c/commits/all)