cmake_minimum_required(VERSION 3.6)
project(EIS)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -L/lib -lnetcdf -lc -lfftw3 -lm")

include_directories ("/usr/include")
link_directories("/usr/lib")

if (WIN32)
    find_package(libnetcdf.dll)
    find_package(libfftw3.dll)
    find_package(libm)
elseif (UNIX)
    find_package(libnetcdf)
    find_package(libfft3)
    find_package(libm)
endif()

set(SOURCE_FILES
        lib/macros.h
        lib/constants.h
        lib/dalloc.c
        lib/dalloc.h
        lib/ialloc.c
        lib/ialloc.h
        lib/structs.h
        lib/slicer.c
        lib/slicer.h
        lib/stats.c
        lib/stats.h
        src/damping.c
        src/damping.h
        src/getdata.c
        src/getdata.h
        src/header.h
        src/lsm.c
        src/lsm.h
        src/solveode.c
        src/solveode.h
        src/save.c
        src/save.h
        src/main.c
        src/divergence.c
        src/divergence.h
        src/wavelength.c
        src/wavelength.h
        src/hybrid.c
        src/hybrid.h
        src/readtxt.c)

add_executable(EIS ${SOURCE_FILES})

target_link_libraries(EIS netcdf fftw3 m)

# set folder to copy executables to
# set(PROJECT_BINARY_DIR ${PROJECT_SOURCE_DIR}/bin)

# get name and location of binary, namecan change with os (cf. ".exe" suffix on windows)
# get_target_property(EIS_BIN_NAME EIS LOCATION)

# copy bin to binary folder
#ADD_CUSTOM_COMMAND(TARGET EIS
#        POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy EIS ${PROJECT_BINARY_DIR}/.)
#