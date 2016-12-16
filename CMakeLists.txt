cmake_minimum_required(VERSION 3.6)
project(EIS)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lnetcdf -lc -lfftw3 -lm")
# include_directories ("/usr/include")

set(SOURCE_FILES
        lib/macros.h
        lib/constants.h
        lib/dalloc.c
        lib/dalloc.h
        lib/ialloc.c
        lib/ialloc.h
        lib/structs.h
        src/damping.c
        src/damping.h
        src/getdata.c
        src/getdata.h
        src/header.h
        src/input.h
        src/lsm.c
        src/lsm.h
        src/solveode.c
        src/solveode.h
        src/save.c
        src/save.h
        src/main.c lib/slicer.c lib/slicer.h src/divergence.c src/divergence.h src/wavelength.c src/wavelength.h lib/stats.c lib/stats.h)

add_executable(EIS ${SOURCE_FILES})

# set folder to copy executables to
# set(PROJECT_BINARY_DIR ${PROJECT_SOURCE_DIR}/bin)

# get name and location of binary, namecan change with os (cf. ".exe" suffix on windows)
# get_target_property(EIS_BIN_NAME EIS LOCATION)

# copy bin to binary folder
#ADD_CUSTOM_COMMAND(TARGET EIS
#        POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy EIS ${PROJECT_BINARY_DIR}/.)
#