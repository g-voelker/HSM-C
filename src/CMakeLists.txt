cmake_minimum_required(VERSION 3.5)
project(src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    damping.c
    getdata.c
    header.h
    input.h
    lsm.c
    lsm.h
    main.c)

add_executable(src ${SOURCE_FILES})