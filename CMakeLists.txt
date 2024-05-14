cmake_minimum_required(VERSION 3.27)
project(CWT)

set(CMAKE_CXX_STANDARD 17)

add_executable(CWT benchmark/benchmark.cpp
        transforms.cpp
        transforms.h
        Scales.cpp
        Scales.h
        Ricker.cpp
        Ricker.h)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
if (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(CWT PRIVATE "$<$<CONFIG:Release>:/Zi>")
    target_link_options(CWT PRIVATE "$<$<CONFIG:Release>:/DEBUG>")
    target_link_options(CWT PRIVATE "$<$<CONFIG:Release>:/OPT:REF>")
    target_link_options(CWT PRIVATE "$<$<CONFIG:Release>:/OPT:ICF>")
endif()