cmake_minimum_required(VERSION 3.25)
project(Bianchi-Maass_Forms)
include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 20)

#/opt/homebrew/Cellar/llvm/16.0.6/bin/clang

include_directories(/usr/local/include)
include_directories(/usr/local/lib)
include_directories(/opt/homebrew/include)
include_directories(/opt/homebrew/lib)
#include_directories(/Users/ericmoss/sage/sage-10.0/local/include)
#include_directories(/Users/ericmoss/sage/sage-10.0/local/lib)
include_directories(/opt/homebrew/Cellar/libomp/16.0.6/lib)
include_directories(/opt/homebrew/Cellar/libomp/16.0.6/include)
include_directories(/Users/ericmoss/opt/anaconda3/pkgs/python-3.9.13-hdfd78df_1/lib)

find_library(FLINT_LIBRARY flint REQUIRED)
find_library(ARB_LIBRARY arb REQUIRED)
find_library(MPFR_LIBRARY mpfr REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(MPFI_LIBRARY mpfi REQUIRED)
find_library(MPI_LIBRARY mpi REQUIRED)
find_library(PARI_LIBRARY pari REQUIRED)
find_library(FFTW_LIBRARY fftw3 REQUIRED)
#find_library(ARCHTKBESSEL_LIBRARY ArchtKBessel REQUIRED)
#find_library(ARCHT_LIBRARY ArchtKBessel REQUIRED)

#variable_watch(ARCHT_LIBRARY)
variable_watch(PYTHON_LIBRARY)

#find_package(FFTW CONFIG REQUIRED)
#find_package(OpenCV REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(${PROJECT_NAME}
        main.cpp
        Index.cpp
        Index.h
        Quaternion.cpp
        Quaternion.h
        Auxiliary.cpp
        Auxiliary.h
        KBesselApproximator.cpp
        KBesselApproximator.h
        SL2C.cpp
        SL2C.h
        FunctionToEvaluate.h
        FunctionToEvaluate.cpp
        #[[PlotWindow.h
        PlotWindow.cpp
        Plotter.h
        Plotter.cpp]]
        tinycolormap.hpp
        BianchiMaassPointwise.cpp
        BianchiMaassPointwise.h
        BianchiMaassSearch.cpp
        BianchiMaassSearch.h
        ImaginaryQuadraticIntegers.cpp
        ImaginaryQuadraticIntegers.h
        KBesselExact.cpp
        KBesselExact.h
        ArchtKBessel.h
        ArchtKBessel.cpp)

target_link_libraries(${PROJECT_NAME} ${FLINT_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${ARB_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${MPFR_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${GMP_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${MPFI_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${FFTW_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${PARI_LIBRARY})
#target_link_libraries(${PROJECT_NAME} ${ARCHTKBESSEL_LIBRARY})

target_link_libraries(${PROJECT_NAME} ${OpenMP})
#target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS})
target_link_libraries(${PROJECT_NAME} ${Boost})

cmake_print_variables(FFTW_LIBRARY)
cmake_print_variables(OpenMP)
#cmake_print_variables(ARCHTKBESSEL_LIBRARY)

target_link_libraries(${PROJECT_NAME} ${ARCHT_LIBRARY})
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mno-omit-leaf-frame-pointer")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

cmake_print_variables(CMAKE_CXX_FLAGS)
cmake_print_variables(CMAKE_C_FLAGS)

#[[
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I /opt/homebrew/Cellar/python@3.11/3.11.5/Frameworks/Python.framework/Versions/3.11/include/python3.11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I /opt/homebrew/Cellar/numpy/1.26.0/lib/python3.11/site-packages/numpy/core/include")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L /opt/homebrew/Cellar/python@3.11/3.11.5/Frameworks/Python.framework/Versions/3.11/lib")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpython3.11")]]
