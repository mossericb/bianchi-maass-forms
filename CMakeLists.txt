cmake_minimum_required(VERSION 3.25)
project(Bianchi-Maass_Forms)
include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 20)

include_directories(/usr/local/include)
include_directories(/usr/local/lib)
include_directories(/opt/homebrew/include)
include_directories(/opt/homebrew/lib)
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
        KBessel.cpp
        KBessel.h
        SL2C.cpp
        SL2C.h
        BianchiMaassSearch.cpp
        BianchiMaassSearch.h
        ImaginaryQuadraticIntegers.cpp
        ImaginaryQuadraticIntegers.h
        archtKBessel.h
        archtKBessel.cpp
        KBesselReal.cpp
        KBesselReal.h)

target_link_libraries(${PROJECT_NAME} ${FLINT_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${ARB_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${MPFR_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${GMP_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${MPFI_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${FFTW_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${PARI_LIBRARY})

target_link_libraries(${PROJECT_NAME} ${OpenMP})
target_link_libraries(${PROJECT_NAME} ${Boost})

cmake_print_variables(FFTW_LIBRARY)
cmake_print_variables(OpenMP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass=loop-vectorize")

cmake_print_variables(CMAKE_CXX_FLAGS)
cmake_print_variables(CMAKE_C_FLAGS)