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

find_library(FLINT_LIBRARY flint REQUIRED)
find_library(ARB_LIBRARY arb REQUIRED)
find_library(MPFR_LIBRARY mpfr REQUIRED)
find_library(GMP_LIBRARY gmp REQUIRED)
find_library(MPFI_LIBRARY mpfi REQUIRED)

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED)

if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(${PROJECT_NAME}
        src/main.cpp
        src/Index.cpp
        include/Index.h
        src/Quaternion.cpp
        include/Quaternion.h
        src/Auxiliary.cpp
        include/Auxiliary.h
        src/KBessel.cpp
        include/KBessel.h
        src/SL2C.cpp
        include/SL2C.h
        src/BianchiMaassSearch.cpp
        include/BianchiMaassSearch.h
        src/ImaginaryQuadraticIntegers.cpp
        include/ImaginaryQuadraticIntegers.h
        include/archtKBessel.h
        src/archtKBessel.cpp
        src/KBesselReal.cpp
        include/KBesselReal.h)

target_link_libraries(${PROJECT_NAME} ${FLINT_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${ARB_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${MPFR_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${GMP_LIBRARY})
target_link_libraries(${PROJECT_NAME} ${MPFI_LIBRARY})

target_link_libraries(${PROJECT_NAME} ${OpenMP})
target_link_libraries(${PROJECT_NAME} ${Boost})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass=loop-vectorize")

cmake_print_variables(CMAKE_CXX_FLAGS)
cmake_print_variables(CMAKE_C_FLAGS)