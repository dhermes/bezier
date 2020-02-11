cmake_minimum_required(VERSION 3.1)

project(
  quadpack
  DESCRIPTION "Library for numerical integration of one-dimensional functions."
  LANGUAGES Fortran)

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Intel|Flang)$")
  message(FATAL_ERROR "gfortran, ifort and flang are the only supported compilers (current compiler ID is ${CMAKE_Fortran_COMPILER_ID})")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release)$")
  message(FATAL_ERROR "Debug and Release are the only supported build types")
endif()

# Unset all Fortran flags (instead rely on our flags).
set(OLD_CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS})
set(OLD_CMAKE_Fortran_FLAGS_DEBUG ${CMAKE_Fortran_FLAGS_DEBUG})
set(OLD_CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
set(CMAKE_Fortran_FLAGS "")
set(CMAKE_Fortran_FLAGS_DEBUG "")
set(CMAKE_Fortran_FLAGS_RELEASE "")

file(GLOB SOURCES *.f)
add_library(quadpack OBJECT ${SOURCES})
unset(SOURCES)

set_target_properties(quadpack PROPERTIES Fortran_FORMAT "FIXED")
if(${BUILD_SHARED_LIBS})
  set_target_properties(quadpack PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "^GNU$")
  # ``-Wextra`` includes ``no-compare-reals``, but we have comparisons like
  # ``value == 0.0_dp``.
  target_compile_options(
    quadpack
    PRIVATE -fno-second-underscore
            -Wall
            -Wextra
            -Wno-compare-reals
            -Wimplicit-interface
            -fmax-errors=1
            -std=f2008)
  if(CMAKE_BUILD_TYPE MATCHES "^Debug$")
    target_compile_options(quadpack PRIVATE -fcheck=all -fbacktrace
                                            -fimplicit-none -pedantic)
    target_compile_options(quadpack PUBLIC -g)
  elseif(CMAKE_BUILD_TYPE MATCHES "^Release$")
    if(${TARGET_NATIVE_ARCH})
      target_compile_options(quadpack PRIVATE -march=native)
    endif()
    target_compile_options(quadpack PUBLIC -DNDEBUG -O3 -funroll-loops)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel$")
  target_compile_options(quadpack PUBLIC -warn all)
  if(CMAKE_BUILD_TYPE MATCHES "^Debug$")
    target_compile_options(quadpack PUBLIC -check all)
  elseif(CMAKE_BUILD_TYPE MATCHES "^Release$")
    target_compile_options(quadpack PUBLIC -fast)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^Flang$")
  target_compile_options(
    quadpack
    PRIVATE -Wall
            -Wextra
            -Wimplicit
            -ferror-limit=1
            -std=f2008)
  if(CMAKE_BUILD_TYPE MATCHES "^Debug$")
    target_compile_options(quadpack PRIVATE -fcheck=all -fbacktrace
                                            -fimplicit-none -pedantic)
    target_compile_options(quadpack PUBLIC -g)
  elseif(CMAKE_BUILD_TYPE MATCHES "^Release$")
    if(${TARGET_NATIVE_ARCH})
      target_compile_options(quadpack PRIVATE -march=native)
    endif()
    target_compile_options(quadpack PUBLIC -DNDEBUG -O3 -funroll-loops)
  endif()
else()
  message(FATAL_ERROR "gfortran, ifort and flang are the only supported compilers (current compiler ID is ${CMAKE_Fortran_COMPILER_ID})")
endif()

# Restore all Fortran flags.
set(CMAKE_Fortran_FLAGS ${OLD_CMAKE_Fortran_FLAGS})
set(CMAKE_Fortran_FLAGS_DEBUG ${OLD_CMAKE_Fortran_FLAGS_DEBUG})
set(CMAKE_Fortran_FLAGS_RELEASE ${OLD_CMAKE_Fortran_FLAGS_RELEASE})
unset(OLD_CMAKE_Fortran_FLAGS)
unset(OLD_CMAKE_Fortran_FLAGS_DEBUG)
unset(OLD_CMAKE_Fortran_FLAGS_RELEASE)
