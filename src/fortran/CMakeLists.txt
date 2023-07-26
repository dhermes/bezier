cmake_minimum_required(VERSION 3.1)

project(
  bezier
  VERSION 2023.7.26
  LANGUAGES Fortran C)
set(AUTHOR "Danny Hermes")
set(AUTHOR_DETAILS "daniel.j.hermes@gmail.com")
set(DESCRIPTION "Library for Bezier curves and triangles.")

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(TARGET_NATIVE_ARCH "Optimize build for host (native) architecture" ON)

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Intel)$")
  message(
    FATAL_ERROR
      "gfortran and ifort are the only supported compilers (current compiler ID is ${CMAKE_Fortran_COMPILER_ID})"
  )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(WIN32 AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static")

  if(NOT CMAKE_IMPORT_LIBRARY_SUFFIX_Fortran)
    set(CMAKE_IMPORT_LIBRARY_SUFFIX_Fortran ".lib")
  endif()
  if(NOT CMAKE_IMPORT_LIBRARY_PREFIX_Fortran)
    set(CMAKE_IMPORT_LIBRARY_PREFIX_Fortran "")
  endif()
  if(NOT CMAKE_SHARED_LIBRARY_PREFIX_Fortran)
    set(CMAKE_SHARED_LIBRARY_PREFIX_Fortran "")
  endif()
endif()

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release)$")
  message(
    FATAL_ERROR
      "Debug and Release are the only supported build types (current build type is ${CMAKE_BUILD_TYPE})"
  )
endif()

# Unset all Fortran flags (instead rely on our flags).
set(OLD_CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS})
set(OLD_CMAKE_Fortran_FLAGS_DEBUG ${CMAKE_Fortran_FLAGS_DEBUG})
set(OLD_CMAKE_Fortran_FLAGS_RELEASE ${CMAKE_Fortran_FLAGS_RELEASE})
set(CMAKE_Fortran_FLAGS "")
set(CMAKE_Fortran_FLAGS_DEBUG "")
set(CMAKE_Fortran_FLAGS_RELEASE "")

file(GLOB SOURCES *.f90)
add_library(bezier ${SOURCES})
unset(SOURCES)

if(${BUILD_SHARED_LIBS})
  set_target_properties(
    bezier
    PROPERTIES VERSION ${bezier_VERSION}
               SOVERSION 2023
               POSITION_INDEPENDENT_CODE ON)
endif()

target_include_directories(bezier PUBLIC $<INSTALL_INTERFACE:include>)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # ``-Wextra`` includes ``no-compare-reals``, but we have comparisons like
  # ``value == 0.0_dp``.
  target_compile_options(
    bezier
    PRIVATE -fno-second-underscore
            -Wall
            -Wextra
            -Wno-compare-reals
            -Wno-conversion
            -Wimplicit-interface
            -fmax-errors=1
            -std=f2008)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(bezier PRIVATE -fcheck=all -fbacktrace
                                          -fimplicit-none -pedantic)
    target_compile_options(bezier PUBLIC -g)
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(bezier PRIVATE -Werror)
    if(${TARGET_NATIVE_ARCH})
      target_compile_options(bezier PRIVATE -march=native)
    endif()
    target_compile_options(bezier PUBLIC -DNDEBUG -O3 -funroll-loops)
  endif()
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  target_compile_options(bezier PUBLIC -warn all)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(bezier PUBLIC -check all)
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(bezier PUBLIC -fast)
  endif()
else()
  message(
    FATAL_ERROR
      "gfortran and ifort are the only supported compilers (current compiler ID is ${CMAKE_Fortran_COMPILER_ID})"
  )
endif()

include(GNUInstallDirs)
install(
  TARGETS bezier
  EXPORT BezierConfig
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT BezierConfig DESTINATION share/bezier/cmake)
export(TARGETS bezier FILE BezierConfig.cmake)

# Restore all Fortran flags.
set(CMAKE_Fortran_FLAGS ${OLD_CMAKE_Fortran_FLAGS})
set(CMAKE_Fortran_FLAGS_DEBUG ${OLD_CMAKE_Fortran_FLAGS_DEBUG})
set(CMAKE_Fortran_FLAGS_RELEASE ${OLD_CMAKE_Fortran_FLAGS_RELEASE})
unset(OLD_CMAKE_Fortran_FLAGS)
unset(OLD_CMAKE_Fortran_FLAGS_DEBUG)
unset(OLD_CMAKE_Fortran_FLAGS_RELEASE)
