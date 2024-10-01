cmake_minimum_required(VERSION 3.17)

project(
  pybezier
  LANGUAGES C
  VERSION ${SKBUILD_PROJECT_VERSION})

find_package(Python 3.10 REQUIRED COMPONENTS Interpreter Development.Module
                                             NumPy)

add_subdirectory(src/fortran)

python_add_library(_speedup MODULE src/python/bezier/_speedup.c)
target_link_libraries(_speedup PRIVATE bezier Python::NumPy)
target_compile_definitions(
  _speedup
  PRIVATE NPY_NO_DEPRECATED_API
  PRIVATE NPY_1_7_API_VERSION)

if(APPLE)
  set_target_properties(_speedup PROPERTIES INSTALL_RPATH @loader_path/lib)
else()
  set_target_properties(_speedup PROPERTIES INSTALL_RPATH $ORIGIN/../lib)
endif()

install(TARGETS _speedup DESTINATION .)
