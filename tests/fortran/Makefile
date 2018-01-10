# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

CURR_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
SRC_DIR := $(realpath $(CURR_DIR)/../../src/bezier)
QUADPACK_DIR := $(realpath $(SRC_DIR)/quadpack)

FC = gfortran
F77 = .f
F90 = .f90
OBJ = .o
FCFLAGS = -Wall -Wextra -Wno-compare-reals -Wno-conversion -Wimplicit-interface -Werror -fPIC -fmax-errors=1 -std=f2008 -O0 -g -fcheck=all -fbacktrace -fimplicit-none -pedantic
FC77FLAGS = -fPIC -O0
COV_FLAGS = -fprofile-arcs -ftest-coverage

ifdef BRANCH_COVERAGE
LCOV_OPTIONS = --rc lcov_branch_coverage=1
GENHTML_OPTIONS = --branch-coverage
else
LCOV_OPTIONS =
GENHTML_OPTIONS =
endif

# NOTE: **Must** specify the order for source files.
F77_SOURCES := \
	$(QUADPACK_DIR)/d1mach$(F77) \
	$(QUADPACK_DIR)/dqelg$(F77) \
	$(QUADPACK_DIR)/dqpsrt$(F77) \
	$(QUADPACK_DIR)/dqk21$(F77) \
	$(QUADPACK_DIR)/dqagse$(F77)
F90_SOURCES := \
	$(SRC_DIR)/types$(F90) \
	$(SRC_DIR)/status$(F90) \
	$(SRC_DIR)/helpers$(F90) \
	$(SRC_DIR)/curve$(F90) \
	$(SRC_DIR)/surface$(F90) \
	$(SRC_DIR)/curve_intersection$(F90) \
	$(SRC_DIR)/surface_intersection$(F90)
TEST_SOURCES := \
	$(CURR_DIR)/unit_test_helpers$(F90) \
	$(CURR_DIR)/test_helpers$(F90) \
	$(CURR_DIR)/test_curve$(F90) \
	$(CURR_DIR)/test_surface$(F90) \
	$(CURR_DIR)/test_curve_intersection$(F90) \
	$(CURR_DIR)/test_surface_intersection$(F90)
FNL_TEST_SOURCES := \
	$(CURR_DIR)/functional_test_helpers$(F90) \
	$(CURR_DIR)/functional_curve$(F90)
F77_OBJS := $(patsubst $(QUADPACK_DIR)/%$(F77), $(CURR_DIR)/%$(OBJ), $(F77_SOURCES))
F90_OBJS := $(patsubst $(SRC_DIR)/%$(F90), $(CURR_DIR)/%$(OBJ), $(F90_SOURCES))
TEST_OBJS := $(patsubst $(CURR_DIR)/%$(F90), $(CURR_DIR)/%$(OBJ), $(TEST_SOURCES))
FNL_TEST_OBJS := $(patsubst $(CURR_DIR)/%$(F90), $(CURR_DIR)/%$(OBJ), $(FNL_TEST_SOURCES))

help:
	@echo 'Makefile for Fortran unit and functional tests'
	@echo ''
	@echo 'Usage:'
	@echo '   make unit          Build and run unit tests'
	@echo '   make lcov          Run unit tests with coverage'
	@echo '   make functional    Build and run functional tests'
	@echo '   make clean         Clean generated files'
	@echo ''

$(CURR_DIR)/%$(OBJ): $(QUADPACK_DIR)/%$(F77)
	$(FC) $(FC77FLAGS) -c $< -o $@

$(CURR_DIR)/%$(OBJ): $(SRC_DIR)/%$(F90)
	$(FC) $(FCFLAGS) $(COV_FLAGS) -c $< -o $@

$(CURR_DIR)/%$(OBJ): $(CURR_DIR)/%$(F90)
	$(FC) $(FCFLAGS) $(COV_FLAGS) -c $< -o $@

$(CURR_DIR)/unit-test-bin: $(F77_OBJS) $(F90_OBJS) $(TEST_OBJS) $(CURR_DIR)/unit_test$(F90)
	$(FC) $(FCFLAGS) $(COV_FLAGS) \
	  -o $(CURR_DIR)/unit-test-bin \
	  $(F77_OBJS) \
	  $(F90_OBJS) \
	  $(TEST_OBJS) \
	  $(CURR_DIR)/unit_test$(F90)

unit: $(CURR_DIR)/unit-test-bin
	$(CURR_DIR)/unit-test-bin

lcov: $(CURR_DIR)/unit-test-bin
	$(CURR_DIR)/unit-test-bin
	lcov --capture \
	  --directory $(CURR_DIR) \
	  --output-file $(CURR_DIR)/coverage.info \
	  $(LCOV_OPTIONS)
	genhtml \
	  $(CURR_DIR)/coverage.info \
	  --output-directory $(CURR_DIR)/lcov-html \
	  $(GENHTML_OPTIONS)

$(CURR_DIR)/functional-test-bin: $(F77_OBJS) $(F90_OBJS) $(FNL_TEST_OBJS) $(CURR_DIR)/functional$(F90)
	$(FC) $(FCFLAGS) $(COV_FLAGS) \
	  -o $(CURR_DIR)/functional-test-bin \
	  $(F77_OBJS) \
	  $(F90_OBJS) \
	  $(FNL_TEST_OBJS) \
	  $(CURR_DIR)/functional$(F90)

functional: $(CURR_DIR)/functional-test-bin
	$(CURR_DIR)/functional-test-bin

clean:
	rm -f $(CURR_DIR)/*.gcno
	rm -f $(CURR_DIR)/*.gcda
	rm -f $(CURR_DIR)/*.mod
	rm -f $(CURR_DIR)/*$(OBJ)
	rm -f $(CURR_DIR)/functional-test-bin
	rm -f $(CURR_DIR)/unit-test-bin
	rm -f $(CURR_DIR)/coverage.info
	rm -fr $(CURR_DIR)/unit-test-bin.dSYM
	rm -fr $(CURR_DIR)/functional-test-bin.dSYM
	rm -fr $(CURR_DIR)/lcov-html

.PHONY: help unit lcov functional clean
