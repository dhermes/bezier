CURR_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
SRC_DIR := $(realpath $(CURR_DIR)/../../src/bezier)
QUADPACK_DIR := $(realpath $(SRC_DIR)/quadpack)

FC = gfortran
F77 = .f
F90 = .f90
OBJ = .o
FCFLAGS = -Wall -Wextra -Wno-compare-reals -Wimplicit-interface -Werror -fPIC -fmax-errors=1 -std=f2008 -O0
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
F77_OBJS := $(patsubst $(QUADPACK_DIR)/%$(F77), $(CURR_DIR)/%$(OBJ), $(F77_SOURCES))
F90_OBJS := $(patsubst $(SRC_DIR)/%$(F90), $(CURR_DIR)/%$(OBJ), $(F90_SOURCES))
TEST_OBJS := $(patsubst $(CURR_DIR)/%$(F90), $(CURR_DIR)/%$(OBJ), $(TEST_SOURCES))

test: $(CURR_DIR)/test-bin
	$(CURR_DIR)/test-bin

$(CURR_DIR)/%$(OBJ): $(QUADPACK_DIR)/%$(F77)
	$(FC) $(FC77FLAGS) -c $< -o $@

$(CURR_DIR)/%$(OBJ): $(SRC_DIR)/%$(F90)
	$(FC) $(FCFLAGS) $(COV_FLAGS) -c $< -o $@

$(CURR_DIR)/%$(OBJ): $(CURR_DIR)/%$(F90)
	$(FC) $(FCFLAGS) $(COV_FLAGS) -c $< -o $@

$(CURR_DIR)/test-bin: $(F77_OBJS) $(F90_OBJS) $(TEST_OBJS) $(CURR_DIR)/test.f90
	$(FC) $(FCFLAGS) $(COV_FLAGS) \
	  -o $(CURR_DIR)/test-bin \
	  $(F77_OBJS) \
	  $(F90_OBJS) \
	  $(TEST_OBJS) \
	  $(CURR_DIR)/test.f90

lcov: $(CURR_DIR)/test-bin
	$(CURR_DIR)/test-bin
	lcov --capture \
	  --directory $(CURR_DIR) \
	  --output-file $(CURR_DIR)/coverage.info \
	  $(LCOV_OPTIONS)
	genhtml \
	  $(CURR_DIR)/coverage.info \
	  --output-directory $(CURR_DIR)/lcov-html \
	  $(GENHTML_OPTIONS)

clean:
	rm -f $(CURR_DIR)/*.gcno
	rm -f $(CURR_DIR)/*.gcda
	rm -f $(CURR_DIR)/*.mod
	rm -f $(CURR_DIR)/*$(OBJ)
	rm -f $(CURR_DIR)/test-bin
	rm -f $(CURR_DIR)/coverage.info
	rm -fr $(CURR_DIR)/lcov-html

.PHONY: test lcov clean
