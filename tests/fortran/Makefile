CURR_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
SRC_DIR := $(realpath $(CURR_DIR)/../../src/bezier)
QUADPACK_DIR := $(realpath $(SRC_DIR)/quadpack)

FC = gfortran
F77 = .f
F90 = .f90
OBJ = .o
FCFLAGS = -fPIC -O3

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
	$(SRC_DIR)/curve_intersection$(F90)
F77_OBJS := $(patsubst $(QUADPACK_DIR)/%$(F77), $(QUADPACK_DIR)/%$(OBJ), $(F77_SOURCES))
F90_OBJS := $(patsubst $(SRC_DIR)/%$(F90), $(SRC_DIR)/%$(OBJ), $(F90_SOURCES))

test: $(CURR_DIR)/test-bin
	$(CURR_DIR)/test-bin

$(SRC_DIR)/%$(OBJ): $(SRC_DIR)/%$(F90)
	$(FC) $(FCFLAGS) -c $< -o $@

$(QUADPACK_DIR)/%$(OBJ): $(QUADPACK_DIR)/%$(F77)
	$(FC) $(FCFLAGS) -c $< -o $@

$(CURR_DIR)/test-bin: $(F77_OBJS) $(F90_OBJS) $(CURR_DIR)/test.f90
	$(FC) $(FCFLAGS) \
	  -o $(CURR_DIR)/test-bin \
	  $(F77_OBJS) \
	  $(F90_OBJS) \
	  $(CURR_DIR)/test.f90

clean:
	rm -f $(CURR_DIR)/*.mod
	rm -f $(SRC_DIR)/*$(OBJ)
	rm -f $(QUADPACK_DIR)/*$(OBJ)
	rm -f $(CURR_DIR)/test-bin

.PHONY: all test clean
