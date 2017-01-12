all: bezier/_speedup.so

bezier/_speedup.so: bezier/speedup.f90 .f2py_f2cmap
	f2py -c -m _speedup bezier/speedup.f90
	mv _speedup.so bezier

clean:
	rm -f bezier/_speedup.so

.PHONY: all clean
