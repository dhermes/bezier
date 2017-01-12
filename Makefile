all: src/bezier/_speedup.so

src/bezier/_speedup.so: src/bezier/speedup.f90 .f2py_f2cmap
	f2py -c -m _speedup src/bezier/speedup.f90
	mv _speedup.so src/bezier

clean:
	rm -f src/bezier/_speedup.so

.PHONY: all clean
