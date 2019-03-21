# -*- makefile -*-

MPIEXEC=
PYTHON=python

.PHONY:test
test: run clean

.PHONY:run
run:
	${MPIEXEC} ${PYTHON} petsc_hdf5_save.py

.PHONY:clean
clean:
	${RM} -r __pycache__
	${RM} grid.h5
	
