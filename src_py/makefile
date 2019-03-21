# -*- makefile -*-

MPIEXEC=mpirun
PYTHON=python

.PHONY:test
test: run clean

.PHONY:run
run:
	${MPIEXEC} -n 2 ${PYTHON} petsc_hdf5_save.py

.PHONY:clean
clean:
	${RM} -r __pycache__
	${RM} grid.h5
	
