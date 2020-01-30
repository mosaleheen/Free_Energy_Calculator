#! /bin/bash -x
rm Free-energy-calculator
gfortran -c Constants.f90
gfortran -c FEGeneralProcedures.f90 
gfortran -c FESpecificProcedures.f90
gfortran -o Free-energy-calculator FEMain.f90 FESpecificProcedures.o FEGeneralProcedures.o Constants.o
rm *.o *.mod
