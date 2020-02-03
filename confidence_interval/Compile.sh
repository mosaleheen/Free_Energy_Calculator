#! /bin/bash -x
rm *.o *.mod *.exe
gfortran -c ConfidenceProcedures.f90
gfortran -o ConfidenceCalculator.exe ConfidenceMain.f90 ConfidenceProcedures.o 
rm *.o *.mod
