#! /bin/bash -x
rm PartitionCalculator
gfortran -c Constants.f90
gfortran -c PartitionProcedures.f90 
gfortran -o PartitionCalculator PartitionMain.f90 PartitionProcedures.o Constants.o
rm *.o *.mod
