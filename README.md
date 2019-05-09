# petsc_wrapper
This contains a wrapper for the scientific software library contained within PETSc that enables the user to simply place use statements in their code and call Fortran procedures. Linking in makefile with PETSc is required.

Statements to add to .bashrc

export PETSC_DIR=/home/usr/Documents/Codes/Fortran-Codes/PETSc/linux_opt2
export PETSC_DIR_DEBUG=/home/usr/Documents/Codes/Fortran-Codes/PETSc/linux_optD

Ensure that the linking is done correctly and that the PETSc is compiled with the correct optional flags that you want. An example of one such setup is as follows:

./configure PETSC_ARCH=linux_opt
--prefix=/home/usr/Documents/Codes/Fortran-Codes/PETSc/linux_opt
--with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 --with-debugging=0 --with-x=1
--with-valgrind=0 --with-mpi=1 COPTFLAGS=’-g -O2’ FOPTFLAGS=’-g -O2’.

In .bashrc change PETSC_DIR to be address of where the uncompiled PETSc is, then either alter this makefile to go from there or alternatively change the .bashrc back to the location of the compiled version of PETSc
