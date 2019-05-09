!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Title:   PETSc matrix wrapper for fortran                Date: 09/05/2019 !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Description: This module contains all the necessary include and use       !!
!!              statements for an interface between the PETSc library and    !!
!!              other code that is written that uses it. The following class !!
!!              stores a single PETSc matrix type object, and performs a     !!
!!              range of operations on the allocated sections of the memory  !!
!!              that represent the vector. Information about each procedure  !!
!!              can found above where it is declared.                        !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Dependencies: The PETSc c++ headers that link to library files, and their !!
!!               Fortran counterparts.                                       !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Revision: 1.0 (08/05/2019)                                                !!
!!           Original functionality implemented in wrapper includes ...
!!
!!
!!
!!
!!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
MODULE PETScMatrix

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscmat.h>

  USE petscsys
  USE petscmat

  IMPLICIT NONE

  TYPE PETScMatrixClass
    TYPE( tMat )    ::    mat
  CONTAINS

  END TYPE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Declare any required parameters for the module
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  INTEGER, PARAMETER    ::    pres = selected_real_kind(8)

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Declare any required subroutine interfaces
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



CONTAINS


END MODULE
