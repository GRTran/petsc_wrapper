!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Title:   PETSc numerical methods wrapper for Fortran     Date: 10/05/2019 !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Description: This module contains all the necessary include and use       !!
!!              statements for an interface between the PETSc library and    !!
!!              other code that is written that uses it. The following class !!
!!              performs numerical methods on other PETSc vector and matrix  !!
!!              wrappers.                                                    !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Dependencies: The PETSc c++ headers that link to library files, and their !!
!!               Fortran counterparts.                                       !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Revision: 1.0 (10/05/2019)                                                !!
!!           Original functionality implemented in wrapper includes ...
!!
!!
!!
!!
!!
!!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
MODULE PETScNumMethods

#include <petsc/finclude/petscksp.h>
  USE petscksp

  USE PETScMatrix
  USE PETScVector


  IMPLICIT NONE

  TYPE PETScNumMethodsClass
    TYPE(tKSP)          ::    solver
    REAL(KIND=8)        ::    tol
    INTEGER             ::    max_its
    CHARACTER(LEN=10)   ::    tprecon = 'PCCHOLESKY'
  CONTAINS
    PROCEDURE, PUBLIC, PASS   ::    petscInit => petscInitS
    PROCEDURE, PUBLIC, PASS   ::    petscInitKSP
    PROCEDURE, PUBLIC, PASS   ::    petscCreateKSP
    PROCEDURE, PUBLIC, PASS   ::    petscSolveKSP
    PROCEDURE, PUBLIC, PASS   ::    petscDestroyKSP
  END TYPE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Declare any required parameters for the module
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  INTEGER, PRIVATE, PARAMETER    ::    pres = selected_real_kind(8)

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Declare any required subroutine interfaces
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CONTAINS

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Initialise PETSc to be used
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscInitS ( this )
    CLASS(PETScNumMethodsClass)   ::    this
    PetscErrorCode                      ierr

    call PetscInitialize( PETSC_NULL_CHARACTER, ierr )
    if ( ierr .NE. 0 ) then
      stop "Unable to initialize PETSc"
    endif
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Initialise KSP data such as the tolerances, maximum iterations and the
  !! required preconditioner
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscInitKSP ( this )
    CLASS(PETScNumMethodsClass)   ::    this
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Builds the Kyrlov subspace solver method with the required parameters
  !! specified in the KSP initialise subroutine
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscCreateKSP ( this, Ahold )
    CLASS(PETScNumMethodsClass)   ::    this
    TYPE(PETScMatrixClass)        ::    Ahold
    TYPE(tPC)                     ::    precon
    PetscErrorCode                      ierr

    ASSOCIATE( solver => this%solver, A => Ahold%mat )
      !! create the Kyrlov subspace method solver
      call KSPCreate( MPI_COMM_SELF, solver, ierr );  CHKERRA( ierr )
      !! define matrix to be solved and matrix that is basis of preconditioner
      call KSPSetOperators( solver, A, A, ierr );  CHKERRA( ierr )
      !! builds the ksp for a particular solver (PCG in this case)
      call KSPSetType( solver, KSPCG, ierr );  CHKERRA( ierr )
      !! tell the solver that the initial guess is not zero
      call KSPSetInitialGuessNonZero( solver, PETSC_TRUE, ierr ); CHKERRA( ierr )
      !! get the preconditioner
      call KSPGetPC(solver, precon, ierr); CHKERRA( ierr )
      !! set relevant tolerances for the KSP solver
      call KSPSetTolerances( solver, this%tol, this%tol, 10000D0, this%max_its, ierr ); CHKERRA( ierr )
      !! set the preconditioner type
      call PCSetType( precon, TRIM(this%tprecon), ierr ); CHKERRA( ierr )
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Builds the Kyrlov subspace solver method with the required parameters
  !! specified in the KSP initialise subroutine
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscSolveKSP ( this, xhold, bhold, r )
    CLASS(PETScNumMethodsClass)   ::    this
    TYPE(PETScVectorClass)        ::    xhold
    TYPE(PETScVectorClass)        ::    bhold
    REAL(KIND=pres), OPTIONAL     ::    r
    PetscErrorCode                      ierr

    ASSOCIATE( solver => this%solver, x => xhold%vec, b => bhold%vec )
      !! solve Ax=b
      call KSPSolve( solver, b, x, ierr); CHKERRA( ierr )
      !! get the last computed residual norm if requested
      if ( PRESENT( r ) ) then
        call KSPGetResidualNorm( solver, r, ierr ); CHKERRA( ierr )
      endif
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Destroys the Kyrlov subspace solver object
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscdestroyKSP ( this )
    CLASS(PETScNumMethodsClass)   ::    this
    PetscErrorCode                      ierr
    ASSOCIATE( solver => this%solver )
      call KSPDestroy( solver, ierr ); CHKERRA( ierr )
    END ASSOCIATE
  END SUBROUTINE

END MODULE
