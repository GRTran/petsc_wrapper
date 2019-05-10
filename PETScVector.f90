!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Title:   PETSc vector wrapper for fortran                Date: 08/05/2019 !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Description: This module contains all the necessary include and use       !!
!!              statements for an interface between the PETSc library and    !!
!!              other code that is written that uses it. The following class !!
!!              stores a single PETSc vector type, and performs a range of   !!
!!              operations on the allocated sections of the memory that      !!
!!              represent the vector. Information about each procedure can   !!
!!              found above where it is declared.                            !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Dependencies: The PETSc c++ headers that link to library files, and their !!
!!               Fortran counterparts.                                       !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Revision: 1.0 (08/05/2019)                                                !!
!!           Original functionality implemented in wrapper includes the      !!
!!           creation and destruction operations of the vector (both must be !!
!!           used), copy and duplication features. Setting elements of the   !!
!!           vector, printing the vector and performing the dot product      !!
!!           between two vectors.                                            !!
!!           1.01 (10/05/2019)                                               !!
!!           Adjustment for the zero-indexed nature of C++ language that     !!
!!           petsc has been written in.                                      !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
MODULE PETScVector

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>

  USE petscsys
  USE petscvec

  IMPLICIT NONE

  TYPE PETScVectorClass
    TYPE( tVec )    ::    vec
  CONTAINS
    PROCEDURE, PUBLIC , PASS   ::    petscInit => petscInitFV
    PROCEDURE, PUBLIC , PASS   ::    petscCreateVec
    PROCEDURE, PUBLIC , PASS   ::    petscCopyVec
    PROCEDURE, PUBLIC , PASS   ::    petscDupeVec
    PROCEDURE, PUBLIC , PASS   ::    petscViewVec
    GENERIC  , PUBLIC          ::    petscSetVec => petscSetVecVals, petscSetVecVal
    PROCEDURE, PRIVATE, PASS   ::    petscSetVecVals, petscSetVecVal
    PROCEDURE, PUBLIC , PASS   ::    petscSetVecAll
    PROCEDURE, PUBLIC , PASS   ::    petscDestroyVec
    PROCEDURE, PUBLIC , NOPASS ::    petscVecDot
  END TYPE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Declare any required parameters for the module
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  INTEGER, PRIVATE, PARAMETER    ::    pres = selected_real_kind(8)

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Declare any required subroutine interfaces
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



CONTAINS

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!!
  !!              Vector creation and element alterations                   !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!!


  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Initialise PETSc to be used
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscInitFV( this )
    CLASS(PETScVectorClass)   ::    this
    PetscErrorCode                  ierr

    call PetscInitialize( PETSC_NULL_CHARACTER, ierr )
    if ( ierr .NE. 0 ) then
      stop "Unable to initialize PETSc"
    endif

  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Create an empty PETSc vector, of unknown size
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscCreateVec ( this, n )
    CLASS(PETScVectorClass)   ::    this
    INTEGER, INTENT( IN )     ::    n
    PetscErrorCode                  ierr

    ASSOCIATE ( vec => this%vec )
      call VecCreate( PETSC_COMM_WORLD, vec, ierr ); CHKERRA(ierr)
      call VecSetSizes( vec, PETSC_DECIDE, n, ierr ); CHKERRA(ierr)
      call VecSetFromOptions( vec, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Copies a vector and the contents of elements to a new location in memory.
  !! The new vector is created within the copy class, no need to do it externally
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscCopyVec ( this, vec_in )
    CLASS(PETScVectorClass), INTENT( INOUT )   ::    this
    TYPE(PETScVectorClass) , INTENT( IN )      ::    vec_in
    INTEGER                                    ::    n
    PetscErrorCode                                   ierr

    ASSOCIATE ( copy => this%vec, vec => vec_in%vec )
      call VecGetSize( vec, n, ierr ); CHKERRA( ierr )
      call petscCreateVec( this, n )
      call VecCopy( vec, copy, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Allocates memory for a new vector of the same size as the input vector,
  !! values within elements are not copied across. The new vector is created
  !! within the duplicate class, no need to do it externally.
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscDupeVec ( this, vec_in )
    CLASS(PETScVectorClass), INTENT( INOUT )   ::    this
    TYPE(PETScVectorClass) , INTENT( IN )      ::    vec_in
    INTEGER                                    ::    n
    PetscErrorCode                                   ierr

    ASSOCIATE ( dupe => this%vec, vec => vec_in%vec )
      call VecGetSize( vec, n, ierr ); CHKERRA( ierr )
      call petscCreateVec( this, n )
      call VecDuplicate( vec, dupe, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Insert array of values (v) to a vector at indices array locations (i)
  !! or add values (v) to values already in the inidices array locations (i)
  !! adjustment is carried out for the zero-indexed nature of C
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscSetVecVals ( this, i, v, c )
    CLASS(PETScVectorClass)            ::    this
    INTEGER                            ::    i( : )
    REAL(KIND=pres) , INTENT( IN )     ::    v( : )
    CHARACTER(LEN=*), INTENT( IN )     ::    c
    PetscErrorCode                          ierr

    ASSOCIATE ( vec => this%vec )
      if ( c == 'insert' ) then
        call VecSetValues( vec, size(i), i-1, v, INSERT_VALUES, ierr ); CHKERRA(ierr)
      elseif ( c == 'add' ) then
        call VecSetValues( vec, size(i), i-1, v, ADD_VALUES, ierr ); CHKERRA(ierr)
      else
        stop 'Error petscSetVecVal(): incorrect chatacter input argument (add/insert)'
      end if
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Insert value (v) to a vector at index (i)
  !! or add value (v) to value already within cell at index (i)
  !! adjustment is carried out for the zero-indexed nature of C
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscSetVecVal ( this, i, v, c )
    CLASS(PETScVectorClass)            ::    this
    INTEGER                            ::    i
    REAL(KIND=pres) , INTENT( IN )     ::    v
    CHARACTER(LEN=*), INTENT( IN )     ::    c
    PetscErrorCode                          ierr

    ASSOCIATE ( vec => this%vec )
      if ( c == 'insert' ) then
        call VecSetValue( vec, i-1, v, INSERT_VALUES, ierr ); CHKERRA(ierr)
      elseif ( c == 'add' ) then
        call VecSetValue( vec, i-1, v, ADD_VALUES, ierr ); CHKERRA(ierr)
      else
        stop 'Error petscSetVecVal(): incorrect chatacter input argument (add/insert)'
      end if
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Set all elements of a vector to a value (v)
  !! or add value (v) to value already within cell at index (i)
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscSetVecAll ( this, v )
    CLASS(PETScVectorClass)            ::    this
    REAL(KIND=pres) , INTENT( IN )     ::    v
    PetscErrorCode                          ierr

    ASSOCIATE ( vec => this%vec )
        call VecSet( vec, v, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! View the PETSc vector
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscViewVec ( this )
    CLASS(PETScVectorClass)   ::    this
    PetscErrorCode                  ierr

    ASSOCIATE ( vec => this%vec )
      call VecView( vec, PETSC_VIEWER_STDOUT_SELF, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Destroy the created vector
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscDestroyVec ( this )
    CLASS(PETScVectorClass)   ::    this
    PetscErrorCode                  ierr

    ASSOCIATE ( vec => this%vec )
      call VecDestroy( vec, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!!
  !!                      Perform vector operations                         !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!!

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Perform the dot product between two vectors, no pass specified so input
  !! both of the vectors
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  FUNCTION petscVecDot ( x, y ) RESULT ( val_out )
    CLASS(PETScVectorClass)   ::    x, y
    PetscReal                       dot
    REAL(KIND=pres)           ::    val_out
    PetscErrorCode                  ierr

    ASSOCIATE ( x_vec => x%vec , y_vec => y%vec )
      call VecDot( x_vec, y_vec, dot, ierr ); CHKERRA(ierr)
    END ASSOCIATE
    val_out = dot
  END FUNCTION

END MODULE
