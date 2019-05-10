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
!! Revision: 1.0 (10/05/2019)                                                !!
!!           Original functionality implemented in wrapper includes the      !!
!!           creation mechanism used for the sparse matrix for serial        !!
!!           operation as well as the destruction operation to ensure no     !!
!!           memory leaks. setting matrix values as either a single value,   !!
!!           an array of values representing a row, or population via a      !!
!!           dense matrix with row and column vectors indicating their       !!
!!           positional indices. Assembly of the matrix must then be carried !!
!!           out. A feature to view the matrix has also been included. Copy  !!
!!           and duplicate functionality has also been added.                !!
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
    PROCEDURE, PUBLIC , PASS   ::    petscInit => petscInitFM
    PROCEDURE, PUBLIC , PASS   ::    petscCreateMat
    PROCEDURE, PUBLIC , PASS   ::    petscDupeMat
    PROCEDURE, PUBLIC , PASS   ::    petscCopyMat
    GENERIC  , PUBLIC          ::    petscSetMatVals => petscSetMatRowVals, petscSetMatValsFromMat, petscSetMatRowVal
    PROCEDURE, PRIVATE, PASS   ::    petscSetMatRowVals
    PROCEDURE, PRIVATE, PASS   ::    petscSetMatValsFromMat
    PROCEDURE, PRIVATE, PASS   ::    petscSetMatRowVal
    PROCEDURE, PUBLIC , PASS   ::    petscViewMat
    PROCEDURE, PUBLIC , PASS   ::    petscAssembleMat
    PROCEDURE, PUBLIC , PASS   ::    petscDestroyMat
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
  !!              Matrix creation and element alterations                   !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!!

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Initialise PETSc to be used
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscInitFM ( this )
    CLASS(PETScMatrixClass)   ::    this
    PetscErrorCode                  ierr

    call PetscInitialize( PETSC_NULL_CHARACTER, ierr )
    if ( ierr .NE. 0 ) then
      stop "Unable to initialize PETSc"
    endif
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Create an empty PETSc matrix and state the number of rows (nrows) and
  !! columns (ncols) required as well as the number of nonzeros per row (nzprow)
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscCreateMat ( this, nrows, ncols, nzprow )
    CLASS(PETScMatrixClass), INTENT( INOUT )   ::   this
    INTEGER                , INTENT( IN )      ::   ncols
    INTEGER                , INTENT( IN )      ::   nrows
    INTEGER                , INTENT( IN )      ::   nzprow
    PetscErrorCode                  ierr

    ASSOCIATE ( mat => this%mat )
      call MatCreate( PETSC_COMM_WORLD, mat, ierr ); CHKERRA( ierr )
      call MatSetSizes( mat, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols, ierr ); CHKERRA( ierr )
      call MatSetFromOptions( mat, ierr ); CHKERRA(ierr)
      call MatSeqAIJSetPreallocation( mat, nzprow, PETSC_NULL_INTEGER, ierr ); CHKERRA( ierr )
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Allocates memory for a new matrix of the same size as the input vector,
  !! including the nonzero pattern of the original matrix. Nonzero values are
  !! not copied across. The new matrix is created within the duplicate class,
  !! no need to do it externally.
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscDupeMat ( this, mat_in )
    CLASS(PETScMatrixClass), INTENT( INOUT )   ::    this
    TYPE(PETScMatrixClass) , INTENT( IN )      ::    mat_in
    INTEGER                                    ::    n
    PetscErrorCode                                   ierr

    ASSOCIATE ( dupe => this%mat, mat => mat_in%mat )
      call MatDuplicate( mat, MAT_DO_NOT_COPY_VALUES, dupe, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Allocates memory for a new matrix of the same size as the input vector,
  !! including the nonzero pattern of the original matrix. Nonzero values are
  !! copied across. The new matrix is created within the duplicate class,
  !! no need to do it externally.
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscCopyMat ( this, mat_in )
    CLASS(PETScMatrixClass), INTENT( INOUT )   ::    this
    TYPE(PETScMatrixClass) , INTENT( IN )      ::    mat_in
    INTEGER                                    ::    n
    PetscErrorCode                                   ierr

    ASSOCIATE ( copy => this%mat, mat => mat_in%mat )
      call MatDuplicate( mat, MAT_COPY_VALUES, copy, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Calls the insertion function to insert a row of values into the petsc
  !! matrix population procedure
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscSetMatRowVals ( this, irow, icols, vals, c )
    CLASS(PETScMatrixClass)            ::    this
    INTEGER                            ::    irow
    INTEGER                            ::    icols( : )
    REAL(KIND=pres) , INTENT( IN )     ::    vals( : )
    CHARACTER(LEN=*), INTENT( IN )     ::    c
    INTEGER                            ::    trow(1)
    PetscErrorCode                           ierr

    trow = irow

    call insertValsToMatrix( this, trow, icols, vals, c )
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Calls the insertion function to insert a value at a particular row and
  !! column into the petsc matrix data structure
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscSetMatRowVal ( this, irow, icol, val, c )
    CLASS(PETScMatrixClass)            ::    this
    INTEGER                            ::    irow
    INTEGER                            ::    icol
    REAL(KIND=pres) , INTENT( IN )     ::    val
    CHARACTER(LEN=*), INTENT( IN )     ::    c
    PetscErrorCode                           ierr

    call insertValsToMatrix( this, (/irow/), (/icol/), (/val/), c )
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Prepare a dense matrix for insertion into the PETSc matrix class by
  !! creating a row indexed vector out of the input matrix, call the insertion
  !! function.
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscSetMatValsFromMat ( this, irows, icols, vals, c )
    CLASS(PETScMatrixClass)            ::    this
    INTEGER                            ::    irows( : )
    INTEGER                            ::    icols( : )
    REAL(KIND=pres) , INTENT( IN )     ::    vals( :, : )
    CHARACTER(LEN=*), INTENT( IN )     ::    c
    REAL(KIND=pres)                    ::    tvals( size(vals) )
    INTEGER                            ::    i, j, k
    PetscErrorCode                           ierr

    !! populate the temporary vals vector using a row indexed approach that C
    !! implements
    k = 0
    do i = 1, size(vals(:,1))
      do j = 1, size(vals(1,:))
        k = k + 1
        tvals( k ) = vals(i, j)
      enddo
    enddo
    call insertValsToMatrix( this, irows, icols, tvals, c )
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Insert a vector of values (vals) into a particular row (irow) at column
  !! indices specified by icols. the character, c determines whether to insert
  !! the values (replace existing) or add the values (simple addition).
  !! Adjustment is carried out for the zero-indexed nature of C
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE insertValsToMatrix( mat_wrap, irows, icols, vals, c )
    CLASS(PETScMatrixClass)            ::    mat_wrap
    INTEGER                            ::    irows( : )
    INTEGER                            ::    icols( : )
    REAL(KIND=pres) , INTENT( IN )     ::    vals( : )
    CHARACTER(LEN=*), INTENT( IN )     ::    c
    PetscErrorCode                           ierr

    ASSOCIATE ( mat => mat_wrap%mat )
      if ( c == 'insert' ) then
        call MatSetValues( mat, size(irows), irows-1, size(icols), icols-1, vals, INSERT_VALUES, ierr ); CHKERRA(ierr)
      elseif ( c == 'add' ) then
        call MatSetValues( mat, size(irows), irows-1, size(icols), icols-1, vals, ADD_VALUES, ierr  ); CHKERRA(ierr)
      else
        stop 'Error petscSetVecVal(): incorrect chatacter input argument (add/insert)'
      end if
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Assemble the matrix once all values have been set.
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscAssembleMat ( this )
    CLASS(PETScMatrixClass)            ::    this
    PetscErrorCode                           ierr

    ASSOCIATE ( mat => this%mat )
      call MatAssemblyBegin( mat, MAT_FINAL_ASSEMBLY, ierr ); CHKERRA(ierr)
      call MatAssemblyEnd( mat, MAT_FINAL_ASSEMBLY, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! View the PETSc vector
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscViewMat ( this )
    CLASS(PETScMatrixClass), INTENT( INOUT )   ::    this
    PetscErrorCode                  ierr

    ASSOCIATE ( mat => this%mat )
      call MatView( mat, PETSC_VIEWER_STDOUT_SELF, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Destroy the created vector
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  SUBROUTINE petscDestroyMat ( this )
    CLASS(PETScMatrixClass)   ::    this
    PetscErrorCode                  ierr

    ASSOCIATE ( mat => this%mat )
      call MatDestroy( mat, ierr ); CHKERRA(ierr)
    END ASSOCIATE
  END SUBROUTINE


END MODULE
