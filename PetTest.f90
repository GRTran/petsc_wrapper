PROGRAM PetTest
  USE PETScUseAll


  IMPLICIT NONE

  TYPE( PETScVectorClass )                   ::   vec(4)
  TYPE( PETScMatrixClass )                   ::   mat(4)

  call PETScVectorClassTest( vec )
  call PETScMatrixClassTest( mat )

CONTAINS

  SUBROUTINE PETScVectorClassTest( vec )
    TYPE( PETScVectorClass )                   ::   vec(:)
    real(kind=8), allocatable, dimension(:)    ::   vec_vals !-- input elements to vector --!
    integer, allocatable, dimension(:)         ::   vec_index !-- input index to vector --!
    integer                                    ::   i !-- counter --!
    integer                                    ::   n !-- number elements in vec --!
    real(kind=8)                               ::   dot_prod

      !! declare the number of elements in the vector
      n = 10

      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !! create the input vector to be used in the test
      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      allocate(vec_vals(n))
      allocate(vec_index(n))
      do i = 1, 10
        vec_vals(i) = 2*i
        vec_index(i) = i
      enddo

      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !! initialise PETSc
      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call vec(1)%petscInit()
      call vec(1)%petscCreateVec(n)
      call vec(1)%petscViewVec()

      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !! set values to the vector and visualise using both the insert and add options
      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call vec(1)%petscSetVec ( vec_index, vec_vals, 'insert' ) !!  insert a new array
      call vec(1)%petscViewVec()
      call vec(1)%petscSetVec ( vec_index, vec_vals, 'add' )  !! add an array of values
      call vec(1)%petscViewVec()
      call vec(1)%petscSetVec ( 1, vec_vals(10), 'add' )  !! add individual value
      call vec(1)%petscViewVec()

      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !! create and set a second vector and perform the dot product between the two
      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call vec(2)%petscCreateVec( n )
      call vec(2)%petscSetVecAll( 2.0D0 ) !!  insert second array
      write(*,*) 'FIRST VECTOR: '
      call vec(1)%petscViewVec()
      write(*,*) 'SECOND VECTOR: '
      call vec(2)%petscViewVec()
      dot_prod = petscVecDot( vec(1), vec(2) )
      write(*,*) 'Dot Product:', dot_prod

      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !! Copy the vector and print, then duplicate the same vector and print
      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call vec(3)%petscCopyVec( vec(1) )
      call vec(3)%petscViewVec()
      call vec(4)%petscDupeVec( vec(1) )
      call vec(4)%petscViewVec()

      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !! destroy the vector
      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      call vec(1)%petscDestroyVec()



      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      !! destroy any other allocatables
      !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      deallocate(vec_vals)
      deallocate(vec_index)
    END SUBROUTINE

    SUBROUTINE PETScMatrixClassTest( mat )
      TYPE( PETScMatrixClass )                   ::   mat(:)
      real(kind=8), allocatable, dimension(:)    ::   vec_vals !-- input elements to vector --!
      real(kind=8), allocatable, dimension(:,:)  ::   dmat_vals !-- input elements as dense matrix --!
      integer, allocatable, dimension(:)         ::   vec_index !-- input index to vector --!
      integer                                    ::   i, j !-- counters --!
      integer                                    ::   n !-- number elements in mat rows and columns --!
      integer                                    ::   nzprow !-- number of non-zero elements in each row --!
      real(kind=8)                               ::   dot_prod

        !! state the number of elements in the vector
        n = 10
        !! state the number of non-zero elements per row
        nzprow = 5

        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !! create the input vector to be used in the test
        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(vec_vals(n))
        allocate(vec_index(n))
        do i = 1, n
          vec_vals(i) = 2.0d0*i
          vec_index(i) = i
        enddo

        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !! create the input matrix to be used in the test
        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        allocate(dmat_vals(n,nzprow))
        do i = 1, nzprow
          do j = 1, n
            dmat_vals(j, i) = j*2.0d0
          enddo
        enddo

        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !! initialise PETSc and create matrix, stating the number of nonzero
        !! elements per row
        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call mat(1)%petscInit()
        call mat(1)%petscCreateMat(n, n, nzprow)


        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !! set the first 5 values in the first row of the index using a vector
        !! of values, then populate next three rows, 5 columns by a dense matrix
        !! then we populate a single value in the matrix
        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call mat(1)%petscSetMatVals( 1, vec_index(1:nzprow), vec_vals(1:nzprow), 'insert' )
        call mat(1)%petscSetMatVals( (/3,4,5/), vec_index(2:nzprow+1), dmat_vals(1:3,:), 'insert' )
        call mat(1)%petscSetMatVals( 7, vec_index(7), vec_vals(10), 'insert' )
        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !! all matrix values must be set before matrix assembly
        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call mat(1)%petscAssembleMat()

        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !! only then can the matrix be viewed
        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        write(*,*)
        write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        write(*,*)
        write(*,*) 'Original: '
        write(*,*)
        call mat(1)%petscViewMat()
        write(*,*)
        write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !! Duplicate and copy the original matrix and print
        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call mat(2)%petscDupeMat( mat(1) )
        call mat(3)%petscCopyMat( mat(1) )
        write(*,*)
        write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        write(*,*)
        write(*,*) 'Duplicate: '
        write(*,*)
        call mat(2)%petscViewMat()
        write(*,*)
        write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        write(*,*)
        write(*,*) 'Copy: '
        write(*,*)
        call mat(3)%petscViewMat()
        write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        write(*,*) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! destroy the vector
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        call mat(1)%petscDestroyMat()
        call mat(2)%petscDestroyMat()
        call mat(3)%petscDestroyMat()



        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !! destroy any other allocatables
        !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        deallocate(vec_vals)
        deallocate(vec_index)
        deallocate(dmat_vals)
      END SUBROUTINE
END PROGRAM
