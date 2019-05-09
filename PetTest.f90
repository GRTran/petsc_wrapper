PROGRAM PetTest
  USE PETScVector


  IMPLICIT NONE

  TYPE( PETScVectorClass )                   ::   vec(4)

  call PETScVectorClassTest( vec )

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
        vec_index(i) = i-1
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
      call vec(1)%petscSetVec ( 0, vec_vals(10), 'add' )  !! add individual value
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
END PROGRAM
