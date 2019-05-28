!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Title:   PETSc use all fortran header                    Date: 09/05/2019 !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Description: This module contains all the use statements required to      !!
!!              access classes and methods of the following PETSc wrappers:  !!
!!                - PETScVector                                              !!
!!                - PETScMatrix                                              !!
!!                - PETScNumMethods                                          !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Dependencies: The PETSc wrapper files mentioned above                     !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! Revision: 1.0 (10/05/2019)                                                !!
!!           Original use statements implemented are the following:          !!
!!                - PETScVector                                              !!
!!                - PETScMatrix                                              !!
!!                - PETScNumMethods                                          !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
!! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
MODULE PETScUseAll

  USE PETScVector
  USE PETScMatrix
  USE PETScNumMethods

END MODULE
