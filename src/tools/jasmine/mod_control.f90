!=========================================================================!
! project : clematis
! program : atomic_control.f90
! history : 09/28/2011
! authors : xidai and duliang (email:duleung@gmail.com)
! purpose : some impartant common control variables
! comment : 
!=========================================================================!
  module control
     use constants, only: dp

     implicit none

!-------------------------------------------------------------------------!
!>>> atomic Hamiltonian parameters
!-------------------------------------------------------------------------!
! number of bands
     integer, public, save :: nband

! number of spins
     integer, public, save :: nspin

! number of orbits 
     integer, public, save :: norbs

! number of total electrons
     integer, public, save :: ntots

! number of many-body configurations
     integer, public, save :: ncfgs

! intraorbital Coulomb interaction
     real(dp), public, save :: Uc  

! interorbital Coulomb interaction
     real(dp), public, save :: Uv  

! Hund's exchange interaction
     real(dp), public, save :: Jz  

! spin-flip interaction
     real(dp), public, save :: Js

! pair-hopping interaction
     real(dp), public, save :: Jp

! spin-orbit coupling interaction
     real(dp), public, save :: lamb

! the minimal occupancy number
     integer, public, save :: nmin

! the maximal occupancy number 
     integer, public, save :: nmax

!-------------------------------------------------------------------------!
!>>> MPI related common variables
!-------------------------------------------------------------------------!
! number of processors
     integer, public, save :: nprocs

! the rank of the controller process
     integer, public, save :: master

! the rank of the current process
     integer, public, save :: myrank

  end module control

