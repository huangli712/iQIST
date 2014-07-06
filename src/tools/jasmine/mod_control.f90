!=========================================================================!
! project : jasmne
! program : mod_control.f90
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
! type of task
! 1: model calculation
! 2: material calculation
     integer, public, save :: itask

! type of CTQMC algorithm
! 1: general matrices multiplication
! 2: good quantum number: N
! 3: good quantum number: N, Sz, PS
! 4: good quantum number: N, Jz 
     integer, public, save :: ictqmc

! type of crystal field
! 0: no crystal field
! 1: diagonal crystal field
! 2: non-diagonal crystal field 
     integer, public, save :: icf

! type of spin-orbital coupling (SOC) 
! 0: no spin-orbital coupling
! 1: on-site atomic spin-orbital coupling (SOC), H_soc = \lambda * L*S
     integer, public, save :: isoc

! the type of Coulomb interaction U
! 1: Kanamori parameters (Uc, Uv, Jz, Js, Jp)
! 2: Slater-Cordon parameters (F0, F2, F4, F6)
     integer, public, save :: icu

! whether to truncate the Hilbert space for CTQMC's good quantum number algorithm
! at present, we just support truncate the total number of electrons N
! 0: don't truncate it
! 1: truncate it, should set variables min_N and max_N
     integer, public, save :: itrunc

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
! useful when icu = 1
     real(dp), public, save :: Uc  

! interorbital Coulomb interaction
     real(dp), public, save :: Uv  

! Hund's exchange interaction
     real(dp), public, save :: Jz  

! spin-flip interaction
     real(dp), public, save :: Js

! pair-hopping interaction
     real(dp), public, save :: Jp

! Slater-Cordon parameters
! useful when icu = 2
     real(dp), public, save :: F0
     real(dp), public, save :: F2
     real(dp), public, save :: F4
     real(dp), public, save :: F6

! spin-orbit coupling interaction
! 
     real(dp), public, save :: lambda

! the minimal and maximal occupancy number N
! useful when itrunc = 1
     integer, public, save :: min_N
     integer, public, save :: max_N

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

