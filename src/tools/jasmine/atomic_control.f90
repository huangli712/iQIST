!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : control
!!! source  : atomic_control.f90
!!! type    : modules 
!!! authors : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014
!!! purpose : control parameters
!!! input   :
!!! output  :
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> control parameters
  module control
     use constants, only: dp
  
     implicit none
  
! type of task
! 1: make natural basis inside of the program
! 2: make natural basis outside of the program
     integer, public, save :: itask
  
! type of CTQMC trace algorithm
! 1: general matrices multiplication
! 2: good quantum numbers: N
! 3: good quantum numbers: N, Sz, PS
! 4: good quantum numbers: N, Jz 
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
  
! type of Coulomb interaction U
! 1: Kanamori parameters (Uc, Uv, Jz, Js, Jp)
! 2: Slater-Cordon parameters (F0, F2, F4, F6)
     integer, public, save :: icu
  
! number of bands
     integer, public, save :: nband
  
! number of spins
     integer, public, save :: nspin
  
! number of orbitals
     integer, public, save :: norbs
  
! number of many-body configurations, the dimension of Hilbert space
     integer, public, save :: ncfgs
  
! the following are useful when icu = 1
!--------------------------------------
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
!--------------------------------------
  
! the following are useful when icu = 2
!--------------------------------------
! Coulomb parameters 
     real(dp), public, save :: Ud
  
! Hund's exchange parameters
     real(dp), public, save :: JH
  
! Slater-Cordon parameters
     real(dp), public, save :: F0
     real(dp), public, save :: F2
     real(dp), public, save :: F4
     real(dp), public, save :: F6
!--------------------------------------
     
! spin-orbit coupling interaction
     real(dp), public, save :: lambda
  
! chemical potential
     real(dp), public, save :: mune
  
! MPI related common variables
! number of processors
     integer, public, save :: nprocs
  
! the rank of the controller process
     integer, public, save :: master
  
! the rank of the current process
     integer, public, save :: myrank
  
  end module control
 
