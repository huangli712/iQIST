
  module df_control
     use constants, only : dp

     implicit none

!!========================================================================
!!>>> character variables                                              <<<
!!========================================================================

! the code name of the current dual fermion framework
     character(len = 06), public, save :: cname = 'AZALEA'

!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

! number of correlated bands
     integer, public, save :: nband  = 1

! number of spin projection
     integer, public, save :: nspin  = 2

! number of correlated orbitals (= nband * nspin)
     integer, public, save :: norbs  = 2

! number of matsubara frequency for the two-particle green's function
     integer, public, save :: nffrq  = 32

! number of bosonic frequncy for the two-particle green's function
     integer, public, save :: nbfrq  = 8

     integer, public, save :: nkp_x  = 8
     integer, public, save :: nkp_y  = 8
     integer, public, save :: nkp_z  = 8

! number of k-points
     integer, public, save :: nkpts  = 64

! number of dual fermion iteration
     integer, public, save :: ndfit  = 100

! number of BSE iteration
     integer, public, save :: nbsit  = 100

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

! chemical potential or fermi level
!
! note: it should/can be replaced with eimp
     real(dp), public, save :: mune  = 2.00_dp

! inversion of temperature
     real(dp), public, save :: beta  = 8.00_dp

! coupling parameter t for Hubbard model
     real(dp), public, save :: part  = 0.50_dp

! mixing parameter for dynamical mean field theory self-consistent engine
     real(dp), public, save :: bsmix = 0.70_dp
     real(dp), public, save :: dfmix = 0.70_dp

!!========================================================================
!!>>> MPI related common variables                                     <<<
!!========================================================================

! number of processors: default value 1
     integer, public, save :: nprocs = 1

! the id of current process: default value 0
     integer, public, save :: myid   = 0

! denote as the controller process: default value 0
     integer, public, save :: master = 0

! the id of current process in cartesian topology (cid == myid)
     integer, public, save :: cid    = 0

! the x coordinates of current process in cartesian topology
     integer, public, save :: cx     = 0

! the y coordinates of current process in cartesian topology
     integer, public, save :: cy     = 0

  end module df_control
