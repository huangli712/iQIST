!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : control module
!!! source  : dmft_control.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/15/2009 by li huang (created)
!!!           01/10/2018 by li huang (last modified)
!!! purpose : define global control parameters for diagrammatic framework
!!!           for dynamical mean field theory.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module control
     use constants, only : dp

     implicit none

!!========================================================================
!!>>> character variables                                              <<<
!!========================================================================

!!
!! @var cname
!!
!! code name of the diagrammatic framework
!!
     character(len = 06), public, save :: cname = 'AZALEA'

!!========================================================================
!!>>> integer variables                                                <<<
!!========================================================================

!!
!! @var nband
!!
!! number of correlated bands
!!
     integer, public, save :: nband  = 1

!!
!! @var nspin
!!
!! number of spin projections
!!
     integer, public, save :: nspin  = 2

!!
!! @var norbs
!!
!! number of correlated orbitals (= nband * nspin)
!!
     integer, public, save :: norbs  = 2

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var nffrq
!!
!! number of fermionic frequencies for the two-particle green's function
!!
     integer, public, save :: nffrq  = 32

!!
!! @var nbfrq
!!
!! number of bosonic frequncies for the two-particle green's function
!!
     integer, public, save :: nbfrq  = 8

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var nkpts
!!
!! number of k-points (totally)
!!
     integer, public, save :: nkpts  = 64

!!
!! @var nkp_x
!!
!! number of k-points (x-axis)
!!
     integer, public, save :: nkp_x  = 8

!!
!! @var nkp_y
!!
!! number of k-points (y-axis)
!!
     integer, public, save :: nkp_y  = 8

!!
!! @var nkp_z
!!
!! number of k-points (z-axis)
!!
     integer, public, save :: nkp_z  = 8

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var ndfit
!!
!! number of dual fermion iteration
!!
     integer, public, save :: ndfit  = 1

!!
!! @var nbsit
!!
!! number of Bethe-Salpter equation iteration
!!
     integer, public, save :: nbsit  = 10

!!========================================================================
!!>>> real variables                                                   <<<
!!========================================================================

!!
!! @var mune
!!
!! chemical potential or fermi level
!!
     real(dp), public, save :: mune  = 2.00_dp

!!
!! @var beta
!!
!! inversion of temperature
!!
     real(dp), public, save :: beta  = 8.00_dp

!!
!! @var part
!!
!! hopping parameter t for Hubbard model
!!
     real(dp), public, save :: part  = 0.50_dp

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!!
!! @var dfmix
!!
!! mixing parameter for dual fermion iteration
!!
     real(dp), public, save :: dfmix = 0.70_dp

!!
!! @var bsmix
!!
!! mixing parameter for Bethe-Salpter equation iteration
!!
     real(dp), public, save :: bsmix = 0.70_dp

!!========================================================================
!!>>> MPI related common variables                                     <<<
!!========================================================================

!!
!! @var nprocs
!!
!! number of processors: default value 1
!!
     integer, public, save :: nprocs = 1

!!
!! @var myid
!!
!! the id of current process: default value 0
!!
     integer, public, save :: myid   = 0

!!
!! @var master
!!
!! denote as the controller process: default value 0
!!
     integer, public, save :: master = 0

!!
!! @var cid
!!
!! the id of current process in cartesian topology (cid == myid)
!!
     integer, public, save :: cid    = 0

!!
!! @var cx
!!
!! the x coordinates of current process in cartesian topology
!!
     integer, public, save :: cx     = 0

!!
!! @var cy
!!
!! the y coordinates of current process in cartesian topology
!!
     integer, public, save :: cy     = 0

  end module control
