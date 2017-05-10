!!!-----------------------------------------------------------------------
!!! project : manjushaka
!!! program : ctqmc_core module
!!!           ctqmc_clur module
!!!           ctqmc_flvr module
!!!           ctqmc_mesh module
!!!           ctqmc_meat module
!!!           ctqmc_umat module
!!!           ctqmc_mmat module
!!!           ctqmc_gmat module
!!!           ctqmc_wmat module
!!!           ctqmc_smat module
!!!           context    module
!!! source  : ctqmc_context.f90
!!! type    : modules
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           05/10/2017 by li huang (last modified)
!!! purpose : define the key data structure and global arrays/variables
!!!           for hybridization expansion version continuous time quantum
!!!           Monte Carlo (CTQMC) quantum impurity solver and dynamical
!!!           mean field theory (DMFT) self-consistent engine
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module ctqmc_core                                                <<<
!!========================================================================

!!
!! @mod ctqmc_core
!!
!! containing core (internal) variables used by continuous time quantum
!! Monte Carlo quantum impurity solver
!!
  module ctqmc_core
     use constants, only : dp, zero

     implicit none

!!
!! @var ckink
!!
!! current perturbation expansion order
!!
     integer, public, save  :: ckink = 0

! sign change related with current diagram update operation
     integer, public, save  :: csign = 1

! counter for negative sign, used to measure the sign problem
     integer, public, save  :: cnegs = 0

! averaged sign values, used to measure the sign problem
     integer, public, save  :: caves = 0

! current status of spin-orbital coupling
! if cssoc = 0, no spin-orbital coupling,
! if cssoc = 1, atomic spin-orbital coupling
! note: this variable is determined by atom.cix, do not setup it manually
     integer, public, save  :: cssoc = 0

!-------------------------------------------------------------------------
!::: core variables: real, matrix trace                                :::
!-------------------------------------------------------------------------

! matrix trace of flavor part, current value
     real(dp), public, save :: matrix_ptrace = zero

! matrix trace of flavor part, proposed value
     real(dp), public, save :: matrix_ntrace = zero

!-------------------------------------------------------------------------
!::: core variables: real, insert action counter                       :::
!-------------------------------------------------------------------------

!!
!! @var ins_t
!!
!! insert kink (operators pair) statistics: total insert count
!!
     real(dp), public, save :: ins_t = zero

!!
!! @var ins_a
!!
!! insert kink (operators pair) statistics: total accepted insert count
!!
     real(dp), public, save :: ins_a = zero

!!
!! @var ins_r
!!
!! insert kink (operators pair) statistics: total rejected insert count
!!
     real(dp), public, save :: ins_r = zero

!-------------------------------------------------------------------------
!::: core variables: real, remove action counter                       :::
!-------------------------------------------------------------------------

!!
!! @var rmv_t
!!
!! remove kink (operators pair) statistics: total remove count
!!
     real(dp), public, save :: rmv_t = zero

!!
!! @var rmv_a
!!
!! remove kink (operators pair) statistics: total accepted remove count
!!
     real(dp), public, save :: rmv_a = zero

!!
!! @var rmv_r
!!
!! remove kink (operators pair) statistics: total rejected remove count
!!
     real(dp), public, save :: rmv_r = zero

!-------------------------------------------------------------------------
!::: core variables: real, lshift action counter                       :::
!-------------------------------------------------------------------------

!!
!! @var lsh_t
!!
!! lshift kink (operators pair) statistics: total lshift count
!!
     real(dp), public, save :: lsh_t = zero

!!
!! @var lsh_a
!!
!! lshift kink (operators pair) statistics: total accepted lshift count
!!
     real(dp), public, save :: lsh_a = zero

!!
!! @var lsh_r
!!
!! lshift kink (operators pair) statistics: total rejected lshift count
!!
     real(dp), public, save :: lsh_r = zero

!-------------------------------------------------------------------------
!::: core variables: real, rshift action counter                       :::
!-------------------------------------------------------------------------

!!
!! @var rsh_t
!!
!! rshift kink (operators pair) statistics: total rshift count
!!
     real(dp), public, save :: rsh_t = zero

!!
!! @var rsh_a
!!
!! rshift kink (operators pair) statistics: total accepted rshift count
!!
     real(dp), public, save :: rsh_a = zero

!!
!! @var rsh_r
!!
!! rshift kink (operators pair) statistics: total rejected rshift count
!!
     real(dp), public, save :: rsh_r = zero

!-------------------------------------------------------------------------
!::: core variables: real, reflip action counter                       :::
!-------------------------------------------------------------------------

! reflip kink (operators pair) statistics: total reflip count
     real(dp), public, save :: reflip_t = zero

! reflip kink (operators pair) statistics: total accepted reflip count
     real(dp), public, save :: reflip_a = zero

! reflip kink (operators pair) statistics: total rejected reflip count
     real(dp), public, save :: reflip_r = zero

  end module ctqmc_core

!!========================================================================
!!>>> module ctqmc_clur                                                <<<
!!========================================================================

!!
!! @mod ctqmc_clur
!!
!! containing perturbation expansion series related arrays (colour part)
!! used by continuous time quantum Monte Carlo quantum impurity solver
!!
  module ctqmc_clur
     use constants, only : dp

     use stack, only : istack
     use stack, only : istack_create
     use stack, only : istack_destroy

     implicit none

!!
!! @var index_s
!!
!! memory address index for the imaginary time \tau_s
!!
     integer, public, save, allocatable  :: index_s(:,:)

!!
!! @var index_e
!!
!! memory address index for the imaginary time \tau_e
!!
     integer, public, save, allocatable  :: index_e(:,:)

!!
!! @var time_s
!!
!! imaginary time \tau_s of creation operators
!!
     real(dp), public, save, allocatable :: time_s(:,:)

!!
!! @var time_e
!!
!! imaginary time \tau_e of annihilation operators
!!
     real(dp), public, save, allocatable :: time_e(:,:)

!!
!! @var exp_s
!!
!! exp(i\omega t), s means creation operators
!!
     complex(dp), public, save, allocatable :: exp_s(:,:,:)

!!
!! @var exp_e
!!
!! exp(i\omega t), e means annihilation operators
!!
     complex(dp), public, save, allocatable :: exp_e(:,:,:)

!!
!! @var empty_s
!!
!! container for the empty (unoccupied) memory address index
!!
     type (istack), public, save, allocatable :: empty_s(:)

!!
!! @var empty_e
!!
!! container for the empty (unoccupied) memory address index
!!
     type (istack), public, save, allocatable :: empty_e(:)

  end module ctqmc_clur

!!========================================================================
!!>>> module ctqmc_flvr                                                <<<
!!========================================================================

!!>>> containing perturbation expansion series related arrays (flavor part)
!!>>> used by continuous time quantum Monte Carlo quantum impurity solver
  module ctqmc_flvr
     use constants, only : dp
     use stack, only : istack, istack_create, istack_destroy

     implicit none

! container for the empty (unoccupied) memory address index of operators
     type (istack), public, save :: empty_v

! memory address index for the imaginary time \tau (auxiliary)
     integer, public, save, allocatable  :: index_t(:)

! memory address index for the imaginary time \tau
     integer, public, save, allocatable  :: index_v(:)

! to record type of operators, 1 means create operators, 0 means destroy operators
     integer, public, save, allocatable  :: type_v(:)

! to record flavor of operators, from 1 to norbs
     integer, public, save, allocatable  :: flvr_v(:)

! imaginary time \tau for create and destroy operators
     real(dp), public, save, allocatable :: time_v(:)

! exp(-H\tau), exponent matrix for local hamiltonian multiply \tau (the last point)
     real(dp), public, save, allocatable :: expt_t(:,:)

! exp(-H\tau), exponent matrix for local hamiltonian multiply \tau
     real(dp), public, save, allocatable :: expt_v(:,:)

  end module ctqmc_flvr

!!========================================================================
!!>>> module ctqmc_mesh                                                <<<
!!========================================================================

!!
!! @mod ctqmc_mesh
!!
!! containing mesh related arrays used by continuous time quantum Monte
!! Carlo quantum impurity solver
!!
  module ctqmc_mesh
     use constants, only : dp

     implicit none

!!
!! @var tmesh
!!
!! imaginary time mesh
!!
     real(dp), public, save, allocatable :: tmesh(:)

!!
!! @var rmesh
!!
!! real matsubara frequency mesh
!!
     real(dp), public, save, allocatable :: rmesh(:)

!!
!! @var lmesh
!!
!! uniform mesh on interval [-1,1] for legendre orthogonal polynomial
!!
     real(dp), public, save, allocatable :: lmesh(:)

!!
!! @var rep_l
!!
!! legendre orthogonal polynomial defined on [-1,1]
!!
     real(dp), public, save, allocatable :: rep_l(:,:)

  end module ctqmc_mesh

!!========================================================================
!!>>> module ctqmc_meat                                                <<<
!!========================================================================

!!
!! @mod ctqmc_meat
!!
!! containing physical observables related arrays used by continuous time
!! quantum Monte Carlo quantum impurity solver
!!
  module ctqmc_meat ! to tell you a truth, meat means MEAsuremenT
     use constants, only : dp

     implicit none

!!
!! @var hist
!!
!! histogram for perturbation expansion series
!!
     real(dp), public, save, allocatable :: hist(:)

! auxiliary physical observables
! paux(01) : total energy, Etot
! paux(02) : potential engrgy, Epot
! paux(03) : kinetic energy, Ekin
! paux(04) : magnetic moment, < Sz >
! paux(05) : average of occupation, < N > = < N^1 > = < N1 >
! paux(06) : average of occupation square, < N^2 > = < N2 >
! paux(07) : high order of K, < K^2 > = < K2 >
! paux(08) : high order of K, < K^3 > = < K3 >
! paux(09) : high order of K, < K^4 > = < K4 >
!
! note: K = current perturbation expansion order X 2. The < K2 >, < K3 >,
! and < K4 > can be used to calculate the skewness and kurtosis of the
! perturbation expansion order. Of course, < K1 > is essential. It can be
! calculated from Ekin.
     real(dp), public, save, allocatable :: paux(:)

! probability of eigenstates of local hamiltonian matrix
     real(dp), public, save, allocatable :: prob(:)

! impurity occupation number, < n_i >
     real(dp), public, save, allocatable :: nmat(:)

! impurity double occupation number matrix, < n_i n_j >
     real(dp), public, save, allocatable :: nnmat(:,:)

! number of operators, < k >
     real(dp), public, save, allocatable :: kmat(:)

! square of number of operators, < k^2 >
     real(dp), public, save, allocatable :: kkmat(:,:)

! number of operators at left half axis, < k_l >
     real(dp), public, save, allocatable :: lmat(:)

! number of operators at right half axis, < k_r >
     real(dp), public, save, allocatable :: rmat(:)

! used to evaluate fidelity susceptibility, < k_l k_r >
     real(dp), public, save, allocatable :: lrmat(:,:)

! used to calculate two-particle green's function, real part
     real(dp), public, save, allocatable :: g2_re(:,:,:,:,:)

! used to calculate two-particle green's function, imaginary part
     real(dp), public, save, allocatable :: g2_im(:,:,:,:,:)

! particle-particle pair susceptibility, real part
     real(dp), public, save, allocatable :: ps_re(:,:,:,:,:)

! particle-particle pair susceptibility, imaginary part
     real(dp), public, save, allocatable :: ps_im(:,:,:,:,:)

  end module ctqmc_meat

!!========================================================================
!!>>> module ctqmc_umat                                                <<<
!!========================================================================

!!
!! @mod ctqmc_umat
!!
!! containing auxiliary arrays used by continuous time quantum Monte
!! Carlo quantum impurity solver
!!
  module ctqmc_umat
     use constants, only : dp

     implicit none

!-------------------------------------------------------------------------
!::: ctqmc status variables                                            :::
!-------------------------------------------------------------------------

! current perturbation expansion order for different flavor channel
     integer,  public, save, allocatable :: rank(:)

! diagonal elements of current matrix product of flavor part
! it is used to calculate the probability of eigenstates
     real(dp), public, save, allocatable :: diag(:,:)

!-------------------------------------------------------------------------
!::: input data variables                                              :::
!-------------------------------------------------------------------------

! symmetry properties for correlated orbitals
     integer,  public, save, allocatable :: symm(:)

! impurity level for correlated orbitals
     real(dp), public, save, allocatable :: eimp(:)

! eigenvalues for local hamiltonian matrix
     real(dp), public, save, allocatable :: eigs(:)

! occupation number for the eigenstates of local hamiltonian matrix
     real(dp), public, save, allocatable :: naux(:)

! total spin for the eigenstates of local hamiltonian matrix
     real(dp), public, save, allocatable :: saux(:)

  end module ctqmc_umat

!!========================================================================
!!>>> module ctqmc_mmat                                                <<<
!!========================================================================

!!
!! @mod ctqmc_mmat
!!
!! containing M-matrix and G-matrix related arrays used by continuous
!! time quantum Monte Carlo quantum impurity solver
!!
  module ctqmc_mmat
     use constants, only : dp

     implicit none

!!
!! @var lspace
!!
!! helper matrix for evaluating M & G matrices
!!
     real(dp), public, save, allocatable :: lspace(:,:)

!!
!! @var rspace
!!
!! helper matrix for evaluating M & G matrices
!!
     real(dp), public, save, allocatable :: rspace(:,:)

!!
!! @var mmat
!!
!! M matrix, $ \mathscr{M} $
!!
     real(dp), public, save, allocatable :: mmat(:,:,:)

!!
!! @var lsaves
!!
!! helper matrix for evaluating G matrix
!!
     complex(dp), public, save, allocatable :: lsaves(:,:)

!!
!! @var rsaves
!!
!! helper matrix for evaluating G matrix
!!
     complex(dp), public, save, allocatable :: rsaves(:,:)

!!
!! @var gmat
!!
!! G matrix, $ \mathscr{G} $
!!
     complex(dp), public, save, allocatable :: gmat(:,:,:)

  end module ctqmc_mmat

!!========================================================================
!!>>> module ctqmc_gmat                                                <<<
!!========================================================================

!!
!! @mod ctqmc_gmat
!!
!! containing green's function matrix related arrays used by continuous
!! time quantum Monte Carlo quantum impurity solver
!!
  module ctqmc_gmat
     use constants, only : dp

     implicit none

!!
!! @var gtau
!!
!! impurity green's function in imaginary time axis
!!
     real(dp), public, save, allocatable :: gtau(:,:,:)

!!
!! @var grnf
!!
!! impurity green's function in matsubara frequency axis
!!
     complex(dp), public, save, allocatable :: grnf(:,:,:)

  end module ctqmc_gmat

!!========================================================================
!!>>> module ctqmc_wmat                                                <<<
!!========================================================================

!!
!! @mod ctqmc_wmat
!!
!! containing weiss's function and hybridization function related arrays
!! used by continuous time quantum Monte Carlo quantum impurity solver
!!
  module ctqmc_wmat
     use constants, only : dp

     implicit none

!!
!! @var wtau
!!
!! bath weiss's function in imaginary time axis
!!
     real(dp), public, save, allocatable :: wtau(:,:,:)

!!
!! @var htau
!!
!! hybridization function in imaginary time axis
!!
     real(dp), public, save, allocatable :: htau(:,:,:)

!!
!! @var hsed
!!
!! second order derivates for hybridization function
!!
     real(dp), public, save, allocatable :: hsed(:,:,:)

!!
!! @var wssf
!!
!! bath weiss's function in matsubara frequency axis
!!
     complex(dp), public, save, allocatable :: wssf(:,:,:)

!!
!! @var hybf
!!
!! hybridization function in matsubara frequency axis
!!
     complex(dp), public, save, allocatable :: hybf(:,:,:)

  end module ctqmc_wmat

!!========================================================================
!!>>> module ctqmc_smat                                                <<<
!!========================================================================

!!
!! @mod ctqmc_smat
!!
!! containing self-energy function matrix related arrays used by
!! continuous time quantum Monte Carlo quantum impurity solver
!!
  module ctqmc_smat
     use constants, only : dp

     implicit none

!!
!! @var sig1
!!
!! self-energy function in matsubara frequency axis
!!
     complex(dp), public, save, allocatable :: sig1(:,:,:)

!!
!! @var sig2
!!
!! self-energy function in matsubara frequency axis
!!
     complex(dp), public, save, allocatable :: sig2(:,:,:)

  end module ctqmc_smat

!!========================================================================
!!>>> module context                                                   <<<
!!========================================================================

!!>>> containing memory management subroutines and define global variables
  module context
     use constants
     use control

     use ctqmc_core
     use ctqmc_clur
     use ctqmc_flvr
     use ctqmc_mesh
     use ctqmc_meat
     use ctqmc_umat
     use ctqmc_mmat
     use ctqmc_gmat
     use ctqmc_wmat
     use ctqmc_smat

     implicit none

!!========================================================================
!!>>> declare global variables                                         <<<
!!========================================================================

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

! declaration of module procedures: allocate memory
     public :: ctqmc_allocate_memory_clur
     public :: ctqmc_allocate_memory_flvr
     public :: ctqmc_allocate_memory_mesh
     public :: ctqmc_allocate_memory_meat
     public :: ctqmc_allocate_memory_umat
     public :: ctqmc_allocate_memory_mmat
     public :: ctqmc_allocate_memory_gmat
     public :: ctqmc_allocate_memory_wmat
     public :: ctqmc_allocate_memory_smat

! declaration of module procedures: deallocate memory
     public :: ctqmc_deallocate_memory_clur
     public :: ctqmc_deallocate_memory_flvr
     public :: ctqmc_deallocate_memory_mesh
     public :: ctqmc_deallocate_memory_meat
     public :: ctqmc_deallocate_memory_umat
     public :: ctqmc_deallocate_memory_mmat
     public :: ctqmc_deallocate_memory_gmat
     public :: ctqmc_deallocate_memory_wmat
     public :: ctqmc_deallocate_memory_smat

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!>>> ctqmc_allocate_memory_clur: allocate memory for clur-related variables
  subroutine ctqmc_allocate_memory_clur()
     implicit none

! local variables
! loop index
     integer :: i

! allocate memory
     allocate(index_s(mkink,norbs),     stat=istat)
     allocate(index_e(mkink,norbs),     stat=istat)

     allocate(time_s(mkink,norbs),      stat=istat)
     allocate(time_e(mkink,norbs),      stat=istat)

     allocate(exp_s(nfreq,mkink,norbs), stat=istat)
     allocate(exp_e(nfreq,mkink,norbs), stat=istat)

     allocate(empty_s(norbs),           stat=istat)
     allocate(empty_e(norbs),           stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_clur','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     index_s = 0
     index_e = 0

     time_s  = zero
     time_e  = zero

     exp_s   = czero
     exp_e   = czero

     do i=1,norbs
         call istack_create(empty_s(i), mkink)
         call istack_create(empty_e(i), mkink)
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_allocate_memory_clur

!!>>> ctqmc_allocate_memory_flvr: allocate memory for flvr-related variables
  subroutine ctqmc_allocate_memory_flvr()
     implicit none

! allocate memory
     allocate(index_t(mkink),      stat=istat)
     allocate(index_v(mkink),      stat=istat)

     allocate(type_v(mkink),       stat=istat)
     allocate(flvr_v(mkink),       stat=istat)

     allocate(time_v(mkink),       stat=istat)

     allocate(expt_t(ncfgs,  4  ), stat=istat)
     allocate(expt_v(ncfgs,mkink), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_flvr','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     index_t = 0
     index_v = 0

     type_v  = 1
     flvr_v  = 1

     time_v  = zero

     expt_t  = zero
     expt_v  = zero

     call istack_create(empty_v, mkink)

     return
  end subroutine ctqmc_allocate_memory_flvr

!!>>> ctqmc_allocate_memory_mesh: allocate memory for mesh-related variables
  subroutine ctqmc_allocate_memory_mesh()
     implicit none

! allocate memory
     allocate(tmesh(ntime),       stat=istat)
     allocate(rmesh(mfreq),       stat=istat)

     allocate(pmesh(legrd),       stat=istat)
     allocate(qmesh(chgrd),       stat=istat)

     allocate(ppleg(legrd,lemax), stat=istat)
     allocate(qqche(chgrd,chmax), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_mesh','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     tmesh = zero
     rmesh = zero

     pmesh = zero
     qmesh = zero

     ppleg = zero
     qqche = zero

     return
  end subroutine ctqmc_allocate_memory_mesh

!!>>> ctqmc_allocate_memory_meat: allocate memory for meat-related variables
  subroutine ctqmc_allocate_memory_meat()
     implicit none

! allocate memory
     allocate(hist(mkink),        stat=istat)

     allocate(paux(  9  ),        stat=istat)
     allocate(prob(ncfgs),        stat=istat)

     allocate(nmat(norbs),        stat=istat)
     allocate(nnmat(norbs,norbs), stat=istat)
     allocate(kmat(norbs),        stat=istat)
     allocate(kkmat(norbs,norbs), stat=istat)
     allocate(lmat(norbs),        stat=istat)
     allocate(rmat(norbs),        stat=istat)
     allocate(lrmat(norbs,norbs), stat=istat)

     allocate(g2_re(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(g2_im(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(ps_re(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(ps_im(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_meat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     hist  = zero

     paux  = zero
     prob  = zero

     nmat  = zero
     nnmat = zero
     kmat  = zero
     kkmat = zero
     lmat  = zero
     rmat  = zero
     lrmat = zero

     g2_re = zero
     g2_im = zero
     ps_re = zero
     ps_im = zero

     return
  end subroutine ctqmc_allocate_memory_meat

!!>>> ctqmc_allocate_memory_umat: allocate memory for umat-related variables
  subroutine ctqmc_allocate_memory_umat()
     implicit none

! allocate memory
     allocate(rank(norbs),        stat=istat)

     allocate(diag(ncfgs,  2  ),  stat=istat)

     allocate(symm(norbs),        stat=istat)

     allocate(eimp(norbs),        stat=istat)
     allocate(eigs(ncfgs),        stat=istat)
     allocate(naux(ncfgs),        stat=istat)
     allocate(saux(ncfgs),        stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_umat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     rank  = 0

     diag  = zero

     symm  = 0

     eimp  = zero
     eigs  = zero
     naux  = zero
     saux  = zero

     return
  end subroutine ctqmc_allocate_memory_umat

!!>>> ctqmc_allocate_memory_mmat: allocate memory for mmat-related variables
  subroutine ctqmc_allocate_memory_mmat()
     implicit none

! allocate memory
     allocate(lspace(mkink,norbs),     stat=istat)
     allocate(rspace(mkink,norbs),     stat=istat)

     allocate(mmat(mkink,mkink,norbs), stat=istat)

     allocate(lsaves(nfreq,norbs),     stat=istat)
     allocate(rsaves(nfreq,norbs),     stat=istat)

     allocate(gmat(nfreq,norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_mmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     lspace = zero
     rspace = zero

     mmat   = zero

     lsaves = czero
     rsaves = czero

     gmat   = czero

     return
  end subroutine ctqmc_allocate_memory_mmat

!!>>> ctqmc_allocate_memory_gmat: allocate memory for gmat-related variables
  subroutine ctqmc_allocate_memory_gmat()
     implicit none

! allocate memory
     allocate(gtau(ntime,norbs,norbs), stat=istat)

     allocate(grnf(mfreq,norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_gmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     gtau = zero

     grnf = czero

     return
  end subroutine ctqmc_allocate_memory_gmat

!!>>> ctqmc_allocate_memory_wmat: allocate memory for wmat-related variables
  subroutine ctqmc_allocate_memory_wmat()
     implicit none

! allocate memory
     allocate(wtau(ntime,norbs,norbs), stat=istat)
     allocate(htau(ntime,norbs,norbs), stat=istat)
     allocate(hsed(ntime,norbs,norbs), stat=istat)

     allocate(wssf(mfreq,norbs,norbs), stat=istat)
     allocate(hybf(mfreq,norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_wmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     wtau = zero
     htau = zero
     hsed = zero

     wssf = czero
     hybf = czero

     return
  end subroutine ctqmc_allocate_memory_wmat

!!>>> ctqmc_allocate_memory_smat: allocate memory for smat-related variables
  subroutine ctqmc_allocate_memory_smat()
     implicit none

! allocate memory
     allocate(sig1(mfreq,norbs,norbs), stat=istat)
     allocate(sig2(mfreq,norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_smat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     sig1 = czero
     sig2 = czero

     return
  end subroutine ctqmc_allocate_memory_smat

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!>>> ctqmc_deallocate_memory_clur: deallocate memory for clur-related variables
  subroutine ctqmc_deallocate_memory_clur()
     implicit none

! local variables
! loop index
     integer :: i

     do i=1,norbs
         call istack_destroy(empty_s(i))
         call istack_destroy(empty_e(i))
     enddo ! over i={1,norbs} loop

     if ( allocated(index_s) ) deallocate(index_s)
     if ( allocated(index_e) ) deallocate(index_e)

     if ( allocated(time_s)  ) deallocate(time_s )
     if ( allocated(time_e)  ) deallocate(time_e )

     if ( allocated(exp_s)   ) deallocate(exp_s  )
     if ( allocated(exp_e)   ) deallocate(exp_e  )

     if ( allocated(empty_s) ) deallocate(empty_s)
     if ( allocated(empty_e) ) deallocate(empty_e)

     return
  end subroutine ctqmc_deallocate_memory_clur

!!>>> ctqmc_deallocate_memory_flvr: deallocate memory for flvr-related variables
  subroutine ctqmc_deallocate_memory_flvr()
     implicit none

     call istack_destroy(empty_v)

     if ( allocated(index_t) ) deallocate(index_t)
     if ( allocated(index_v) ) deallocate(index_v)

     if ( allocated(type_v)  ) deallocate(type_v )
     if ( allocated(flvr_v)  ) deallocate(flvr_v )

     if ( allocated(time_v)  ) deallocate(time_v )

     if ( allocated(expt_t)  ) deallocate(expt_t )
     if ( allocated(expt_v)  ) deallocate(expt_v )

     return
  end subroutine ctqmc_deallocate_memory_flvr

!!>>> ctqmc_deallocate_memory_mesh: deallocate memory for mesh-related variables
  subroutine ctqmc_deallocate_memory_mesh()
     implicit none

     if ( allocated(tmesh) )   deallocate(tmesh)
     if ( allocated(rmesh) )   deallocate(rmesh)

     if ( allocated(pmesh) )   deallocate(pmesh)
     if ( allocated(qmesh) )   deallocate(qmesh)

     if ( allocated(ppleg) )   deallocate(ppleg)
     if ( allocated(qqche) )   deallocate(qqche)

     return
  end subroutine ctqmc_deallocate_memory_mesh

!!>>> ctqmc_deallocate_memory_meat: deallocate memory for meat-related variables
  subroutine ctqmc_deallocate_memory_meat()
     implicit none

     if ( allocated(hist)  )   deallocate(hist )

     if ( allocated(paux)  )   deallocate(paux )
     if ( allocated(prob)  )   deallocate(prob )

     if ( allocated(nmat)  )   deallocate(nmat )
     if ( allocated(nnmat) )   deallocate(nnmat)
     if ( allocated(kmat)  )   deallocate(kmat )
     if ( allocated(kkmat) )   deallocate(kkmat)
     if ( allocated(lmat)  )   deallocate(lmat )
     if ( allocated(rmat)  )   deallocate(rmat )
     if ( allocated(lrmat) )   deallocate(lrmat)

     if ( allocated(g2_re) )   deallocate(g2_re)
     if ( allocated(g2_im) )   deallocate(g2_im)

     if ( allocated(ps_re) )   deallocate(ps_re)
     if ( allocated(ps_im) )   deallocate(ps_im)

     return
  end subroutine ctqmc_deallocate_memory_meat

!!>>> ctqmc_deallocate_memory_umat: deallocate memory for umat-related variables
  subroutine ctqmc_deallocate_memory_umat()
     implicit none

     if ( allocated(rank)  )   deallocate(rank )

     if ( allocated(diag)  )   deallocate(diag )

     if ( allocated(symm)  )   deallocate(symm )

     if ( allocated(eimp)  )   deallocate(eimp )
     if ( allocated(eigs)  )   deallocate(eigs )
     if ( allocated(naux)  )   deallocate(naux )
     if ( allocated(saux)  )   deallocate(saux )

     return
  end subroutine ctqmc_deallocate_memory_umat


!!>>> ctqmc_deallocate_memory_mmat: deallocate memory for mmat-related variables
  subroutine ctqmc_deallocate_memory_mmat()
     implicit none

     if ( allocated(lspace) )  deallocate(lspace)
     if ( allocated(rspace) )  deallocate(rspace)

     if ( allocated(mmat)   )  deallocate(mmat  )

     if ( allocated(lsaves) )  deallocate(lsaves)
     if ( allocated(rsaves) )  deallocate(rsaves)

     if ( allocated(gmat)   )  deallocate(gmat  )

     return
  end subroutine ctqmc_deallocate_memory_mmat

!!>>> ctqmc_deallocate_memory_gmat: deallocate memory for gmat-related variables
  subroutine ctqmc_deallocate_memory_gmat()
     implicit none

     if ( allocated(gtau) )    deallocate(gtau)

     if ( allocated(grnf) )    deallocate(grnf)

     return
  end subroutine ctqmc_deallocate_memory_gmat

!!>>> ctqmc_deallocate_memory_wmat: deallocate memory for wmat-related variables
  subroutine ctqmc_deallocate_memory_wmat()
     implicit none

     if ( allocated(wtau) )    deallocate(wtau)
     if ( allocated(htau) )    deallocate(htau)
     if ( allocated(hsed) )    deallocate(hsed)

     if ( allocated(wssf) )    deallocate(wssf)
     if ( allocated(hybf) )    deallocate(hybf)

     return
  end subroutine ctqmc_deallocate_memory_wmat

!!>>> ctqmc_deallocate_memory_smat: deallocate memory for smat-related variables
  subroutine ctqmc_deallocate_memory_smat()
     implicit none

     if ( allocated(sig1) )    deallocate(sig1)
     if ( allocated(sig2) )    deallocate(sig2)

     return
  end subroutine ctqmc_deallocate_memory_smat

  end module context
