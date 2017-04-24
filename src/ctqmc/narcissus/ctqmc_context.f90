!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_core module
!!!           ctqmc_clur module
!!!           ctqmc_mesh module
!!!           ctqmc_meat module
!!!           ctqmc_umat module
!!!           ctqmc_mmat module
!!!           ctqmc_gmat module
!!!           ctqmc_wmat module
!!!           ctqmc_smat module
!!!           context    module
!!! source  : ctqmc_context.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           04/24/2017 by li huang (last modified)
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

!!
!! @var cstat
!!
!! status indicator for the current flavor, used to sync with stts vector
!!
!! if cstat = 0:
!!     it means null occupation case
!!
!! if cstat = 1:
!!     it means partial occupation case, segment scheme
!!
!! if cstat = 2:
!!     it means partial occupation case, anti-segment scheme
!!
!! if cstat = 3:
!!     it means full occupation case
!!
     integer, public, save  :: cstat = 0

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

!!
!! @var rfl_t
!!
!! reflip kink (operators pair) statistics: total reflip count
!!
     real(dp), public, save :: rfl_t = zero

!!
!! @var rfl_a
!!
!! reflip kink (operators pair) statistics: total accepted reflip count
!!
     real(dp), public, save :: rfl_a = zero

!!
!! @var rfl_r
!!
!! reflip kink (operators pair) statistics: total rejected reflip count
!!
     real(dp), public, save :: rfl_r = zero

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
     integer, public, save, allocatable :: index_s(:,:)

!!
!! @var index_e
!!
!! memory address index for the imaginary time \tau_e
!!
     integer, public, save, allocatable :: index_e(:,:)

!!
!! @var time_s
!!
!! imaginary time \tau_s of create operators
!!
     real(dp), public, save, allocatable :: time_s(:,:)

!!
!! @var time_e
!!
!! imaginary time \tau_e of destroy operators
!!
     real(dp), public, save, allocatable :: time_e(:,:)

! exp(i\omega t), s means create  operators
     complex(dp), public, save, allocatable :: exp_s(:,:,:)

! exp(i\omega t), e means destroy operators
     complex(dp), public, save, allocatable :: exp_e(:,:,:)

! container for the empty (unoccupied) memory address index
     type (istack), public, save, allocatable :: empty_s(:)

! container for the empty (unoccupied) memory address index
     type (istack), public, save, allocatable :: empty_e(:)

  end module ctqmc_clur

!!========================================================================
!!>>> module ctqmc_mesh                                                <<<
!!========================================================================

!!>>> containing mesh related arrays used by continuous time quantum Monte
!!>>> Carlo quantum impurity solver
  module ctqmc_mesh
     use constants, only : dp

     implicit none

! imaginary time mesh
     real(dp), public, save, allocatable :: tmesh(:)

! real matsubara frequency mesh
     real(dp), public, save, allocatable :: rmesh(:)

! interval [-1,1] on which legendre polynomial is defined
     real(dp), public, save, allocatable :: pmesh(:)

! interval [-1,1] on which chebyshev polynomial is defined
     !!real(dp), public, save, allocatable :: qmesh(:)

! legendre polynomial defined on [-1,1]
     real(dp), public, save, allocatable :: ppleg(:,:)

! chebyshev polynomial defined on [-1,1]
     !!real(dp), public, save, allocatable :: qqche(:,:)

  end module ctqmc_mesh

!!========================================================================
!!>>> module ctqmc_meat                                                <<<
!!========================================================================

!!>>> containing physical observables related arrays used by continuous
!!>>> time quantum Monte Carlo quantum impurity solver
  module ctqmc_meat !!>>> To tell you a truth, meat means MEAsuremenT
     use constants, only : dp

     implicit none

! histogram for perturbation expansion series
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

! powers of the local magnetization < S^n_z>, used to calculate Binder cumulant
     real(dp), public, save, allocatable :: szpow(:,:)

! spin-spin correlation function: < Sz(0) Sz(\tau) >, \chi_{loc}, totally-averaged
     real(dp), public, save, allocatable :: schi(:)

! spin-spin correlation function: < Sz(0) Sz(\tau) >, \chi_{loc}, orbital-resolved
     real(dp), public, save, allocatable :: sschi(:,:)

! spin-spin correlation function: \chi^{s}_{i} (i\omega), orbital-resolved
     real(dp), public, save, allocatable :: ssfom(:,:)

! orbital-orbital correlation function: < N(0) N(\tau) >, totally-averaged
     real(dp), public, save, allocatable :: ochi(:)

! orbital-orbital correlation function: < N(0) N(\tau) >, orbital-resolved
     real(dp), public, save, allocatable :: oochi(:,:,:)

! orbital-orbital correlation function: \chi^{c}_{ij} (i\omega), orbital-resolved
     real(dp), public, save, allocatable :: oofom(:,:,:)

! used to calculate two-particle green's function, real part
     real(dp), public, save, allocatable :: g2_re(:,:,:,:,:)

! used to calculate two-particle green's function, imaginary part
     real(dp), public, save, allocatable :: g2_im(:,:,:,:,:)

! used to calculate two-particle green's function, real part
     real(dp), public, save, allocatable :: h2_re(:,:,:,:,:)

! used to calculate two-particle green's function, imaginary part
     real(dp), public, save, allocatable :: h2_im(:,:,:,:,:)

! particle-particle pair susceptibility, real part
     real(dp), public, save, allocatable :: ps_re(:,:,:,:,:)

! particle-particle pair susceptibility, imaginary part
     real(dp), public, save, allocatable :: ps_im(:,:,:,:,:)

  end module ctqmc_meat

!!========================================================================
!!>>> module ctqmc_umat                                                <<<
!!========================================================================

!!>>> containing auxiliary arrays used by continuous time quantum Monte
!!>>> Carlo quantum impurity solver
  module ctqmc_umat
     use constants, only : dp

     implicit none

!-------------------------------------------------------------------------
!::: ctqmc status variables                                            :::
!-------------------------------------------------------------------------

! current perturbation expansion order for different flavor channel
     integer,  public, save, allocatable :: rank(:)

! current occupation status for different flavor channel
     integer,  public, save, allocatable :: stts(:)

! prefactor for improved estimator
     real(dp), public, save, allocatable :: pref(:,:)

!-------------------------------------------------------------------------
!::: input data variables                                              :::
!-------------------------------------------------------------------------

! symmetry properties for correlated orbitals
     integer,  public, save, allocatable :: symm(:)

! impurity level for correlated orbitals
     real(dp), public, save, allocatable :: eimp(:)

! screening function, used to measure dynamical screening effect, K(\tau)
     real(dp), public, save, allocatable :: ktau(:)

! second order derivates for the screening function, K''(\tau)
     real(dp), public, save, allocatable :: ksed(:)

! first  order derivates for the screening function, K'(\tau)
     real(dp), public, save, allocatable :: ptau(:)

! second order derivates for ptau, K'''(\tau)
     real(dp), public, save, allocatable :: psed(:)

! reduced Coulomb interaction matrix, two-index version
     real(dp), public, save, allocatable :: uumat(:,:)

  end module ctqmc_umat

!!========================================================================
!!>>> module ctqmc_mmat                                                <<<
!!========================================================================

!!>>> containing M-matrix and G-matrix related arrays used by continuous
!!>>> time quantum Monte Carlo quantum impurity solver
  module ctqmc_mmat
     use constants, only : dp

     implicit none

! helper matrix for evaluating M & G matrices
     real(dp), public, save, allocatable    :: lspace(:,:)

! helper matrix for evaluating M & G matrices
     real(dp), public, save, allocatable    :: rspace(:,:)

! M matrix, $ \mathscr{M} $
     real(dp), public, save, allocatable    :: mmat(:,:,:)

! helper matrix for evaluating G matrix
     complex(dp), public, save, allocatable :: lsaves(:,:)

! helper matrix for evaluating G matrix
     complex(dp), public, save, allocatable :: rsaves(:,:)

! G matrix, $ \mathscr{G} $
     complex(dp), public, save, allocatable :: gmat(:,:,:)

  end module ctqmc_mmat

!!========================================================================
!!>>> module ctqmc_gmat                                                <<<
!!========================================================================

!!>>> containing green's function matrix related arrays used by continuous
!!>>> time quantum Monte Carlo quantum impurity solver
  module ctqmc_gmat
     use constants, only : dp

     implicit none

! impurity green's function, in imaginary time axis, matrix form
     real(dp), public, save, allocatable    :: gtau(:,:,:)

! auxiliary correlation function, in imaginary time axis, matrix form
! used to measure self-energy function, F(\tau)
     real(dp), public, save, allocatable    :: ftau(:,:,:)

! impurity green's function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: grnf(:,:,:)

! auxiliary correlation function, in matsubara frequency axis, matrix form
! used to measure self-energy function, F(i\omega)
     complex(dp), public, save, allocatable :: frnf(:,:,:)

  end module ctqmc_gmat

!!========================================================================
!!>>> module ctqmc_wmat                                                <<<
!!========================================================================

!!>>> containing weiss's function and hybridization function matrix related
!!>>> arrays used by continuous time quantum Monte Carlo quantum impurity
!!>>> solver
  module ctqmc_wmat
     use constants, only : dp

     implicit none

! bath weiss's function, in imaginary time axis, matrix form
     real(dp), public, save, allocatable    :: wtau(:,:,:)

! bath weiss's function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: wssf(:,:,:)

! hybridization function, in imaginary time axis, matrix form
     real(dp), public, save, allocatable    :: htau(:,:,:)

! hybridization function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: hybf(:,:,:)

! second order derivates for hybridization function, used to interpolate htau
     real(dp), public, save, allocatable    :: hsed(:,:,:)

  end module ctqmc_wmat

!!========================================================================
!!>>> module ctqmc_smat                                                <<<
!!========================================================================

!!>>> containing self-energy function matrix related arrays used by
!!>>> continuous time quantum Monte Carlo quantum impurity solver
  module ctqmc_smat
     use constants, only : dp

     implicit none

! self-energy function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: sig1(:,:,:)

! self-energy function, in matsubara frequency axis, matrix form
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
     public :: ctqmc_allocate_memory_mesh
     public :: ctqmc_allocate_memory_meat
     public :: ctqmc_allocate_memory_umat
     public :: ctqmc_allocate_memory_mmat
     public :: ctqmc_allocate_memory_gmat
     public :: ctqmc_allocate_memory_wmat
     public :: ctqmc_allocate_memory_smat

! declaration of module procedures: deallocate memory
     public :: ctqmc_deallocate_memory_clur
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

!!>>> ctqmc_allocate_memory_mesh: allocate memory for mesh-related variables
  subroutine ctqmc_allocate_memory_mesh()
     implicit none

! allocate memory
     allocate(tmesh(ntime),       stat=istat)
     allocate(rmesh(mfreq),       stat=istat)

     allocate(pmesh(legrd),       stat=istat)
     !!allocate(qmesh(chgrd),       stat=istat)

     allocate(ppleg(legrd,lemax), stat=istat)
     !!allocate(qqche(chgrd,chmax), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_mesh','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     tmesh = zero
     rmesh = zero

     pmesh = zero
     !!qmesh = zero

     ppleg = zero
     !!qqche = zero

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
     allocate(szpow(  4  ,norbs), stat=istat)
     allocate(schi(ntime),        stat=istat)
     allocate(sschi(ntime,nband), stat=istat)
     allocate(ssfom(nbfrq,nband), stat=istat)
     allocate(ochi(ntime),        stat=istat)
     allocate(oochi(ntime,norbs,norbs), stat=istat)
     allocate(oofom(nbfrq,norbs,norbs), stat=istat)

     allocate(g2_re(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(g2_im(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(h2_re(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(h2_im(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
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
     szpow = zero
     schi  = zero
     sschi = zero
     ssfom = zero
     ochi  = zero
     oochi = zero
     oofom = zero

     g2_re = zero
     g2_im = zero
     h2_re = zero
     h2_im = zero
     ps_re = zero
     ps_im = zero

     return
  end subroutine ctqmc_allocate_memory_meat

!!>>> ctqmc_allocate_memory_umat: allocate memory for umat-related variables
  subroutine ctqmc_allocate_memory_umat()
     implicit none

! allocate memory
     allocate(rank(norbs),        stat=istat)
     allocate(stts(norbs),        stat=istat)

     allocate(pref(mkink,norbs),  stat=istat)

     allocate(symm(norbs),        stat=istat)

     allocate(eimp(norbs),        stat=istat)
     allocate(ktau(ntime),        stat=istat)
     allocate(ksed(ntime),        stat=istat)
     allocate(ptau(ntime),        stat=istat)
     allocate(psed(ntime),        stat=istat)
     allocate(uumat(norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_umat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     rank  = 0
     stts  = 0

     pref  = zero

     symm  = 0

     eimp  = zero
     ktau  = zero
     ksed  = zero
     ptau  = zero
     psed  = zero
     uumat = zero

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
     allocate(ftau(ntime,norbs,norbs), stat=istat)

     allocate(grnf(mfreq,norbs,norbs), stat=istat)
     allocate(frnf(mfreq,norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_gmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     gtau = zero
     ftau = zero

     grnf = czero
     frnf = czero

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

!!>>> ctqmc_deallocate_memory_mesh: deallocate memory for mesh-related variables
  subroutine ctqmc_deallocate_memory_mesh()
     implicit none

     if ( allocated(tmesh) )   deallocate(tmesh)
     if ( allocated(rmesh) )   deallocate(rmesh)

     if ( allocated(pmesh) )   deallocate(pmesh)
     !!if ( allocated(qmesh) )   deallocate(qmesh)

     if ( allocated(ppleg) )   deallocate(ppleg)
     !!if ( allocated(qqche) )   deallocate(qqche)

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
     if ( allocated(szpow) )   deallocate(szpow)
     if ( allocated(schi)  )   deallocate(schi )
     if ( allocated(sschi) )   deallocate(sschi)
     if ( allocated(ssfom) )   deallocate(ssfom)
     if ( allocated(ochi)  )   deallocate(ochi )
     if ( allocated(oochi) )   deallocate(oochi)
     if ( allocated(oofom) )   deallocate(oofom)

     if ( allocated(g2_re) )   deallocate(g2_re)
     if ( allocated(g2_im) )   deallocate(g2_im)
     if ( allocated(h2_re) )   deallocate(h2_re)
     if ( allocated(h2_im) )   deallocate(h2_im)
     if ( allocated(ps_re) )   deallocate(ps_re)
     if ( allocated(ps_im) )   deallocate(ps_im)

     return
  end subroutine ctqmc_deallocate_memory_meat

!!>>> ctqmc_deallocate_memory_umat: deallocate memory for umat-related variables
  subroutine ctqmc_deallocate_memory_umat()
     implicit none

     if ( allocated(rank)  )   deallocate(rank )
     if ( allocated(stts)  )   deallocate(stts )

     if ( allocated(pref)  )   deallocate(pref )

     if ( allocated(symm)  )   deallocate(symm )

     if ( allocated(eimp)  )   deallocate(eimp )
     if ( allocated(ktau)  )   deallocate(ktau )
     if ( allocated(ksed)  )   deallocate(ksed )
     if ( allocated(ptau)  )   deallocate(ptau )
     if ( allocated(psed)  )   deallocate(psed )
     if ( allocated(uumat) )   deallocate(uumat)

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
     if ( allocated(ftau) )    deallocate(ftau)

     if ( allocated(grnf) )    deallocate(grnf)
     if ( allocated(frnf) )    deallocate(frnf)

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
