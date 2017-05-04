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
!!! type    : modules
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           05/04/2017 by li huang (last modified)
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

!!
!! @var exp_s
!!
!! exp(i\omega t), s means create operators
!!
     complex(dp), public, save, allocatable :: exp_s(:,:,:)

!!
!! @var exp_e
!!
!! exp(i\omega t), e means destroy operators
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
!! interval [-1,1] on which legendre orthogonal polynomial is defined
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

!!
!! @var prob
!!
!! probability of atomic eigenstates of local hamiltonian matrix
!!
     real(dp), public, save, allocatable :: prob(:)

!!
!! @var paux
!!
!! auxiliary physical observables, it is a vector with size = 9
!!
!! paux(01) : total energy, Etot
!! paux(02) : potential engrgy, Epot
!! paux(03) : kinetic energy, Ekin
!! paux(04) : local magnetic moment, < Sz >
!! paux(05) : average of occupation, < N > = < N^1 > = < N1 >
!! paux(06) : average of occupation square, < N^2 > = < N2 >
!! paux(07) : high order of K, < K^2 > = < K2 >
!! paux(08) : high order of K, < K^3 > = < K3 >
!! paux(09) : high order of K, < K^4 > = < K4 >
!!
!! note: K = current perturbation expansion order * 2. the < K2 >, < K3 >,
!! and < K4 > can be used to calculate the skewness and kurtosis of the
!! perturbation expansion order. of course, < K1 > is essential. It can be
!! calculated from Ekin.
!!
     real(dp), public, save, allocatable :: paux(:)

!!
!! @var nimp
!!
!! impurity occupation number, < n_i >
!!
     real(dp), public, save, allocatable :: nimp(:)

!!
!! @var nmat
!!
!! impurity double occupation number matrix, < n_i n_j >
!!
     real(dp), public, save, allocatable :: nmat(:,:)

!!
!! @var knop
!!
!! number of operators, < k >
!!
     real(dp), public, save, allocatable :: knop(:)

!!
!! @var kmat
!!
!! square of number of operators, < k^2 >
!!
     real(dp), public, save, allocatable :: kmat(:,:)

!!
!! @var lnop
!!
!! number of operators at left half axis, < k_l >
!!
     real(dp), public, save, allocatable :: lnop(:)

!!
!! @var rnop
!!
!! number of operators at right half axis, < k_r >
!!
     real(dp), public, save, allocatable :: rnop(:)

!!
!! @var lrmat
!!
!! used to evaluate fidelity susceptibility, < k_l k_r >
!!
     real(dp), public, save, allocatable :: lrmat(:,:)

!!
!! @var szpow
!!
!! powers of the local magnetization < S^n_z>, which is used to calculate
!! the Binder cumulant
!!
     real(dp), public, save, allocatable :: szpow(:,:)

!!
!! @var schi
!!
!! spin-spin correlation function: < Sz(0) Sz(\tau) >,
!! totally-averaged
!!
     real(dp), public, save, allocatable :: schi(:)

!!
!! @var sschi
!!
!! spin-spin correlation function: < Sz(0) Sz(\tau) >,
!! orbital-resolved
!!
     real(dp), public, save, allocatable :: sschi(:,:)

!!
!! @var ssfom
!!
!! spin-spin correlation function: \chi^{s}_{i} (i\omega),
!! orbital-resolved
!!
     real(dp), public, save, allocatable :: ssfom(:,:)

!!
!! @var ochi
!!
!! orbital-orbital correlation function: < N(0) N(\tau) >,
!! totally-averaged
!!
     real(dp), public, save, allocatable :: ochi(:)

!!
!! @var oochi
!!
!! orbital-orbital correlation function: < N(0) N(\tau) >,
!! orbital-resolved
!!
     real(dp), public, save, allocatable :: oochi(:,:,:)

!!
!! @var oofom
!!
!! orbital-orbital correlation function: \chi^{c}_{ij} (i\omega),
!! orbital-resolved
!!
     real(dp), public, save, allocatable :: oofom(:,:,:)

!!
!! @var g2pw
!!
!! two-particle green's function
!!
     complex(dp), public, save, allocatable :: g2pw(:,:,:,:,:)

!!
!! @var h2pw
!!
!! used to calculate two-particle irreducible vertex function
!!
     complex(dp), public, save, allocatable :: h2pw(:,:,:,:,:)

!!
!! @var p2pw
!!
!! particle-particle pairing susceptibility
!!
     complex(dp), public, save, allocatable :: p2pw(:,:,:,:,:)

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

!!
!! @var rank
!!
!! current perturbation expansion order for different flavor channel
!!
     integer,  public, save, allocatable :: rank(:)

!!
!! @var stts
!!
!! current occupation status for different flavor channel
!!
     integer,  public, save, allocatable :: stts(:)

!!
!! @var pref
!!
!! prefactor for improved estimator for self-energy
!!
     real(dp), public, save, allocatable :: pref(:,:)

!-------------------------------------------------------------------------
!::: input data variables                                              :::
!-------------------------------------------------------------------------

!!
!! @var symm
!!
!! symmetry properties for correlated orbitals
!!
     integer,  public, save, allocatable :: symm(:)

!!
!! @var eimp
!!
!! impurity level for correlated orbitals
!!
     real(dp), public, save, allocatable :: eimp(:)

!!
!! @var ktau
!!
!! screening function, used to determine dynamic interaction, K(\tau)
!!
     real(dp), public, save, allocatable :: ktau(:)

!!
!! @var ksed
!!
!! second order derivates for the screening function, K''(\tau)
!!
     real(dp), public, save, allocatable :: ksed(:)

!!
!! @var ptau
!!
!! first order derivates for the screening function, K'(\tau)
!!
     real(dp), public, save, allocatable :: ptau(:)

!!
!! @var psed
!!
!! second order derivates for ptau, K'''(\tau)
!!
     real(dp), public, save, allocatable :: psed(:)

!!
!! @var uumat
!!
!! reduced Coulomb interaction matrix, two-index version
!!
     real(dp), public, save, allocatable :: uumat(:,:)

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
     real(dp), public, save, allocatable    :: lspace(:,:)

!!
!! @var rspace
!!
!! helper matrix for evaluating M & G matrices
!!
     real(dp), public, save, allocatable    :: rspace(:,:)

!!
!! @var mmat
!!
!! M matrix, $ \mathscr{M} $
!!
     real(dp), public, save, allocatable    :: mmat(:,:,:)

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
     real(dp), public, save, allocatable    :: gtau(:,:,:)

!!
!! @var ftau
!!
!! auxiliary correlation function in imaginary time axis, used to measure
!! self-energy function, F(\tau)
!!
     real(dp), public, save, allocatable    :: ftau(:,:,:)

!!
!! @var grnf
!!
!! impurity green's function in matsubara frequency axis
!!
     complex(dp), public, save, allocatable :: grnf(:,:,:)

!!
!! @var frnf
!!
!! auxiliary correlation function in matsubara frequency axis, used to
!! measure self-energy function, F(i\omega)
!!
     complex(dp), public, save, allocatable :: frnf(:,:,:)

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
     real(dp), public, save, allocatable    :: wtau(:,:,:)

!!
!! @var wssf
!!
!! bath weiss's function in matsubara frequency axis
!!
     complex(dp), public, save, allocatable :: wssf(:,:,:)

!!
!! @var htau
!!
!! hybridization function in imaginary time axis
!!
     real(dp), public, save, allocatable    :: htau(:,:,:)

!!
!! @var hybf
!!
!! hybridization function in matsubara frequency axis
!!
     complex(dp), public, save, allocatable :: hybf(:,:,:)

!!
!! @var hsed
!!
!! second order derivates for hybridization function, it should be used
!! to interpolate htau
!!
     real(dp), public, save, allocatable    :: hsed(:,:,:)

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

!!
!! @mod context
!!
!! containing memory management subroutines and initialize all of the
!! global variables and arrays
!!
  module context
     use constants, only : dp, zero, czero

     use control, only : nband, norbs, ncfgs
     use control, only : lemax, legrd
     use control, only : mkink, mfreq
     use control, only : nffrq, nbfrq
     use control, only : nfreq
     use control, only : ntime

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

!!
!! @sub ctqmc_allocate_memory_clur
!!
!! allocate memory for clur-related variables
!!
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

!!
!! @sub ctqmc_allocate_memory_mesh
!!
!! allocate memory for mesh-related variables
!!
  subroutine ctqmc_allocate_memory_mesh()
     implicit none

! allocate memory
     allocate(tmesh(ntime),       stat=istat)
     allocate(rmesh(mfreq),       stat=istat)

     allocate(lmesh(legrd),       stat=istat)
     allocate(rep_l(legrd,lemax), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_mesh','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     tmesh = zero
     rmesh = zero

     lmesh = zero
     rep_l = zero

     return
  end subroutine ctqmc_allocate_memory_mesh

!!
!! @sub ctqmc_allocate_memory_meat
!!
!! allocate memory for meat-related variables
!!
  subroutine ctqmc_allocate_memory_meat()
     implicit none

! allocate memory
     allocate(hist(mkink),        stat=istat)
     allocate(prob(ncfgs),        stat=istat)
     allocate(paux(  9  ),        stat=istat)
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
     prob  = zero
     paux  = zero
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

!!
!! @sub ctqmc_allocate_memory_umat
!!
!! allocate memory for umat-related variables
!!
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

!!
!! @sub ctqmc_allocate_memory_mmat
!!
!! allocate memory for mmat-related variables
!!
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

!!
!! @sub ctqmc_allocate_memory_gmat
!!
!! allocate memory for gmat-related variables
!!
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

!!
!! @sub ctqmc_allocate_memory_wmat
!!
!! allocate memory for wmat-related variables
!!
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

!!
!! @sub ctqmc_allocate_memory_smat
!!
!! allocate memory for smat-related variables
!!
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

!!
!! @sub ctqmc_deallocate_memory_clur
!!
!! deallocate memory for clur-related variables
!!
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

!!
!! @sub ctqmc_deallocate_memory_mesh
!!
!! deallocate memory for mesh-related variables
!!
  subroutine ctqmc_deallocate_memory_mesh()
     implicit none

     if ( allocated(tmesh) )   deallocate(tmesh)
     if ( allocated(rmesh) )   deallocate(rmesh)

     if ( allocated(lmesh) )   deallocate(lmesh)
     if ( allocated(rep_l) )   deallocate(rep_l)

     return
  end subroutine ctqmc_deallocate_memory_mesh

!!
!! @sub ctqmc_deallocate_memory_meat
!!
!! deallocate memory for meat-related variables
!!
  subroutine ctqmc_deallocate_memory_meat()
     implicit none

     if ( allocated(hist)  )   deallocate(hist )
     if ( allocated(prob)  )   deallocate(prob )
     if ( allocated(paux)  )   deallocate(paux )
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

!!
!! @sub ctqmc_deallocate_memory_umat
!!
!! deallocate memory for umat-related variables
!!
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

!!
!! @sub ctqmc_deallocate_memory_mmat
!!
!! deallocate memory for mmat-related variables
!!
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

!!
!! @sub ctqmc_deallocate_memory_gmat
!!
!! deallocate memory for gmat-related variables
!!
  subroutine ctqmc_deallocate_memory_gmat()
     implicit none

     if ( allocated(gtau) )    deallocate(gtau)
     if ( allocated(ftau) )    deallocate(ftau)

     if ( allocated(grnf) )    deallocate(grnf)
     if ( allocated(frnf) )    deallocate(frnf)

     return
  end subroutine ctqmc_deallocate_memory_gmat

!!
!! @sub ctqmc_deallocate_memory_wmat
!!
!! deallocate memory for wmat-related variables
!!
  subroutine ctqmc_deallocate_memory_wmat()
     implicit none

     if ( allocated(wtau) )    deallocate(wtau)
     if ( allocated(htau) )    deallocate(htau)
     if ( allocated(hsed) )    deallocate(hsed)

     if ( allocated(wssf) )    deallocate(wssf)
     if ( allocated(hybf) )    deallocate(hybf)

     return
  end subroutine ctqmc_deallocate_memory_wmat

!!
!! @sub ctqmc_deallocate_memory_smat
!!
!! deallocate memory for smat-related variables
!!
  subroutine ctqmc_deallocate_memory_smat()
     implicit none

     if ( allocated(sig1) )    deallocate(sig1)
     if ( allocated(sig2) )    deallocate(sig2)

     return
  end subroutine ctqmc_deallocate_memory_smat

  end module context
