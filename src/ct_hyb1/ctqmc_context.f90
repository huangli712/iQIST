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
!!!           08/06/2017 by li huang (last modified)
!!! purpose : define the key data structure and global arrays/variables
!!!           for hybridization expansion version continuous time quantum
!!!           Monte Carlo (CTQMC) quantum impurity solver and dynamical
!!!           mean field theory (DMFT) self-consistent engine.
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
     use constants, only : dp
     use constants, only : zero

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
!! if cstat == 0:
!!     it means null occupation case
!!
!! if cstat == 1:
!!     it means partial occupation case (segment scheme)
!!
!! if cstat == 2:
!!     it means partial occupation case (anti-segment scheme)
!!
!! if cstat == 3:
!!     it means full occupation case
!!
     integer, public, save  :: cstat = 0

!-------------------------------------------------------------------------
!::: core variables: real, insert action counter                       :::
!-------------------------------------------------------------------------

!!
!! @var ins_t
!!
!! insert kink (operators pair) statistics: total count
!!
     real(dp), public, save :: ins_t = zero

!!
!! @var ins_a
!!
!! insert kink (operators pair) statistics: accepted count
!!
     real(dp), public, save :: ins_a = zero

!!
!! @var ins_r
!!
!! insert kink (operators pair) statistics: rejected count
!!
     real(dp), public, save :: ins_r = zero

!-------------------------------------------------------------------------
!::: core variables: real, remove action counter                       :::
!-------------------------------------------------------------------------

!!
!! @var rmv_t
!!
!! remove kink (operators pair) statistics: total count
!!
     real(dp), public, save :: rmv_t = zero

!!
!! @var rmv_a
!!
!! remove kink (operators pair) statistics: accepted count
!!
     real(dp), public, save :: rmv_a = zero

!!
!! @var rmv_r
!!
!! remove kink (operators pair) statistics: rejected count
!!
     real(dp), public, save :: rmv_r = zero

!-------------------------------------------------------------------------
!::: core variables: real, lshift action counter                       :::
!-------------------------------------------------------------------------

!!
!! @var lsh_t
!!
!! lshift kink (operators pair) statistics: total count
!!
     real(dp), public, save :: lsh_t = zero

!!
!! @var lsh_a
!!
!! lshift kink (operators pair) statistics: accepted count
!!
     real(dp), public, save :: lsh_a = zero

!!
!! @var lsh_r
!!
!! lshift kink (operators pair) statistics: rejected count
!!
     real(dp), public, save :: lsh_r = zero

!-------------------------------------------------------------------------
!::: core variables: real, rshift action counter                       :::
!-------------------------------------------------------------------------

!!
!! @var rsh_t
!!
!! rshift kink (operators pair) statistics: total count
!!
     real(dp), public, save :: rsh_t = zero

!!
!! @var rsh_a
!!
!! rshift kink (operators pair) statistics: accepted count
!!
     real(dp), public, save :: rsh_a = zero

!!
!! @var rsh_r
!!
!! rshift kink (operators pair) statistics: rejected count
!!
     real(dp), public, save :: rsh_r = zero

!-------------------------------------------------------------------------
!::: core variables: real, reflip action counter                       :::
!-------------------------------------------------------------------------

!!
!! @var rfl_t
!!
!! reflip kink (operators pair) statistics: total count
!!
     real(dp), public, save :: rfl_t = zero

!!
!! @var rfl_a
!!
!! reflip kink (operators pair) statistics: accepted count
!!
     real(dp), public, save :: rfl_a = zero

!!
!! @var rfl_r
!!
!! reflip kink (operators pair) statistics: rejected count
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
!! exp(i\omega \tau_s), s means creation operators
!!
     complex(dp), public, save, allocatable :: exp_s(:,:,:)

!!
!! @var exp_e
!!
!! exp(i\omega \tau_e), e means annihilation operators
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
!! uniform mesh on interval [-1,1] for legendre orthogonal polynomial
!!
     real(dp), public, save, allocatable :: lmesh(:)

!!
!! @var smesh
!!
!! uniform mesh on interval [-1,1] for svd-type orthogonal polynomial
!!
     real(dp), public, save, allocatable :: smesh(:)

!!
!! @var rep_l
!!
!! legendre orthogonal polynomial defined on [-1,1]
!!
     real(dp), public, save, allocatable :: rep_l(:,:)

!!
!! @var rep_s
!!
!! svd-type orthogonal polynomial defined on [-1,1]
!!
     real(dp), public, save, allocatable :: rep_s(:,:)

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
!! @var ac_v
!!
!! a sequence of specified observable which will be used to measure the
!! autocorrelation function. here the total occupation number is the
!! chosen observable
!!
     real(dp), public, save, allocatable :: ac_v(:)

!!
!! @var ac_f
!!
!! autocorrelation function
!!
     real(dp), public, save, allocatable :: ac_f(:)

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
!! auxiliary physical observables, it is a vector with fixed size = 9
!!
!! if p == 01:
!!     paux(p) -> total energy, Etot
!!
!! if p == 02:
!!     paux(p) -> potential engrgy, Epot
!!
!! if p == 03:
!!     paux(p) -> kinetic energy, Ekin
!!
!! if p == 04:
!!     paux(p) -> local magnetic moment, < Sz >
!!
!! if p == 05:
!!     paux(p) -> occupation number, < N > = < N^1 > = < N1 >
!!
!! if p == 06:
!!     paux(p) -> square of occupation number, < N^2 > = < N2 >
!!
!! if p == 07:
!!     paux(p) -> high order of K, < K^2 > = < K2 >
!!
!! if p == 08:
!!     paux(p) -> high order of K, < K^3 > = < K3 >
!!
!! if p == 09:
!!     paux(p) -> high order of K, < K^4 > = < K4 >
!!
!! K = current perturbation expansion order * 2. the < K^2 >, < K^3 >,
!! and < K^4 > can be used to calculate the skewness and kurtosis of the
!! perturbation expansion order. of course, < K^1 > is essential. it can
!! be evaluated from Ekin
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
!! number of operators, < k_i >
!!
     real(dp), public, save, allocatable :: knop(:)

!!
!! @var kmat
!!
!! crossing product of k_i and k_j, < k_i k_j >, i and j \in [1,norbs]
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
!! @var lrmm
!!
!! crossing product of k_l and k_r, < k_l k_r >
!!
     real(dp), public, save, allocatable :: lrmm(:,:)

!!
!! @var szpw
!!
!! powers of the local magnetization, < S^n_z>
!!
     real(dp), public, save, allocatable :: szpw(:,:)

!!
!! @var schi
!!
!! spin-spin correlation function, \chi_{sp}(\tau), totally-averaged
!!
     real(dp), public, save, allocatable :: schi(:)

!!
!! @var sp_t
!!
!! spin-spin correlation function, \chi_{sp}(\tau), orbital-resolved
!!
     real(dp), public, save, allocatable :: sp_t(:,:)

!!
!! @var sp_w
!!
!! spin-spin correlation function, \chi_{sp}(i\omega), orbital-resolved
!!
     real(dp), public, save, allocatable :: sp_w(:,:)

!!
!! @var cchi
!!
!! charge-charge correlation function, \chi_{ch}(\tau), totally-averaged
!!
     real(dp), public, save, allocatable :: cchi(:)

!!
!! @var ch_t
!!
!! charge-charge correlation function, \chi_{ch}(\tau), orbital-resolved
!!
     real(dp), public, save, allocatable :: ch_t(:,:,:)

!!
!! @var ch_w
!!
!! charge-charge correlation function, \chi_{ch}(i\omega), orbital-resolved
!!
     real(dp), public, save, allocatable :: ch_w(:,:,:)

!!
!! @var g2ph
!!
!! two-particle green's function, particle-hole channel
!!
     complex(dp), public, save, allocatable :: g2ph(:,:,:,:,:)

!!
!! @var h2ph
!!
!! used to calculate two-particle vertex function, particle-hole channel
!!
     complex(dp), public, save, allocatable :: h2ph(:,:,:,:,:)

!!
!! @var g2pp
!!
!! two-particle green's function, particle-particle channel
!!
     complex(dp), public, save, allocatable :: g2pp(:,:,:,:,:)

!!
!! @var h2pp
!!
!! used to calculate two-particle vertex function, particle-particle channel
!!
     complex(dp), public, save, allocatable :: h2pp(:,:,:,:,:)

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
     integer, public, save, allocatable  :: rank(:)

!!
!! @var stts
!!
!! current occupation status for different flavor channel
!!
     integer, public, save, allocatable  :: stts(:)

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
     integer, public, save, allocatable  :: symm(:)

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
!! @var umat
!!
!! reduced Coulomb interaction matrix, two-fermions version
!!
     real(dp), public, save, allocatable :: umat(:,:)

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
!! @var ftau
!!
!! auxiliary correlation function in imaginary time axis, used to measure
!! self-energy function, F(\tau)
!!
     real(dp), public, save, allocatable :: ftau(:,:,:)

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

!!
!! @mod context
!!
!! containing memory management subroutines, which initialize all of the
!! global variables and arrays
!!
  module context
     use constants, only : dp
     use constants, only : zero, czero

     use control, only : nband, norbs, ncfgs
     use control, only : lemax, legrd
     use control, only : svmax, svgrd
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
!!>>> declare private variables                                        <<<
!!========================================================================

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

! declaration of module procedures: allocate memory
     public :: cat_alloc_clur
     public :: cat_alloc_mesh
     public :: cat_alloc_meat
     public :: cat_alloc_umat
     public :: cat_alloc_mmat
     public :: cat_alloc_gmat
     public :: cat_alloc_wmat
     public :: cat_alloc_smat

! declaration of module procedures: deallocate memory
     public :: cat_free_clur
     public :: cat_free_mesh
     public :: cat_free_meat
     public :: cat_free_umat
     public :: cat_free_mmat
     public :: cat_free_gmat
     public :: cat_free_wmat
     public :: cat_free_smat

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!
!! @sub cat_alloc_clur
!!
!! allocate memory for clur-related variables
!!
  subroutine cat_alloc_clur()
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
         call s_print_error('cat_alloc_clur','can not allocate enough memory')
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
  end subroutine cat_alloc_clur

!!
!! @sub cat_alloc_mesh
!!
!! allocate memory for mesh-related variables
!!
  subroutine cat_alloc_mesh()
     implicit none

! allocate memory
     allocate(tmesh(ntime),       stat=istat)
     allocate(rmesh(mfreq),       stat=istat)

     allocate(lmesh(legrd),       stat=istat)
     allocate(smesh(svgrd),       stat=istat)
     allocate(rep_l(legrd,lemax), stat=istat)
     allocate(rep_s(svgrd,svmax), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_mesh','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     tmesh = zero
     rmesh = zero

     lmesh = zero
     smesh = zero
     rep_l = zero
     rep_s = zero

     return
  end subroutine cat_alloc_mesh

!!
!! @sub cat_alloc_meat
!!
!! allocate memory for meat-related variables
!!
  subroutine cat_alloc_meat()
     implicit none

! allocate memory
     allocate(ac_v(ntime + 2),   stat=istat)
     allocate(ac_f(ntime + 2),   stat=istat)

     allocate(hist(mkink),       stat=istat)
     allocate(prob(ncfgs),       stat=istat)
     allocate(paux(  9  ),       stat=istat)
     allocate(nimp(norbs),       stat=istat)
     allocate(nmat(norbs,norbs), stat=istat)

     allocate(knop(norbs),       stat=istat)
     allocate(kmat(norbs,norbs), stat=istat)
     allocate(lnop(norbs),       stat=istat)
     allocate(rnop(norbs),       stat=istat)
     allocate(lrmm(norbs,norbs), stat=istat)
     allocate(szpw(  4  ,norbs), stat=istat)

     allocate(schi(ntime),       stat=istat)
     allocate(sp_t(ntime,nband), stat=istat)
     allocate(sp_w(nbfrq,nband), stat=istat)
     allocate(cchi(ntime),       stat=istat)
     allocate(ch_t(ntime,norbs,norbs), stat=istat)
     allocate(ch_w(nbfrq,norbs,norbs), stat=istat)

     allocate(g2ph(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(h2ph(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(g2pp(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)
     allocate(h2pp(nffrq,nffrq,nbfrq,norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_meat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     ac_v = zero
     ac_f = zero

     hist = zero
     prob = zero
     paux = zero
     nimp = zero
     nmat = zero

     knop = zero
     kmat = zero
     lnop = zero
     rnop = zero
     lrmm = zero
     szpw = zero

     schi = zero
     sp_t = zero
     sp_w = zero
     cchi = zero
     ch_t = zero
     ch_w = zero

     g2ph = czero
     h2ph = czero
     g2pp = czero
     h2pp = czero

     return
  end subroutine cat_alloc_meat

!!
!! @sub cat_alloc_umat
!!
!! allocate memory for umat-related variables
!!
  subroutine cat_alloc_umat()
     implicit none

! allocate memory
     allocate(rank(norbs),       stat=istat)
     allocate(stts(norbs),       stat=istat)

     allocate(pref(mkink,norbs), stat=istat)

     allocate(symm(norbs),       stat=istat)

     allocate(eimp(norbs),       stat=istat)
     allocate(ktau(ntime),       stat=istat)
     allocate(ksed(ntime),       stat=istat)
     allocate(ptau(ntime),       stat=istat)
     allocate(psed(ntime),       stat=istat)
     allocate(umat(norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_umat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     rank = 0
     stts = 0

     pref = zero

     symm = 0

     eimp = zero
     ktau = zero
     ksed = zero
     ptau = zero
     psed = zero
     umat = zero

     return
  end subroutine cat_alloc_umat

!!
!! @sub cat_alloc_mmat
!!
!! allocate memory for mmat-related variables
!!
  subroutine cat_alloc_mmat()
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
         call s_print_error('cat_alloc_mmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     lspace = zero
     rspace = zero

     mmat   = zero

     lsaves = czero
     rsaves = czero

     gmat   = czero

     return
  end subroutine cat_alloc_mmat

!!
!! @sub cat_alloc_gmat
!!
!! allocate memory for gmat-related variables
!!
  subroutine cat_alloc_gmat()
     implicit none

! allocate memory
     allocate(gtau(ntime,norbs,norbs), stat=istat)
     allocate(ftau(ntime,norbs,norbs), stat=istat)

     allocate(grnf(mfreq,norbs,norbs), stat=istat)
     allocate(frnf(mfreq,norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_gmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     gtau = zero
     ftau = zero

     grnf = czero
     frnf = czero

     return
  end subroutine cat_alloc_gmat

!!
!! @sub cat_alloc_wmat
!!
!! allocate memory for wmat-related variables
!!
  subroutine cat_alloc_wmat()
     implicit none

! allocate memory
     allocate(wtau(ntime,norbs,norbs), stat=istat)
     allocate(htau(ntime,norbs,norbs), stat=istat)
     allocate(hsed(ntime,norbs,norbs), stat=istat)

     allocate(wssf(mfreq,norbs,norbs), stat=istat)
     allocate(hybf(mfreq,norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_wmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     wtau = zero
     htau = zero
     hsed = zero

     wssf = czero
     hybf = czero

     return
  end subroutine cat_alloc_wmat

!!
!! @sub cat_alloc_smat
!!
!! allocate memory for smat-related variables
!!
  subroutine cat_alloc_smat()
     implicit none

! allocate memory
     allocate(sig1(mfreq,norbs,norbs), stat=istat)
     allocate(sig2(mfreq,norbs,norbs), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_smat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     sig1 = czero
     sig2 = czero

     return
  end subroutine cat_alloc_smat

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!
!! @sub cat_free_clur
!!
!! deallocate memory for clur-related variables
!!
  subroutine cat_free_clur()
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
  end subroutine cat_free_clur

!!
!! @sub cat_free_mesh
!!
!! deallocate memory for mesh-related variables
!!
  subroutine cat_free_mesh()
     implicit none

     if ( allocated(tmesh) )   deallocate(tmesh)
     if ( allocated(rmesh) )   deallocate(rmesh)

     if ( allocated(lmesh) )   deallocate(lmesh)
     if ( allocated(smesh) )   deallocate(smesh)
     if ( allocated(rep_l) )   deallocate(rep_l)
     if ( allocated(rep_s) )   deallocate(rep_s)

     return
  end subroutine cat_free_mesh

!!
!! @sub cat_free_meat
!!
!! deallocate memory for meat-related variables
!!
  subroutine cat_free_meat()
     implicit none

     if ( allocated(ac_v) )    deallocate(ac_v)
     if ( allocated(ac_f) )    deallocate(ac_f)

     if ( allocated(hist) )    deallocate(hist)
     if ( allocated(prob) )    deallocate(prob)
     if ( allocated(paux) )    deallocate(paux)
     if ( allocated(nimp) )    deallocate(nimp)
     if ( allocated(nmat) )    deallocate(nmat)

     if ( allocated(knop) )    deallocate(knop)
     if ( allocated(kmat) )    deallocate(kmat)
     if ( allocated(lnop) )    deallocate(lnop)
     if ( allocated(rnop) )    deallocate(rnop)
     if ( allocated(lrmm) )    deallocate(lrmm)
     if ( allocated(szpw) )    deallocate(szpw)

     if ( allocated(schi) )    deallocate(schi)
     if ( allocated(sp_t) )    deallocate(sp_t)
     if ( allocated(sp_w) )    deallocate(sp_w)
     if ( allocated(cchi) )    deallocate(cchi)
     if ( allocated(ch_t) )    deallocate(ch_t)
     if ( allocated(ch_w) )    deallocate(ch_w)

     if ( allocated(g2ph) )    deallocate(g2ph)
     if ( allocated(h2ph) )    deallocate(h2ph)
     if ( allocated(g2pp) )    deallocate(g2pp)
     if ( allocated(h2pp) )    deallocate(h2pp)

     return
  end subroutine cat_free_meat

!!
!! @sub cat_free_umat
!!
!! deallocate memory for umat-related variables
!!
  subroutine cat_free_umat()
     implicit none

     if ( allocated(rank) )    deallocate(rank)
     if ( allocated(stts) )    deallocate(stts)

     if ( allocated(pref) )    deallocate(pref)

     if ( allocated(symm) )    deallocate(symm)

     if ( allocated(eimp) )    deallocate(eimp)
     if ( allocated(ktau) )    deallocate(ktau)
     if ( allocated(ksed) )    deallocate(ksed)
     if ( allocated(ptau) )    deallocate(ptau)
     if ( allocated(psed) )    deallocate(psed)
     if ( allocated(umat) )    deallocate(umat)

     return
  end subroutine cat_free_umat

!!
!! @sub cat_free_mmat
!!
!! deallocate memory for mmat-related variables
!!
  subroutine cat_free_mmat()
     implicit none

     if ( allocated(lspace) )  deallocate(lspace)
     if ( allocated(rspace) )  deallocate(rspace)

     if ( allocated(mmat)   )  deallocate(mmat  )

     if ( allocated(lsaves) )  deallocate(lsaves)
     if ( allocated(rsaves) )  deallocate(rsaves)

     if ( allocated(gmat)   )  deallocate(gmat  )

     return
  end subroutine cat_free_mmat

!!
!! @sub cat_free_gmat
!!
!! deallocate memory for gmat-related variables
!!
  subroutine cat_free_gmat()
     implicit none

     if ( allocated(gtau) )    deallocate(gtau)
     if ( allocated(ftau) )    deallocate(ftau)

     if ( allocated(grnf) )    deallocate(grnf)
     if ( allocated(frnf) )    deallocate(frnf)

     return
  end subroutine cat_free_gmat

!!
!! @sub cat_free_wmat
!!
!! deallocate memory for wmat-related variables
!!
  subroutine cat_free_wmat()
     implicit none

     if ( allocated(wtau) )    deallocate(wtau)
     if ( allocated(htau) )    deallocate(htau)
     if ( allocated(hsed) )    deallocate(hsed)

     if ( allocated(wssf) )    deallocate(wssf)
     if ( allocated(hybf) )    deallocate(hybf)

     return
  end subroutine cat_free_wmat

!!
!! @sub cat_free_smat
!!
!! deallocate memory for smat-related variables
!!
  subroutine cat_free_smat()
     implicit none

     if ( allocated(sig1) )    deallocate(sig1)
     if ( allocated(sig2) )    deallocate(sig2)

     return
  end subroutine cat_free_smat

  end module context
