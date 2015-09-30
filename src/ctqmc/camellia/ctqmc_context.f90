!!!-----------------------------------------------------------------------
!!! project : camellia
!!! program : ctqmc_core module
!!!           ctqmc_clur module
!!!           ctqmc_flvr module
!!!           ctqmc_mesh module
!!!           ctqmc_meat module
!!!           ctqmc_umat module
!!!           ctqmc_fmat module
!!!           ctqmc_mmat module
!!!           ctqmc_gmat module
!!!           ctqmc_wmat module
!!!           ctqmc_smat module
!!!           context    module
!!! source  : ctqmc_context.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : To define the key data structure and global arrays/variables
!!!           for hybridization expansion version continuous time quantum
!!!           Monte Carlo (CTQMC) quantum impurity solver and dynamical
!!!           mean field theory (DMFT) self-consistent engine
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module ctqmc_core                                                <<<
!!========================================================================

!!>>> containing core (internal) variables used by continuous time quantum
!!>>> Monte Carlo quantum impurity solver
  module ctqmc_core
     use constants, only : dp, zero

     implicit none

! current perturbation expansion order
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

! insert kink (operators pair) statistics: total insert count
     real(dp), public, save :: insert_tcount = zero

! insert kink (operators pair) statistics: total accepted insert count
     real(dp), public, save :: insert_accept = zero

! insert kink (operators pair) statistics: total rejected insert count
     real(dp), public, save :: insert_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, remove action counter                       :::
!-------------------------------------------------------------------------

! remove kink (operators pair) statistics: total remove count
     real(dp), public, save :: remove_tcount = zero

! remove kink (operators pair) statistics: total accepted remove count
     real(dp), public, save :: remove_accept = zero

! remove kink (operators pair) statistics: total rejected remove count
     real(dp), public, save :: remove_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, lshift action counter                       :::
!-------------------------------------------------------------------------

! lshift kink (operators pair) statistics: total lshift count
     real(dp), public, save :: lshift_tcount = zero

! lshift kink (operators pair) statistics: total accepted lshift count
     real(dp), public, save :: lshift_accept = zero

! lshift kink (operators pair) statistics: total rejected lshift count
     real(dp), public, save :: lshift_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, rshift action counter                       :::
!-------------------------------------------------------------------------

! rshift kink (operators pair) statistics: total rshift count
     real(dp), public, save :: rshift_tcount = zero

! rshift kink (operators pair) statistics: total accepted rshift count
     real(dp), public, save :: rshift_accept = zero

! rshift kink (operators pair) statistics: total rejected rshift count
     real(dp), public, save :: rshift_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, reflip action counter                       :::
!-------------------------------------------------------------------------

! reflip kink (operators pair) statistics: total reflip count
     real(dp), public, save :: reflip_tcount = zero

! reflip kink (operators pair) statistics: total accepted reflip count
     real(dp), public, save :: reflip_accept = zero

! reflip kink (operators pair) statistics: total rejected reflip count
     real(dp), public, save :: reflip_reject = zero

  end module ctqmc_core

!!========================================================================
!!>>> module ctqmc_clur                                                <<<
!!========================================================================

!!>>> containing perturbation expansion series related arrays (colour part)
!!>>> used by continuous time quantum Monte Carlo quantum impurity solver
  module ctqmc_clur
     use constants, only : dp
     use stack, only : istack, istack_create, istack_destroy

     implicit none

! memory address index for the imaginary time \tau_s
     integer, public, save, allocatable :: index_s(:,:)

! memory address index for the imaginary time \tau_e
     integer, public, save, allocatable :: index_e(:,:)

! imaginary time \tau_s of create  operators
     real(dp), public, save, allocatable :: time_s(:,:)

! imaginary time \tau_e of destroy operators
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

! imaginary time interval \delta \tau (the last point)
     real(dp), public, save, allocatable :: expt_t(:)

! imaginary time interval \delta \tau
     real(dp), public, save, allocatable :: expt_v(:)

! exp(-H\beta), exponent matrix for local hamiltonian multiply \beta
! it is used to calculate the partition function
     real(dp), public, save, allocatable :: expt_z(:)

  end module ctqmc_flvr

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
     real(dp), public, save, allocatable :: qmesh(:)

! legendre polynomial defined on [-1,1]
     real(dp), public, save, allocatable :: ppleg(:,:)

! chebyshev polynomial defined on [-1,1]
     real(dp), public, save, allocatable :: qqche(:,:)

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
! paux(1) : total energy, Etot
! paux(2) : potential engrgy, Epot
! paux(3) : kinetic energy, Ekin
! paux(4) : magnetic moment, < Sz >
! paux(5) : average of occupation, < N > = < N1 >
! paux(6) : average of occupation square, < N2 >
! paux(7) : reserved
! paux(8) : reserved
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
!!>>> module ctqmc_fmat                                                <<<
!!========================================================================

!!>>> containing F-matrix related arrays used by continuous time quantum
!!>>> Monte Carlo quantum impurity solver
  module ctqmc_fmat
     use constants, only : dp

     implicit none

!-------------------------------------------------------------------------
!::: sparse matrix structure                                           :::
!-------------------------------------------------------------------------

     type T_spmat

! dimension size for original matrix
         integer :: ndim = 0

! maximum number of non-zero elements
         integer :: nval = 0

! row index: element i of it gives the index of the element in the
! vv array that is first non-zero element in a row i
       integer, allocatable  :: iv(:)

! column index: element j of it is the number of the column that contains
! the j-th element in the vv array
       integer, allocatable  :: jv(:)

! a array that contains the non-zero elements for sparse matrix
       real(dp), allocatable :: vv(:)

     end type T_spmat

!-------------------------------------------------------------------------
!::: auxiliary matrix                                                  :::
!-------------------------------------------------------------------------

! eigenvectors of local hamiltonian matrix
     real(dp), public, save, allocatable :: vmat(:,:)

! transposed eigenvectors of local hamiltonian matrix
     real(dp), public, save, allocatable :: wmat(:,:)

! local hamiltonian matrix
     real(dp), public, save, allocatable :: hmat(:,:)

!-------------------------------------------------------------------------
!::: dense matrix style for op                                         :::
!-------------------------------------------------------------------------

! F-matrix <alpha| f^{\dag}_{m} |beta> for create operators
     real(dp), public, save, allocatable :: op_c(:,:,:)

! F-matrix <alpha| f_{m} |beta> for destroy operators
     real(dp), public, save, allocatable :: op_d(:,:,:)

!-------------------------------------------------------------------------
!::: sparse matrix style for op (Compressed Sparse Row (CSR) format)   :::
!-------------------------------------------------------------------------

! spm_c and spm_d are F-matrix, spm_c is for create operator, while spm_d
! is for destroy operator. we need to multiply a series of spm_c, spm_d
! and exponent matrix to get the matrix product trace
     type (T_spmat), public, save, allocatable :: spm_c(:)
     type (T_spmat), public, save, allocatable :: spm_d(:)

! spm_h (hmat) is the local hamiltonian matrix in fock space, it is used
! to calculate exp(tH) |v> by real leja points method
     type (T_spmat), public, save, allocatable :: spm_h(:)

! spm_s is used to calculate matrix product trace efficiently. the final
! matrix product should be stored in spm_s matrix
     type (T_spmat), public, save, allocatable :: spm_s(:)

! spm_n is the precomputed < c^{\dag} c > matrix, it is used to calculate
! impurity occupation number (nmat)
     type (T_spmat), public, save, allocatable :: spm_n(:)

! spm_m is the precomputed < c^{\dag} c c^{\dag} c > matrix, it is used
! to calculate impurity double occupation number (nnmat)
     type (T_spmat), public, save, allocatable :: spm_m(:,:)

  end module ctqmc_fmat

!=========================================================================
!>>> module ctqmc_mmat                                                 <<<
!=========================================================================
!>>> containing M-matrix and G-matrix related arrays used by continuous
! time quantum Monte Carlo quantum impurity solver
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

!=========================================================================
!>>> module ctqmc_gmat                                                 <<<
!=========================================================================
!>>> containing green's function matrix related arrays used by continuous
! time quantum Monte Carlo quantum impurity solver
  module ctqmc_gmat
     use constants, only : dp

     implicit none

! impurity green's function, in imaginary time axis, matrix form
     real(dp), public, save, allocatable    :: gtau(:,:,:)

! impurity green's function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: grnf(:,:,:)

  end module ctqmc_gmat

!=========================================================================
!>>> module ctqmc_wmat                                                 <<<
!=========================================================================
!>>> containing weiss's function and hybridization function matrix related
! arrays used by continuous time quantum Monte Carlo quantum impurity solver
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

!=========================================================================
!>>> module ctqmc_smat                                                 <<<
!=========================================================================
!>>> containing self-energy function matrix related arrays used by
! continuous time quantum Monte Carlo quantum impurity solver
  module ctqmc_smat
     use constants, only : dp

     implicit none

! self-energy function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: sig1(:,:,:)

! self-energy function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: sig2(:,:,:)

  end module ctqmc_smat

!=========================================================================
!>>> module context                                                    <<<
!=========================================================================
!>>> containing memory management subroutines and define global variables
  module context
     use constants
     use control

     use ctqmc_core
     use ctqmc_clur
     use ctqmc_flvr

     use ctqmc_umat
     use ctqmc_fmat
     use ctqmc_mmat

     use ctqmc_gmat
     use ctqmc_wmat
     use ctqmc_smat

     implicit none

! status flag
     integer, private :: istat

! declaration of module procedures: allocate memory
     public :: ctqmc_allocate_memory_clur
     public :: ctqmc_allocate_memory_flvr
     public :: ctqmc_allocate_memory_umat
     public :: ctqmc_allocate_memory_fmat
     public :: ctqmc_allocate_memory_mmat
     public :: ctqmc_allocate_memory_gmat
     public :: ctqmc_allocate_memory_wmat
     public :: ctqmc_allocate_memory_smat

! declaration of module procedures: deallocate memory
     public :: ctqmc_deallocate_memory_clur
     public :: ctqmc_deallocate_memory_flvr
     public :: ctqmc_deallocate_memory_umat
     public :: ctqmc_deallocate_memory_fmat
     public :: ctqmc_deallocate_memory_mmat
     public :: ctqmc_deallocate_memory_gmat
     public :: ctqmc_deallocate_memory_wmat
     public :: ctqmc_deallocate_memory_smat

     contains

!=========================================================================
!>>> allocate memory subroutines                                       <<<
!=========================================================================

!>>> allocate memory for clur-related variables
     subroutine ctqmc_allocate_memory_clur()
         implicit none

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
             call ctqmc_print_error('ctqmc_allocate_memory_clur','can not allocate enough memory')
         endif

! initialize them
         index_s = 0
         index_e = 0

         time_s  = zero
         time_e  = zero

         exp_s   = czero
         exp_e   = czero

         do i=1,norbs
             empty_s(i) = istack_create(mkink)
             empty_e(i) = istack_create(mkink) 
         enddo ! over i={1,norbs} loop

         return
     end subroutine ctqmc_allocate_memory_clur

!>>> allocate memory for flvr-related variables
     subroutine ctqmc_allocate_memory_flvr()
         implicit none

! allocate memory
         allocate(index_t(mkink),      stat=istat)
         allocate(index_v(mkink),      stat=istat)

         allocate(type_v(mkink),       stat=istat)
         allocate(flvr_v(mkink),       stat=istat)

         allocate(time_v(mkink),       stat=istat)

         allocate(expt_t(  2  ),       stat=istat)
         allocate(expt_v(mkink),       stat=istat)
         allocate(expt_z(ncfgs),       stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_flvr','can not allocate enough memory')
         endif

! initialize them
         index_t = 0
         index_v = 0

         type_v  = 1
         flvr_v  = 1

         time_v  = zero

         expt_t  = zero
         expt_v  = zero
         expt_z  = zero

         empty_v = istack_create(mkink)

         return
     end subroutine ctqmc_allocate_memory_flvr

!>>> allocate memory for umat-related variables
     subroutine ctqmc_allocate_memory_umat()
         implicit none

! allocate memory
         allocate(hist(mkink),        stat=istat)
         allocate(rank(norbs),        stat=istat)

         allocate(symm(norbs),        stat=istat)

         allocate(eimp(norbs),        stat=istat)
         allocate(eigs(ncfgs),        stat=istat)
         allocate(naux(ncfgs),        stat=istat)
         allocate(saux(ncfgs),        stat=istat)

         allocate(prob(ncfgs),        stat=istat)
         allocate(paux(  4  ),        stat=istat)
         allocate(schi(ntime),        stat=istat)
         allocate(nmat(norbs),        stat=istat)

         allocate(sschi(ntime,nband), stat=istat)
         allocate(nnmat(norbs,norbs), stat=istat)
         allocate(ddmat(ncfgs,  2  ), stat=istat)

         allocate(ppleg(legrd,lemax), stat=istat)
         allocate(qqche(chgrd,chmax), stat=istat)

         allocate(pmesh(legrd),       stat=istat)
         allocate(qmesh(chgrd),       stat=istat)

         allocate(tmesh(ntime),       stat=istat)
         allocate(rmesh(mfreq),       stat=istat)

         allocate(cmesh(mfreq),       stat=istat)

         allocate(unity(norbs,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_umat','can not allocate enough memory')
         endif

! initialize them
         hist  = 0
         rank  = 0

         symm  = 0

         eimp  = zero
         eigs  = zero
         naux  = zero
         saux  = zero

         prob  = zero
         paux  = zero
         schi  = zero
         nmat  = zero

         sschi = zero
         nnmat = zero
         ddmat = zero

         ppleg = zero
         qqche = zero

         pmesh = zero
         qmesh = zero

         tmesh = zero
         rmesh = zero

         cmesh = czero

         unity = czero

         return
     end subroutine ctqmc_allocate_memory_umat

!>>> allocate memory for fmat-related variables
     subroutine ctqmc_allocate_memory_fmat()
         implicit none

! allocate memory
         allocate(sop_ic(ncfgs+1,norbs),   stat=istat)
         allocate(sop_jc(nfmat,norbs),     stat=istat)
         allocate(sop_c(nfmat,norbs),      stat=istat)

         allocate(sop_id(ncfgs+1,norbs),   stat=istat)
         allocate(sop_jd(nfmat,norbs),     stat=istat)
         allocate(sop_d(nfmat,norbs),      stat=istat)

         allocate(sop_is(ncfgs+1,2),       stat=istat)
         allocate(sop_js(nfmat,2),         stat=istat)
         allocate(sop_s(nfmat,2),          stat=istat)

         allocate(sop_in(ncfgs+1,norbs),   stat=istat)
         allocate(sop_jn(nfmat,norbs),     stat=istat)
         allocate(sop_n(nfmat,norbs),      stat=istat)

         allocate(sop_im(ncfgs+1,norbs,norbs), stat=istat)
         allocate(sop_jm(nfmat,norbs,norbs),   stat=istat)
         allocate(sop_m(nfmat,norbs,norbs),    stat=istat)

         allocate(sop_ih(ncfgs+1),         stat=istat)
         allocate(sop_jh(nhmat),           stat=istat)
         allocate(sop_h(nhmat),            stat=istat)

         allocate(vmat(ncfgs,ncfgs),       stat=istat)
         allocate(wmat(ncfgs,ncfgs),       stat=istat)
         allocate(hmat(ncfgs,ncfgs),       stat=istat)

         allocate(op_c(ncfgs,ncfgs,norbs), stat=istat)
         allocate(op_d(ncfgs,ncfgs,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_fmat','can not allocate enough memory')
         endif

! initialize them
         sop_ic = 0
         sop_jc = 0
         sop_c  = zero

         sop_id = 0
         sop_jd = 0
         sop_d  = zero

         sop_is = 0
         sop_js = 0
         sop_s  = zero

         sop_in = 0
         sop_jn = 0
         sop_n  = zero

         sop_im = 0
         sop_jm = 0
         sop_m  = zero

         sop_ih = 0
         sop_jh = 0
         sop_h  = zero

         vmat   = zero
         wmat   = zero
         hmat   = zero

         op_c   = zero
         op_d   = zero

         return
     end subroutine ctqmc_allocate_memory_fmat

!>>> allocate memory for mmat-related variables
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
             call ctqmc_print_error('ctqmc_allocate_memory_mmat','can not allocate enough memory')
         endif

! initialize them
         lspace = zero
         rspace = zero

         mmat   = zero

         lsaves = czero
         rsaves = czero

         gmat   = czero

         return
     end subroutine ctqmc_allocate_memory_mmat

!>>> allocate memory for gmat-related variables
     subroutine ctqmc_allocate_memory_gmat()
         implicit none

! allocate memory
         allocate(gtau(ntime,norbs,norbs), stat=istat)

         allocate(grnf(mfreq,norbs,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_gmat','can not allocate enough memory')
         endif

! initialize them
         gtau = zero

         grnf = czero

         return
     end subroutine ctqmc_allocate_memory_gmat

!>>> allocate memory for wmat-related variables
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
             call ctqmc_print_error('ctqmc_allocate_memory_wmat','can not allocate enough memory')
         endif

! initialize them
         wtau = zero
         htau = zero
         hsed = zero

         wssf = czero
         hybf = czero

         return
     end subroutine ctqmc_allocate_memory_wmat

!>>> allocate memory for smat-related variables
     subroutine ctqmc_allocate_memory_smat()
         implicit none

! allocate memory
         allocate(sig1(mfreq,norbs,norbs), stat=istat)
         allocate(sig2(mfreq,norbs,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_smat','can not allocate enough memory')
         endif

! initialize them
         sig1 = czero
         sig2 = czero

         return
     end subroutine ctqmc_allocate_memory_smat

!=========================================================================
!>>> deallocate memory subroutines                                     <<<
!=========================================================================

!>>> deallocate memory for clur-related variables
     subroutine ctqmc_deallocate_memory_clur()
         implicit none

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

!>>> deallocate memory for flvr-related variables
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
         if ( allocated(expt_z)  ) deallocate(expt_z )

         return
     end subroutine ctqmc_deallocate_memory_flvr

!>>> deallocate memory for umat-related variables
     subroutine ctqmc_deallocate_memory_umat()
         implicit none

         if ( allocated(hist)  )   deallocate(hist )
         if ( allocated(rank)  )   deallocate(rank )

         if ( allocated(symm)  )   deallocate(symm )

         if ( allocated(eimp)  )   deallocate(eimp )
         if ( allocated(eigs)  )   deallocate(eigs )
         if ( allocated(naux)  )   deallocate(naux )
         if ( allocated(saux)  )   deallocate(saux )

         if ( allocated(prob)  )   deallocate(prob )
         if ( allocated(paux)  )   deallocate(paux )
         if ( allocated(schi)  )   deallocate(schi )
         if ( allocated(nmat)  )   deallocate(nmat )

         if ( allocated(sschi) )   deallocate(sschi)
         if ( allocated(nnmat) )   deallocate(nnmat)
         if ( allocated(ddmat) )   deallocate(ddmat)

         if ( allocated(ppleg) )   deallocate(ppleg)
         if ( allocated(qqche) )   deallocate(qqche)

         if ( allocated(pmesh) )   deallocate(pmesh)
         if ( allocated(qmesh) )   deallocate(qmesh)

         if ( allocated(tmesh) )   deallocate(tmesh)
         if ( allocated(rmesh) )   deallocate(rmesh)

         if ( allocated(cmesh) )   deallocate(cmesh)

         if ( allocated(unity) )   deallocate(unity)

         return
     end subroutine ctqmc_deallocate_memory_umat

!>>> deallocate memory for fmat-related variables
     subroutine ctqmc_deallocate_memory_fmat()
         implicit none

         if ( allocated(sop_ic) )  deallocate(sop_ic)
         if ( allocated(sop_jc) )  deallocate(sop_jc)
         if ( allocated(sop_c)  )  deallocate(sop_c )

         if ( allocated(sop_id) )  deallocate(sop_id)
         if ( allocated(sop_jd) )  deallocate(sop_jd)
         if ( allocated(sop_d)  )  deallocate(sop_d )

         if ( allocated(sop_is) )  deallocate(sop_is)
         if ( allocated(sop_js) )  deallocate(sop_js)
         if ( allocated(sop_s)  )  deallocate(sop_s )

         if ( allocated(sop_in) )  deallocate(sop_in)
         if ( allocated(sop_jn) )  deallocate(sop_jn)
         if ( allocated(sop_n)  )  deallocate(sop_n )

         if ( allocated(sop_im) )  deallocate(sop_im)
         if ( allocated(sop_jm) )  deallocate(sop_jm)
         if ( allocated(sop_m)  )  deallocate(sop_m )

         if ( allocated(sop_ih) )  deallocate(sop_ih)
         if ( allocated(sop_jh) )  deallocate(sop_jh)
         if ( allocated(sop_h)  )  deallocate(sop_h )

         if ( allocated(vmat)   )  deallocate(vmat  )
         if ( allocated(wmat)   )  deallocate(wmat  )
         if ( allocated(hmat)   )  deallocate(hmat  )

         if ( allocated(op_c)   )  deallocate(op_c  )
         if ( allocated(op_d)   )  deallocate(op_d  )

         return
     end subroutine ctqmc_deallocate_memory_fmat

!>>> deallocate memory for mmat-related variables
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

!>>> deallocate memory for gmat-related variables
     subroutine ctqmc_deallocate_memory_gmat()
         implicit none

         if ( allocated(gtau) )    deallocate(gtau)

         if ( allocated(grnf) )    deallocate(grnf)

         return
     end subroutine ctqmc_deallocate_memory_gmat

!>>> deallocate memory for wmat-related variables
     subroutine ctqmc_deallocate_memory_wmat()
         implicit none

         if ( allocated(wtau) )    deallocate(wtau)
         if ( allocated(htau) )    deallocate(htau)
         if ( allocated(hsed) )    deallocate(hsed)

         if ( allocated(wssf) )    deallocate(wssf)
         if ( allocated(hybf) )    deallocate(hybf)

         return
     end subroutine ctqmc_deallocate_memory_wmat

!>>> deallocate memory for smat-related variables
     subroutine ctqmc_deallocate_memory_smat()
         implicit none

         if ( allocated(sig1) )    deallocate(sig1)
         if ( allocated(sig2) )    deallocate(sig2)

         return
     end subroutine ctqmc_deallocate_memory_smat

  end module context
