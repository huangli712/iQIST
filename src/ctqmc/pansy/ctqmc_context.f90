!!!-----------------------------------------------------------------------
!!! project : pansy
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
!!!           m_sect     module
!!!           m_part     module
!!! source  : ctqmc_context.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!!           yilin wang (email:qhwyl2006@126.com)
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

! exp(-H\tau), exponent matrix for local hamiltonian multiply \tau (the last point)
     real(dp), public, save, allocatable :: expt_t(:,:)

! exp(-H\tau), exponent matrix for local hamiltonian multiply \tau
     real(dp), public, save, allocatable :: expt_v(:,:)

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
! paux(5) : average of occupation, < N > = < N^1 > = < N1 >
! paux(6) : average of occupation square, < N^2 > = < N2 >
! paux(7) : K = current perturbation expansion order X 2, < K^2 > = < K2 >
! paux(8) : K = current perturbation expansion order X 2, < K^3 > = < K3 >
! paux(9) : K = current perturbation expansion order X 2, < K^4 > = < K4 >
     real(dp), public, save, allocatable :: paux(:)

! probability of eigenstates of local hamiltonian matrix
     real(dp), public, save, allocatable :: prob(:)

! impurity occupation number, < n_i >
     real(dp), public, save, allocatable :: nmat(:)

! impurity double occupation number matrix, < n_i n_j >
     real(dp), public, save, allocatable :: nnmat(:,:)

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

! impurity green's function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: grnf(:,:,:)

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

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_mesh','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     tmesh = zero
     rmesh = zero

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




!!========================================================================
!!>>> module m_sect                                                    <<<
!!========================================================================

!!>>> define the data structure for good quantum numbers (GQNs) algorithm
  module m_sect
     use constants, only : dp, zero

     use control, only : norbs
     use control, only : mkink
     use context, only : type_v, flvr_v

     implicit none

!!========================================================================
!!>>> declare global structures                                        <<<
!!========================================================================

! data structure for one F-matrix
!-------------------------------------------------------------------------
     public :: t_fmat
     type t_fmat

! the dimension, n x m
         integer :: n
         integer :: m

! the memory space for the matrix
         real(dp), allocatable :: val(:,:)

     end type t_fmat

! data structure for one sector
!-------------------------------------------------------------------------
     public :: t_sector
     type t_sector

! number of states in this sector
         integer :: ndim

! number of fermion operators, it should be equal to norbs
         integer :: nops

! start index of this sector
         integer :: istart

! total number of electrons
         integer :: nele

! z component of spin: Sz
         integer :: sz

! z component of spin-orbit momentum: Jz
         integer :: jz

! PS good quantum number
         integer :: ps

! the next sector when a fermion operator acts on the sector
! next(nops,0) for annihilation and next(nops,1) for creation operators
! -1 means it is outside the Hilbert space,
! otherwise, it is the index of next sector
         integer, allocatable  :: next(:,:)

! the eigenvalues
         real(dp), allocatable :: eval(:)

! final products of matrices
         real(dp), allocatable :: prod(:)

! the F-matrix between this sector and all other sectors
! fmat(nops,0) for annihilation and fmat(nops,1) for creation operators
! if this sector doesn't point to some other sectors, it is not allocated
         type (t_fmat), allocatable :: fmat(:,:)

     end type t_sector

!!========================================================================
!!>>> declare global variables                                         <<<
!!========================================================================

! total number of sectors
     integer, public, save  :: nsect

! maximal dimension of the sectors
     integer, public, save  :: max_dim_sect

! average dimension of the sectors
     real(dp), public, save :: ave_dim_sect

! array of t_sector contains all the sectors
     type (t_sector), public, save, allocatable :: sectors(:)

!!========================================================================
!!>>> declare private variables                                        <<<
!!========================================================================

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: ctqmc_allocate_memory_one_fmat
     public :: ctqmc_allocate_memory_one_sect
     public :: ctqmc_allocate_memory_sect

     public :: ctqmc_deallocate_memory_one_fmat
     public :: ctqmc_deallocate_memory_one_sect
     public :: ctqmc_deallocate_memory_sect

     public :: cat_make_string

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!>>> ctqmc_allocate_memory_one_fmat: allocate memory for one F-matrix
  subroutine ctqmc_allocate_memory_one_fmat(mat)
     implicit none

! external variables
! F-matrix structure
     type (t_fmat), intent(inout) :: mat

! allocate memory
     allocate(mat%val(mat%n,mat%m), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_one_fmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize it
     mat%val = zero

     return
  end subroutine ctqmc_allocate_memory_one_fmat

!!>>> ctqmc_allocate_memory_one_sect: allocate memory for one sector
  subroutine ctqmc_allocate_memory_one_sect(sect)
     implicit none

! external variables
! sector structure
     type (t_sector), intent(inout) :: sect

! local variables
! loop index
     integer :: i
     integer :: j

! allocate memory
     allocate(sect%next(sect%nops,0:1), stat=istat)

     allocate(sect%eval(sect%ndim),     stat=istat)
     allocate(sect%prod(sect%ndim),     stat=istat)

     allocate(sect%fmat(sect%nops,0:1), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_one_sect','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     sect%next = 0

     sect%eval = zero
     sect%prod = zero

! initialize fmat one by one
     do i=1,sect%nops
         do j=0,1
             sect%fmat(i,j)%n = 0
             sect%fmat(i,j)%m = 0
         enddo ! over j={0,1} loop
     enddo ! over i={1,sect%nops} loop

     return
  end subroutine ctqmc_allocate_memory_one_sect

!!>>> ctqmc_allocate_memory_sect: allocate memory for sector related variables
  subroutine ctqmc_allocate_memory_sect()
     implicit none

! local variables
! loop index
     integer :: i

! allocate memory
     allocate(sectors(nsect), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_sect','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     do i=1,nsect
         sectors(i)%ndim   = 0
         sectors(i)%nops   = norbs
         sectors(i)%istart = 0
         sectors(i)%nele   = 0
         sectors(i)%sz     = 0
         sectors(i)%jz     = 0
         sectors(i)%ps     = 0
     enddo ! over i={1,nsect} loop

     return
  end subroutine ctqmc_allocate_memory_sect

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!>>> ctqmc_deallocate_memory_one_fmat: deallocate memory for one F-matrix
  subroutine ctqmc_deallocate_memory_one_fmat(mat)
     implicit none

! external variables
! F-matrix structure
     type (t_fmat), intent(inout) :: mat

     if ( allocated(mat%val) ) deallocate(mat%val)

     return
  end subroutine ctqmc_deallocate_memory_one_fmat

!!>>> ctqmc_deallocate_memory_one_sect: deallocate memory for one sector
  subroutine ctqmc_deallocate_memory_one_sect(sect)
     implicit none

! external variables
! sector structure
     type (t_sector), intent(inout) :: sect

! local variables
! loop index
     integer :: i
     integer :: j

     if ( allocated(sect%next) ) deallocate(sect%next)

     if ( allocated(sect%eval) ) deallocate(sect%eval)
     if ( allocated(sect%prod) ) deallocate(sect%prod)

! deallocate fmat one by one
     if ( allocated(sect%fmat) ) then
         do i=1,sect%nops
             do j=0,1
                 call ctqmc_deallocate_memory_one_fmat(sect%fmat(i,j))
             enddo ! over j={0,1} loop
         enddo ! over i={1,sect%nops} loop
         deallocate(sect%fmat)
     endif ! back if ( allocated(sect%fmat) ) block

     return
  end subroutine ctqmc_deallocate_memory_one_sect

!!>>> ctqmc_deallocate_memory_sect: deallocate memory for sector
!!>>> related variables
  subroutine ctqmc_deallocate_memory_sect()
     implicit none

! local variables
! loop index
     integer :: i

! first, loop over all the sectors and deallocate their component's memory
! then, deallocate memory of the sectors itself
     if ( allocated(sectors) ) then
         do i=1,nsect
             call ctqmc_deallocate_memory_one_sect(sectors(i))
         enddo ! over i={1,nsect} loop
         deallocate(sectors)
     endif ! back if ( allocated(sectors) ) block

     return
  end subroutine ctqmc_deallocate_memory_sect

!!========================================================================
!!>>> core service subroutines                                         <<<
!!========================================================================

!!>>> cat_make_string: it is used to build a time evolution string
  subroutine cat_make_string(csize, vindex, string)
     implicit none

! external variables
! number of fermion operators for the current diagram
     integer, intent(in)  :: csize

! memory address index of fermion operators
     integer, intent(in)  :: vindex(mkink)

! time evolution string, i.e., sequence of sector index
! if it is not a valid string, then all of its values should be -1
     integer, intent(out) :: string(csize+1,nsect)

! local variables
! loop index
     integer :: i
     integer :: j

! flavor and type of fermion operators
     integer :: vf
     integer :: vt

! current sector index and next sector index
     integer :: curr_sect
     integer :: next_sect

! init return array, we assume all of strings are invalid
     string = -1

! we try to build a string from left to right, that is, 0 -> \beta
! we assume the sectors are S1, S2, S3, ..., SM, and the fermion
! operators are F1, F2, F3, F4, .... FN. here, F1 is in \tau_1, F2
! is in \tau_2, F3 is in \tau_3, and so on, and 
!     0 < \tau_1 < \tau_2 < \tau_3 < ... < \beta
! is always guaranteed. then a typical (and also valid) string must
! look like this:
!     F1       F2       F3       F4       F5        FN
! S1 ----> S2 ----> S3 ----> S4 ----> S5 ----> ... ----> S1
! then the sequence of sector indices is the so-called string. if some
! Si are -1 (null sector), this string is invalid. we will enforce all
! elements in it to be -1. it is easy to speculate that if the number
! of fermion operators is csize, the length of string must be csize + 1
     SECTOR_SCAN_LOOP: do i=1,nsect
! setup starting sector
         curr_sect = i
         string(1,i) = curr_sect
         OPERATOR_SCAN_LOOP: do j=1,csize
! determine the type and flavor of current operator
             vt = type_v( vindex(j) )
             vf = flvr_v( vindex(j) )
! get the next sector
             next_sect = sectors(curr_sect)%next(vf,vt)
! meet null sector, it is an invalid string. we will try another
! new string
             if ( next_sect == -1 ) then
                 string(:,i) = -1; EXIT OPERATOR_SCAN_LOOP
! the string is still alive, we record the sector, and set it to
! the current sector
             else
                 string(j+1,i) = next_sect
                 curr_sect = next_sect
             endif ! back if ( next_sect == -1 ) block
         enddo OPERATOR_SCAN_LOOP ! over j={1,csize} loop
! we have to ensure that the first sector is the same with the last
! sector in this string, or else it is invalid
         if ( string(1,i) /= string(csize+1,i) ) then
             string(:,i) = -1
         endif ! back if ( string(1,i) /= string(csize+1,i) ) block
     enddo SECTOR_SCAN_LOOP ! over i={1,nsect} loop

     return
  end subroutine cat_make_string

  end module m_sect




!!========================================================================
!!>>> module m_part                                                    <<<
!!========================================================================

!!>>> contains some key global variables and subroutines for divide and
!!>>> conquer algorithm to speed up the trace evaluation
  module m_part
     use constants, only : dp, zero, one

     use control, only : ncfgs
     use control, only : mkink
     use control, only : npart
     use control, only : beta
     use context, only : type_v, flvr_v, time_v, expt_v

     use m_sect, only : nsect, max_dim_sect
     use m_sect, only : sectors

     implicit none

!!========================================================================
!!>>> declare global variables                                         <<<
!!========================================================================

! number of operators for each part
     integer, public, save, allocatable  :: nop(:)

! start index of operators for each part
     integer, public, save, allocatable  :: ops(:)

! end index of operators for each part
     integer, public, save, allocatable  :: ope(:)

! how to treat each part when calculating trace
! 0: matrices product for this part has been calculated previously
! 1: this part should be recalculated, and the result must be
!    stored in saved_p, if this Monte Caro move has been accepted
     integer, public, save, allocatable  :: renew(:)

! determine which parts of saved_p are unsafe or invalid (we just call
! it asynchronization), and have to be updated (or synchronized) for
! future trace calculations
! 0: synchronous, this part of saved_p is OK
! 1: asynchronous, this part of saved_p is invalid
! Q: why is renew not enough? why do we need async and is_cp?
! A: because string is not always valid. string broken is possible. at
! that time, even renew(j) is 1, some sectors in this part will be not
! updated successfully. of course, saved_p for them will be not updated
! as well. so we have to mark the corresponding saved_p as wrong value.
! this is the role of async. due to the same reason, we cann't use renew
! to control which parts of saved_p should be updated with saved_n only.
! so we need is_cp as well.
     integer, public, save, allocatable  :: async(:,:)

! determine which parts of saved_p should be updated by the corresponding
! parts of saved_n
! 0: do nothing, saved_n and saved_p have the same values, or saved_n is
!    unavailable, we can not use it to update saved_p
! 1: saved_n will be copied to saved_p in ctqmc_make_evolve() subroutine
     integer, public, save, allocatable  :: is_cp(:,:)

! number of columns to be copied, in order to save copy time
     integer, public, save, allocatable  :: nc_cp(:,:)

! saved parts of matrices product, for previous accepted configuration
     real(dp), public, save, allocatable :: saved_p(:,:,:,:)

! saved parts of matrices product, for new proposed configuration
     real(dp), public, save, allocatable :: saved_n(:,:,:,:)

!!========================================================================
!!>>> declare private variables                                        <<<
!!========================================================================

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: ctqmc_allocate_memory_part
     public :: ctqmc_deallocate_memory_part

     public :: cat_make_npart
     public :: cat_make_trace

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!>>> ctqmc_allocate_memory_part: allocate memory for part related variables
  subroutine ctqmc_allocate_memory_part()
     implicit none

! allocate memory
     allocate(nop(npart),         stat=istat)
     allocate(ops(npart),         stat=istat)
     allocate(ope(npart),         stat=istat)

     allocate(renew(npart),       stat=istat)
     allocate(async(npart,nsect), stat=istat)
     allocate(is_cp(npart,nsect), stat=istat)
     allocate(nc_cp(npart,nsect), stat=istat)

     allocate(saved_p(max_dim_sect,max_dim_sect,npart,nsect), stat=istat)
     allocate(saved_n(max_dim_sect,max_dim_sect,npart,nsect), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_part','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     nop   = 0
     ops   = 0
     ope   = 0

     renew = 0
     async = 0
     is_cp = 0
     nc_cp = 0

     saved_p = zero
     saved_n = zero

     return
  end subroutine ctqmc_allocate_memory_part

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!>>> ctqmc_deallocate_memory_part: deallocate memory for part related variables
  subroutine ctqmc_deallocate_memory_part()
     implicit none

     if ( allocated(nop)     ) deallocate(nop    )
     if ( allocated(ops)     ) deallocate(ops    )
     if ( allocated(ope)     ) deallocate(ope    )

     if ( allocated(renew)   ) deallocate(renew  )
     if ( allocated(async)   ) deallocate(async  )
     if ( allocated(is_cp)   ) deallocate(is_cp  )
     if ( allocated(nc_cp)   ) deallocate(nc_cp  )

     if ( allocated(saved_p) ) deallocate(saved_p)
     if ( allocated(saved_n) ) deallocate(saved_n)

     return
  end subroutine ctqmc_deallocate_memory_part

!!========================================================================
!!>>> core service subroutines                                         <<<
!!========================================================================

!!>>> cat_make_npart: it is used to determine renew, which parts should
!!>>> be recalculated, is_cp is also reseted in this subroutine
  subroutine cat_make_npart(cmode, csize, index_loc, tau_s, tau_e)
     implicit none

! external arguments
! mode for different Monte Carlo moves
     integer, intent(in)  :: cmode

! total number of operators for current diagram
     integer, intent(in)  :: csize

! local version of index_t
     integer, intent(in)  :: index_loc(mkink)

! imaginary time value of operator A, only valid in cmode = 1 or 2
     real(dp), intent(in) :: tau_s

! imaginary time value of operator B, only valid in cmode = 1 or 2
     real(dp), intent(in) :: tau_e

! local variables
! loop index
     integer  :: i
     integer  :: j

! position of the operator A and operator B, index of part
     integer  :: tis
     integer  :: tie
     integer  :: tip

! length in imaginary time axis for each part
     real(dp) :: interval

! evaluate interval at first
     interval = beta / real(npart)

! init key arrays
     nop = 0
     ops = 0
     ope = 0

! init global arrays (renew and is_cp)
     renew = 0
     is_cp = 0

! calculate number of operators for each part
     do i=1,csize
         j = ceiling( time_v( index_loc(i) ) / interval )
         nop(j) = nop(j) + 1
     enddo ! over i={1,csize} loop

! calculate the start and end index of operators for each part
     do i=1,npart
         if ( nop(i) > 0 ) then
             ops(i) = 1
             do j=1,i-1
                 ops(i) = ops(i) + nop(j)
             enddo ! over j={1,i-1} loop
             ope(i) = ops(i) + nop(i) - 1
         endif ! back if ( nop(i) > 0 ) block
     enddo ! over i={1,npart} loop

! next we have to figure out which parts should be updated
! case 1: only some parts need to be updated
     if ( cmode == 1 .or. cmode == 2 ) then

! get the position of operator A and operator B
         tis = ceiling( tau_s / interval )
         tie = ceiling( tau_e / interval )

! determine the influence of operator A, which part should be recalculated
         renew(tis) = 1
! special attention: if operator A is on the left or right boundary, then
! the neighbour part should be recalculated as well
         if ( nop(tis) > 0 ) then
             if ( tau_s >= time_v( index_loc( ope(tis) ) ) ) then
                 tip = tis + 1
                 do while ( tip <= npart )
                     if ( nop(tip) > 0 ) then
                         renew(tip) = 1; EXIT
                     endif ! back if ( nop(tip) > 0 ) block
                     tip = tip + 1
                 enddo ! over do while loop
             endif ! back if ( tau_s >= time_v( index_t( ope(tis) ) ) ) block
         else
             tip = tis + 1
             do while ( tip <= npart )
                 if ( nop(tip) > 0 ) then
                     renew(tip) = 1; EXIT
                 endif ! back if ( nop(tip) > 0 ) block
                 tip = tip + 1
             enddo ! over do while loop
         endif ! back if ( nop(tis) > 0 ) block

! determine the influence of operator B, which part should be recalculated
         renew(tie) = 1
! special attention: if operator B is on the left or right boundary, then
! the neighbour part should be recalculated as well
         if ( nop(tie) > 0 ) then
             if ( tau_e >= time_v( index_loc( ope(tie) ) ) ) then
                 tip = tie + 1
                 do while ( tip <= npart )
                     if ( nop(tip) > 0 ) then
                         renew(tip) = 1; EXIT
                     endif ! back if ( nop(tip) > 0 ) block
                     tip = tip + 1
                 enddo ! over do while loop
             endif ! back if ( tau_e >= time_v( index_t( ope(tie) ) ) ) block
         else
             tip = tie + 1
             do while ( tip <= npart )
                 if ( nop(tip) > 0 ) then
                     renew(tip) = 1; EXIT
                 endif ! back if ( nop(tip) > 0 ) block
                 tip = tip + 1
             enddo ! over do while loop
         endif ! back if ( nop(tie) > 0 ) block

! case 2: all parts should be updated
     else
         renew = 1
     endif

     return
  end subroutine cat_make_npart

!!>>> cat_make_trace: calculate the contribution to final trace for
!!>>> a given string
  subroutine cat_make_trace(csize, string, index_loc, expt_loc, trace)
     implicit none

! external variables
! number of total fermion operators
     integer, intent(in)   :: csize

! evolution string for this sector
     integer, intent(in)   :: string(csize+1)

! memory address index of fermion operators
     integer, intent(in)   :: index_loc(mkink)

! diagonal elements of last time-evolution matrices
     real(dp), intent(in)  :: expt_loc(ncfgs)

! the calculated trace of this sector
     real(dp), intent(out) :: trace

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: l

! type for current operator
     integer  :: vt

! flavor channel for current operator
     integer  :: vf

! start index of this sector
     integer  :: indx

! dimension for the sectors
     integer  :: dim1
     integer  :: dim2
     integer  :: dim3
     integer  :: dim4

! index for sectors
     integer  :: isect
     integer  :: sect1
     integer  :: sect2

! the first part with non-zero fermion operators
     integer  :: fpart

! counter for fermion operators
     integer  :: counter

! real(dp) dummy matrices
     real(dp) :: mat_r(max_dim_sect,max_dim_sect)
     real(dp) :: mat_t(max_dim_sect,max_dim_sect)

! initialize dummy arrays
     mat_r = zero
     mat_t = zero

! select the first sector in the string
     isect = string(1)
     dim1  = sectors( string(1) )%ndim

! determine fpart
     fpart = 0
     do i=1,npart
         if ( nop(i) > 0 ) then
             fpart = i; EXIT
         endif ! back if ( nop(i) > 0 ) block
     enddo ! over i={1,npart} loop

! next we perform time evolution from left to right: 0 -> \beta
! loop over all the parts
     do i=1,npart

! empty part, we just skip it
         if ( nop(i) == 0 ) CYCLE

! this part should be recalcuated
         if ( renew(i) == 1 .or. async(i,isect) == 1 ) then
             sect1 = string(ope(i)+1)
             sect2 = string(ops(i))
             dim4 = sectors(sect2)%ndim
             saved_n(:,:,i,isect) = zero

! set its copy status
             is_cp(i,isect) = 1
             nc_cp(i,isect) = dim4

! loop over all the fermion operators in this part
             counter = 0
             do j=ops(i),ope(i)
                 counter = counter + 1
                 indx = sectors( string(j)   )%istart
                 dim2 = sectors( string(j+1) )%ndim
                 dim3 = sectors( string(j)   )%ndim

! multiply the diagonal matrix of time evolution operator
                 if ( counter > 1 ) then
                     do l=1,dim4
                         do k=1,dim3
                             mat_t(k,l) = saved_n(k,l,i,isect) * expt_v(indx+k-1,index_loc(j))
                         enddo ! over k={1,dim3} loop
                     enddo ! over l={1,dim4} loop
                 else
                     mat_t = zero
                     do k=1,dim3
                         mat_t(k,k) = expt_v(indx+k-1,index_loc(j))
                     enddo ! over k={1,dim3} loop
                 endif ! back if ( counter > 1 ) block

! multiply the matrix of fermion operator
                 vt = type_v( index_loc(j) )
                 vf = flvr_v( index_loc(j) )
                 call dgemm( 'N', 'N', dim2, dim4, dim3, &
                                                    one, &
                   sectors( string(j) )%fmat(vf,vt)%val, &
                                            dim2, mat_t, &
                                           max_dim_sect, &
                             zero, saved_n(:,:,i,isect), &
                                           max_dim_sect )
             enddo ! over j={ops(i),ope(i)} loop

! multiply this part with the rest parts
             if ( i > fpart ) then
                 call dgemm( 'N', 'N', dim2, dim1, dim4, &
                              one, saved_n(:,:,i,isect), &
                                           max_dim_sect, &
                                                  mat_r, &
                                           max_dim_sect, &
                                            zero, mat_t, &
                                           max_dim_sect )
                 mat_r(:,1:dim1) = mat_t(:,1:dim1)
             else
                 mat_r(:,1:dim1) = saved_n(:,1:dim1,i,isect)
             endif ! back if ( i > fpart ) block

! this part has been calculated previously, just use its results
         else
             sect1 = string(ope(i)+1)
             sect2 = string(ops(i))
             dim2 = sectors(sect1)%ndim
             dim3 = sectors(sect2)%ndim
             if ( i > fpart ) then
                 call dgemm( 'N', 'N', dim2, dim1, dim3, &
                              one, saved_p(:,:,i,isect), &
                                           max_dim_sect, &
                                                  mat_r, &
                                           max_dim_sect, &
                              zero, mat_t, max_dim_sect )
                 mat_r(:,1:dim1) = mat_t(:,1:dim1)
             else
                 mat_r(:,1:dim1) = saved_p(:,1:dim1,i,isect)
             endif ! back if ( i > fpart ) block

         endif ! back if ( renew(i) == 1 .or. async(i,isect) == 1 )  block

! setup the start sector for next part
         isect = sect1
     enddo ! over i={1,npart} loop

! special treatment of the last time evolution operator
     indx = sectors( string(1) )%istart

! no fermion operators
     if ( csize == 0 ) then
         do k=1,dim1
             mat_r(k,k) = expt_loc(indx+k-1)
         enddo ! over k={1,dim1} loop
! multiply the last time evolution operator
     else
         do l=1,dim1
             do k=1,dim1
                 mat_r(k,l) = mat_r(k,l) * expt_loc(indx+k-1)
             enddo ! over k={1,dim1} loop
         enddo ! over l={1,dim1} loop
     endif ! back if ( csize == 0 ) block

! calculate the trace and store the final product
     trace = zero
     do j=1,sectors( string(1) )%ndim
         trace = trace + mat_r(j,j)
         sectors( string(1) )%prod(j) = mat_r(j,j)
     enddo ! over j={1,sectors( string(1) )%ndim} loop

     return
  end subroutine cat_make_trace

  end module m_part
