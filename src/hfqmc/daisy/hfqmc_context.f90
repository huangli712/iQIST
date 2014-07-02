!-------------------------------------------------------------------------
! project : daisy
! program : hfqmc_core module
!           hfqmc_umat module
!           hfqmc_base module
!           context    module
! source  : hfqmc_context.f90
! type    : module
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 10/24/2008 by li huang
!           10/27/2008 by li huang
!           12/18/2008 by li huang
!           12/30/2008 by li huang
!           01/03/2009 by li huang
!           03/16/2009 by li huang
!           04/18/2009 by li huang
!           08/10/2009 by li huang
!           08/24/2009 by li huang
!           12/24/2009 by li huang
!           02/26/2010 by li huang
!           03/08/2010 by li huang
!           03/27/2010 by li huang
! purpose : define the key data structure and global arrays/variables for
!           Hirsch-Fye quantum Monte Carlo (HFQMC) quantum impurity solver
!           and dynamical mean field theory (DMFT) self-consistent engine
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!=========================================================================
!>>> module hfqmc_core                                                 <<<
!=========================================================================
!>>> containing core (internal) variables used by Hirsch-Fye quantum Monte
! Carlo quantum impurity solver
  module hfqmc_core
     use constants, only : dp

     implicit none

! cyclic steps of different orbitals for delayed update algorithm
     integer, public, save, allocatable  :: ktep(:)

! auxliary orbital index matrix
     integer, public, save, allocatable  :: pmat(:,:)

! reduced Coulomb interaction vector
     real(dp), public, save, allocatable :: umat(:)

! $\lambda$, Hubbard-Stratonovich transformation parameter
     real(dp), public, save, allocatable :: lmat(:)

! auxiliary ising-like fields
     real(dp), public, save, allocatable :: imat(:,:)

! $\sigma$, Pauli matrix
     real(dp), public, save, allocatable :: smat(:,:)

! diagonal elements of green's function matrix, used by delayed update algorithm
     real(dp), public, save, allocatable :: diag(:,:)

! column elements of green's function matrix, used by delayed update algorithm
     real(dp), public, save, allocatable :: atep(:,:,:)

! row elements of green's function matrix, used by delayed update algorithm
     real(dp), public, save, allocatable :: btep(:,:,:)

! green's function matrix
     real(dp), public, save, allocatable :: gmat(:,:,:)

! weiss's function matrix
     real(dp), public, save, allocatable :: wmat(:,:,:)

  end module hfqmc_core

!=========================================================================
!>>> module hfqmc_umat                                                 <<<
!=========================================================================
!>>> containing util-matrix related arrays used by Hirsch-Fye quantum Monte
! Carlo quantum impurity solver
  module hfqmc_umat
     use constants, only : dp

     implicit none

! symmetry properties for correlated orbitals
     integer, public, save, allocatable  :: symm(:)

! impurity level for correlated orbitals
     real(dp), public, save, allocatable :: eimp(:)

! quasiparticle weight, Z
     real(dp), public, save, allocatable :: quas(:)

! impurity occupation number, < n_i >
     real(dp), public, save, allocatable :: nmat(:)

! imaginary time mesh
     real(dp), public, save, allocatable :: tmesh(:)

! real matsubara frequency mesh
     real(dp), public, save, allocatable :: rmesh(:)

! impurity double occupation number matrix, < n_i n_j >
     real(dp), public, save, allocatable :: nnmat(:,:)

! identity matrix
     real(dp), public, save, allocatable :: unity(:,:)

  end module hfqmc_umat

!=========================================================================
!>>> module hfqmc_base                                                 <<<
!=========================================================================
!>>> containing basic arrays used by the dynamical mean field theory self-
! consistent engine
  module hfqmc_base
     use constants, only : dp

     implicit none

! impurity green's function in imaginary time axis
     real(dp), public, save, allocatable    :: gtau(:,:)

! bath weiss's function in imaginary time axis
     real(dp), public, save, allocatable    :: wtau(:,:)

! impurity green's function in matsubara space
     complex(dp), public, save, allocatable :: grnf(:,:)

! bath weiss's function in matsubara space
     complex(dp), public, save, allocatable :: wssf(:,:)

! self-energy function in matsubara space
     complex(dp), public, save, allocatable :: sig1(:,:)

! self-energy function in matsubara space
     complex(dp), public, save, allocatable :: sig2(:,:)

  end module hfqmc_base

!=========================================================================
!>>> module context                                                    <<<
!=========================================================================
!>>> containing memory management subroutines and define global variables
  module context
     use constants
     use control

     use hfqmc_core
     use hfqmc_umat
     use hfqmc_base

     implicit none

! status flag
     integer, private :: istat

! declaration of module procedures: allocate memory
     public :: hfqmc_allocate_memory_core
     public :: hfqmc_allocate_memory_umat
     public :: hfqmc_allocate_memory_base

! declaration of module procedures: deallocate memory
     public :: hfqmc_deallocate_memory_core
     public :: hfqmc_deallocate_memory_umat
     public :: hfqmc_deallocate_memory_base

     contains

!=========================================================================
!>>> allocate memory subroutines                                       <<<
!=========================================================================

!>>> allocate memory for core-related variables
     subroutine hfqmc_allocate_memory_core()
         implicit none

! allocate memory
         allocate(ktep(norbs),             stat=istat)
         allocate(pmat(nsing,nspin),       stat=istat)

         allocate(umat(nsing),             stat=istat)
         allocate(lmat(nsing),             stat=istat)

         allocate(imat(ntime,nsing),       stat=istat)
         allocate(smat(norbs,nsing),       stat=istat)

         allocate(diag(ntime,norbs),       stat=istat)
         allocate(atep(ntime,mstep,norbs), stat=istat)
         allocate(btep(ntime,mstep,norbs), stat=istat)

         allocate(gmat(ntime,ntime,norbs), stat=istat)
         allocate(wmat(ntime,ntime,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call hfqmc_print_error('hfqmc_allocate_memory_core','can not allocate enough memory')
         endif

! initialize them
         ktep = 0
         pmat = 0

         umat = zero
         lmat = zero

         imat = zero
         smat = zero

         diag = zero
         atep = zero
         btep = zero

         gmat = zero
         wmat = zero

         return
     end subroutine hfqmc_allocate_memory_core

!>>> allocate memory for umat-related variables
     subroutine hfqmc_allocate_memory_umat()
         implicit none

! allocate memory
         allocate(symm(norbs),        stat=istat)

         allocate(eimp(norbs),        stat=istat)

         allocate(quas(norbs),        stat=istat)
         allocate(nmat(norbs),        stat=istat)

         allocate(tmesh(ntime),       stat=istat)
         allocate(rmesh(mfreq),       stat=istat)

         allocate(nnmat(norbs,norbs), stat=istat)

         allocate(unity(ntime,ntime), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call hfqmc_print_error('hfqmc_allocate_memory_umat','can not allocate enough memory')
         endif

! initialize them
         symm  = 0

         eimp  = zero

         quas  = zero
         nmat  = zero

         tmesh = zero
         rmesh = zero

         nnmat = zero

         unity = zero

         return
     end subroutine hfqmc_allocate_memory_umat

!>>> allocate memory for base-related variables
     subroutine hfqmc_allocate_memory_base()
         implicit none

! allocate memory
         allocate(gtau(ntime,norbs), stat=istat)
         allocate(wtau(ntime,norbs), stat=istat)

         allocate(grnf(mfreq,norbs), stat=istat)
         allocate(wssf(mfreq,norbs), stat=istat)

         allocate(sig1(mfreq,norbs), stat=istat)
         allocate(sig2(mfreq,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call hfqmc_print_error('hfqmc_allocate_memory_base','can not allocate enough memory')
         endif

! initialize them
         gtau = zero
         wtau = zero

         grnf = czero
         wssf = czero

         sig1 = czero
         sig2 = czero

         return
     end subroutine hfqmc_allocate_memory_base

!=========================================================================
!>>> deallocate memory subroutines                                     <<<
!=========================================================================

!>>> deallocate memory for core-related variables
     subroutine hfqmc_deallocate_memory_core()
         implicit none

         if ( allocated(ktep) )    deallocate(ktep)
         if ( allocated(pmat) )    deallocate(pmat)

         if ( allocated(umat) )    deallocate(umat)
         if ( allocated(lmat) )    deallocate(lmat)

         if ( allocated(imat) )    deallocate(imat)
         if ( allocated(smat) )    deallocate(smat)

         if ( allocated(diag) )    deallocate(diag)
         if ( allocated(atep) )    deallocate(atep)
         if ( allocated(btep) )    deallocate(btep)

         if ( allocated(gmat) )    deallocate(gmat)
         if ( allocated(wmat) )    deallocate(wmat)

         return
     end subroutine hfqmc_deallocate_memory_core

!>>> deallocate memory for umat-related variables
     subroutine hfqmc_deallocate_memory_umat()
         implicit none

         if ( allocated(symm)  )   deallocate(symm )

         if ( allocated(eimp)  )   deallocate(eimp )

         if ( allocated(quas)  )   deallocate(quas )
         if ( allocated(nmat)  )   deallocate(nmat )

         if ( allocated(tmesh) )   deallocate(tmesh)
         if ( allocated(rmesh) )   deallocate(rmesh)

         if ( allocated(nnmat) )   deallocate(nnmat)

         if ( allocated(unity) )   deallocate(unity)

         return
     end subroutine hfqmc_deallocate_memory_umat

!>>> deallocate memory for base-related variables
     subroutine hfqmc_deallocate_memory_base()
         implicit none

         if ( allocated(gtau) )    deallocate(gtau)
         if ( allocated(wtau) )    deallocate(wtau)

         if ( allocated(grnf) )    deallocate(grnf)
         if ( allocated(wssf) )    deallocate(wssf)

         if ( allocated(sig1) )    deallocate(sig1)
         if ( allocated(sig2) )    deallocate(sig2)

         return
     end subroutine hfqmc_deallocate_memory_base

  end module context
