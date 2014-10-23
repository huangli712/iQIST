!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : m_full     module
!!!           m_sector   module
!!!           m_spmat    module
!!! source  : atomic_context.f90
!!! type    : modules
!!! authors : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!!           10/20/2014 by li huang
!!! purpose : define global data structures for the atomic eigenvalue
!!!           problem solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module m_full                                                    <<<
!!========================================================================

!!>>> define Fock basis of full Hilbert space and corresponding eigensystem
  module m_full
     use constants, only : dp, zero, czero

     use control, only : norbs, ncfgs

     implicit none

! dimension of subspace of total electron N
     integer, public, allocatable, save  :: dim_sub_n(:)

! binary form of Fock basis
     integer, public, allocatable, save  :: bin_basis(:,:)

! decimal form of Fock basis
     integer, public, allocatable, save  :: dec_basis(:)

! index of Fock basis, given their decimal number
     integer, public, allocatable, save  :: index_basis(:)

! eigenvalues of hmat
     real(dp), public, allocatable, save :: eval(:)

! eigenvectors of hmat
     real(dp), public, allocatable, save :: evec(:,:)

! F-matrix for annihilation fermion operators
     real(dp), public, allocatable, save :: fmat(:,:,:)

! occupany number for atomic eigenstates
     real(dp), public, allocatable, save :: occu(:,:)

! atomic Hamiltonian
     complex(dp), public, allocatable, save :: hmat(:,:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: alloc_m_full_basis
     public :: dealloc_m_full_basis

     public :: alloc_m_full
     public :: dealloc_m_full

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!>>> alloc_m_full_basis: allocate memory for Fock basis matrices
  subroutine alloc_m_full_basis()
     implicit none

! the status flag
     integer :: istat

! allocate memory
     allocate(dim_sub_n(0:norbs),     stat=istat)
     allocate(bin_basis(norbs,ncfgs), stat=istat)
     allocate(dec_basis(ncfgs),       stat=istat)
     allocate(index_basis(0:ncfgs-1), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('alloc_m_full_basis','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     dim_sub_n = 0
     bin_basis = 0
     dec_basis = 0
     index_basis = 0

     return
  end subroutine alloc_m_full_basis

!!>>> alloc_m_full: allocate memory for eigensystem defined in Fock basis
  subroutine alloc_m_full()
     implicit none

! the status flag
     integer :: istat

! allocate memory
     allocate(eval(ncfgs),             stat=istat)
     allocate(evec(ncfgs,ncfgs),       stat=istat)
     allocate(fmat(ncfgs,ncfgs,norbs), stat=istat)
     allocate(occu(ncfgs,ncfgs),       stat=istat)

     allocate(hmat(ncfgs,ncfgs),       stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('alloc_m_full','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     eval = zero
     evec = zero
     fmat = zero
     occu = zero

     hmat = czero

     return
  end subroutine alloc_m_full

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!>>> dealloc_m_full_basis: deallocate memory for these matrices
  subroutine dealloc_m_full_basis()
     implicit none

     if ( allocated(dim_sub_n)   ) deallocate(dim_sub_n  )
     if ( allocated(bin_basis)   ) deallocate(bin_basis  )
     if ( allocated(dec_basis)   ) deallocate(dec_basis  )
     if ( allocated(index_basis) ) deallocate(index_basis)

     return
  end subroutine dealloc_m_full_basis

!!>>> dealloc_m_full: deallocate memory for these matrices
  subroutine dealloc_m_full()
     implicit none

     if ( allocated(eval) ) deallocate(eval)
     if ( allocated(evec) ) deallocate(evec)
     if ( allocated(fmat) ) deallocate(fmat)
     if ( allocated(occu) ) deallocate(occu)
     if ( allocated(hmat) ) deallocate(hmat)

     return
  end subroutine dealloc_m_full

  end module m_full

!!========================================================================
!!>>> module m_sector                                                  <<<
!!========================================================================

!!>>> data structure for good quantum numbers (GQNs) algorithm
  module m_sector
     use constants, only : dp, zero, czero

     implicit none

! the F-matrix between any two sectors, it is just a matrix
     private :: T_fmat
     type T_fmat

! the dimension, n X m
         integer :: n
         integer :: m

! the items of the matrix
         real(dp), pointer :: item(:,:)

     end type T_fmat

! data structure for one sector
     public :: T_sector
     type :: T_sector

! the dimension of this sector
         integer :: ndim

! total number of electrons N
         integer :: nele

! number of fermion operators
         integer :: nops

! the start index of this sector
         integer :: istart

! the Fock basis index of this sector
         integer, pointer  :: basis(:)

! the next sector it points to when a fermion operator acts on this sector
! -1: outside of the Hilbert space, otherwise, it is the index of next sector
! next(nops,0:1), 0 for annihilation and 1 for creation operators, respectively
         integer, pointer  :: next(:,:)

! the eigenvalues
         real(dp), pointer :: eval(:)

! the eigenvectors, since Hamiltonian must be real, then it is real as well
         real(dp), pointer :: evec(:,:)

! the Hamiltonian of this sector
         complex(dp), pointer :: hmat(:,:)

! the F-matrix between this sector and all other sectors
! if this sector doesn't point to some other sectors, the pointer is null
! fmat(nops,0:1), 0 for annihilation and 1 for creation operators, respectively
         type (T_fmat), pointer :: fmat(:,:)

     end type T_sector

! number of sectors
     integer, public, save  :: nsectors

! maximum dimension of sectors
     integer, public, save  :: max_dim_sect

! average dimension of sectors
     real(dp), public, save :: ave_dim_sect

! all the sectors
     type(T_sector), public, allocatable, save :: sectors(:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: alloc_one_fmat
     public :: nullify_one_fmat
     public :: dealloc_one_fmat

     public :: alloc_one_sector
     public :: nullify_one_sector
     public :: dealloc_one_sector

     public :: alloc_m_sector
     public :: dealloc_m_sector

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!>>> alloc_one_fmat: allocate one fmat
  subroutine alloc_one_fmat(one_fmat)
     implicit none

! external arguments
! the fmat
     type(T_fmat), intent(inout) :: one_fmat

! local variables
! the status flag
     integer :: istat

! allocate memory
     allocate(one_fmat%item(one_fmat%n,one_fmat%m), stat=istat)

! check status
     if ( istat /= 0 ) then
         call s_print_error('alloc_one_fmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize it
     one_fmat%item = zero

     return
  end subroutine alloc_one_fmat

!>>> alloc_one_sector: allocate memory for one sector
  subroutine alloc_one_sector(one_sector)
     implicit none

! external arguments
! the sector
     type(T_sector), intent(inout) :: one_sector

! local variables
! the status flag
     integer :: istat

! loop index
     integer :: i
     integer :: j

! allocate memory
     allocate(one_sector%basis(one_sector%ndim),                stat=istat)
     allocate(one_sector%next(one_sector%nops,0:1),             stat=istat)
     allocate(one_sector%eval(one_sector%ndim),                 stat=istat)
     allocate(one_sector%evec(one_sector%ndim,one_sector%ndim), stat=istat)
     allocate(one_sector%hmat(one_sector%ndim,one_sector%ndim), stat=istat)
     allocate(one_sector%fmat(one_sector%nops,0:1),             stat=istat)

! check status
     if ( istat /= 0 ) then
         call s_print_error('alloc_one_sector','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     one_sector%basis = 0
     one_sector%next  = 0
     one_sector%eval  = zero
     one_sector%evec  = zero
     one_sector%hmat  = czero

! initialize fmat one by one
     do i=1,one_sector%nops
        do j=0,1
            one_sector%fmat(i,j)%n = 0
            one_sector%fmat(i,j)%m = 0
            call nullify_one_fmat(one_sector%fmat(i,j))
        enddo ! over j={0,1} loop
     enddo ! over i={1,one_sector%nops} loop

     return
  end subroutine alloc_one_sector

!!>>> alloc_m_sector: allocate memory for sectors
  subroutine alloc_m_sector()
     implicit none

! local variables
! the status flag
     integer :: istat

! loop index
     integer :: i

! allocate memory
     allocate(sectors(nsectors), stat=istat)

! check status
     if ( istat /= 0 ) then
         call s_print_error('alloc_m_sector','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! nullify each sector one by one
     do i=1,nsectors
         call nullify_one_sector(sectors(i))
     enddo ! over i={1,nsectors} loop

     return
  end subroutine alloc_m_sector

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!>>> nullify_one_fmat: nullify one fmat
  subroutine nullify_one_fmat(one_fmat)
     implicit none

! external arguments
! the fmat
     type(T_fmat), intent(inout) :: one_fmat

     nullify(one_fmat%item)

     return
  end subroutine nullify_one_fmat

!>>> dealloc_one_fmat: deallocate one fmat
  subroutine dealloc_one_fmat(one_fmat)
     implicit none

! external arguments
! the fmat
     type(T_fmat), intent(inout) :: one_fmat

     if ( associated(one_fmat%item) ) deallocate(one_fmat%item)

     return
  end subroutine dealloc_one_fmat

!!>>> nullify_one_sector: nullify one sector
  subroutine nullify_one_sector(one_sector)
     implicit none

! external arguments
! the sector
     type(T_sector), intent(inout) :: one_sector

     nullify(one_sector%basis)
     nullify(one_sector%next )
     nullify(one_sector%eval )
     nullify(one_sector%evec )
     nullify(one_sector%hmat )
     nullify(one_sector%fmat )

     return
  end subroutine nullify_one_sector

!!>>> dealloc_one_sector: deallocate memory for one sector
  subroutine dealloc_one_sector(one_sector)
     implicit none

! external arguments
! the sector
     type(T_sector), intent(inout) :: one_sector

! local variables
! loop index
     integer :: i
     integer :: j

     if ( associated(one_sector%basis) ) deallocate(one_sector%basis)
     if ( associated(one_sector%next)  ) deallocate(one_sector%next )
     if ( associated(one_sector%eval)  ) deallocate(one_sector%eval )
     if ( associated(one_sector%evec)  ) deallocate(one_sector%evec )
     if ( associated(one_sector%hmat)  ) deallocate(one_sector%hmat )

! deallocate fmat one by one
     do i=1,one_sector%nops
         do j=0,1
             call dealloc_one_fmat(one_sector%fmat(i,j))
         enddo ! over j={0,1} loop
     enddo ! over i={1,one_sector%nops} loop

     return
  end subroutine dealloc_one_sector

!!>>> dealloc_m_sector: deallocate memory of sectors
  subroutine dealloc_m_sector()
     implicit none

! local variables
! loop index
     integer :: i

! deallocate memory for pointers in T_sector
! before deallocating sectors to avoid memory leak
     do i=1,nsectors
         call dealloc_one_sector(sectors(i))
     enddo ! over i={1,nsectors} loop

     if ( allocated(sectors) ) deallocate(sectors)

     return
  end subroutine dealloc_m_sector

  end module m_sector

!!========================================================================
!!>>> module m_spmat                                                   <<<
!!========================================================================

!!>>> single particle related matrices, including:
!!>>> crystal field, spin-orbital coupling, Coulomb interaction U tensor
  module m_spmat
     use constants, only : dp, czero

     use control, only : norbs

     implicit none

! Coulomb interaction U tensor
     complex(dp), public, allocatable, save :: umat(:,:,:,:)

! crystal field (CF)
     complex(dp), public, allocatable, save :: cmat(:,:)

! spin-orbital coupling (SOC)
     complex(dp), public, allocatable, save :: smat(:,:)

! on-site energy (CF + SOC) of impurity
     complex(dp), public, allocatable, save :: emat(:,:)

! the transformation matrix from origional basis to natural basis
     complex(dp), public, allocatable, save :: tmat(:,:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: alloc_m_spmat
     public :: dealloc_m_spmat

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!>>> alloc_m_spmat: allocate memory for these matrices
  subroutine alloc_m_spmat()
     implicit none

! local variables
! the status flag
     integer :: istat

! allocate memory
     allocate(umat(norbs,norbs,norbs,norbs), stat=istat)
     allocate(cmat(norbs,norbs),             stat=istat)
     allocate(smat(norbs,norbs),             stat=istat)
     allocate(emat(norbs,norbs),             stat=istat)
     allocate(tmat(norbs,norbs),             stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('alloc_m_spmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     umat = czero
     cmat = czero
     smat = czero
     emat = czero
     tmat = czero

     return
  end subroutine alloc_m_spmat

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!>>> dealloc_m_spmat: deallocate memory for these matrices
  subroutine dealloc_m_spmat()
     implicit none

     if ( allocated(umat) ) deallocate(umat)
     if ( allocated(cmat) ) deallocate(cmat)
     if ( allocated(smat) ) deallocate(smat)
     if ( allocated(emat) ) deallocate(emat)
     if ( allocated(tmat) ) deallocate(tmat)

     return
  end subroutine dealloc_m_spmat

  end module m_spmat
