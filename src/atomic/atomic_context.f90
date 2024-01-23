!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : m_fock   module
!!!           m_sector module
!!!           m_spmat  module
!!! source  : atomic_context.f90
!!! type    : modules
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           01/23/2024 by li huang (last modified)
!!! purpose : define global data structures, arrays, and variables for
!!!           the atomic eigenvalue problem solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> module m_fock                                                    <<<
!!========================================================================

!!
!! @mod m_fock
!!
!! define the Fock basis, density matrix, spin matrix, atomic Hamiltonian
!! and the corresponding eigensystem
!!
  module m_fock
     use constants, only : dp
     use constants, only : zero, czero

     use control, only : norbs, ncfgs

     implicit none

!!
!! @var dim_sub_n
!!
!! dimension of subspace with total electron N
!! if i is the number of electrons, then dim_sub_n(i) will tell you
!! how many Fock states there are in the subspace with i electrons
!!
     integer, public, save, allocatable  :: dim_sub_n(:)

!!
!! @var bin_basis
!!
!! binary form of Fock state (something like |110110>)
!! bin_basis(:,i) denotes the binary form of the i-th Fock state
!!
     integer, public, save, allocatable  :: bin_basis(:,:)

!!
!! @var dec_basis
!!
!! decimal form of Fock state (a decimal number)
!! dec_basis(i) denotes the decimal form of the i-th Fock state
!!
     integer, public, save, allocatable  :: dec_basis(:)

!!
!! @var ind_basis
!!
!! index of Fock state, given their decimal number
!! for a given decimal number, ind_basis(i) will tell you what
!! the corresponding Fock state index is
!!
     integer, public, save, allocatable  :: ind_basis(:)

!!
!! @var eval
!!
!! eigenvalues of atomic Hamiltonian
!!
     real(dp), public, save, allocatable :: eval(:)

!!
!! @var evec
!!
!! eigenvectors of atomic Hamiltonian
!!
     real(dp), public, save, allocatable :: evec(:,:)

!!
!! @var occu
!!
!! density matrix (occupancy) for the atomic eigenstates
!!
     real(dp), public, save, allocatable :: occu(:,:)

!!
!! @var spin
!!
!! magnetic moment (Sz) for the atomic eigenstates
!!
     real(dp), public, save, allocatable :: spin(:,:)

!!
!! @var fmat
!!
!! annihilation operator matrix, < alpha | f | beta >
!! where | alpha > and | beta > are the atomic eigenstates 
!!
     real(dp), public, save, allocatable :: fmat(:,:,:)

!!
!! @var hmat
!!
!! atomic Hamiltonian in the Fock basis
!!
     complex(dp), public, save, allocatable :: hmat(:,:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     ! declaration of module procedures: allocate memory
     public :: cat_alloc_fock_basis
     public :: cat_alloc_fock_eigen

     ! declaration of module procedures: deallocate memory
     public :: cat_free_fock_basis
     public :: cat_free_fock_eigen

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!
!! @sub cat_alloc_fock_basis
!!
!! allocate memory for the Fock basis
!!
  subroutine cat_alloc_fock_basis()
     implicit none

!! local variables
     ! the status flag
     integer :: istat

!! [body

     ! allocate memory
     allocate(dim_sub_n(0:norbs),     stat=istat)
     allocate(bin_basis(norbs,ncfgs), stat=istat)
     allocate(dec_basis(ncfgs),       stat=istat)
     allocate(ind_basis(0:ncfgs-1),   stat=istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_fock_basis','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     dim_sub_n = 0
     bin_basis = 0
     dec_basis = 0
     ind_basis = 0

!! body]

     return
  end subroutine cat_alloc_fock_basis

!!
!! @sub cat_alloc_fock_eigen
!!
!! allocate memory for atomic eigensystem that
!! is initialized in the Fock basis
!!
  subroutine cat_alloc_fock_eigen()
     implicit none

!! local variables
     ! the status flag
     integer :: istat

!! [body

     ! allocate memory
     allocate(eval(ncfgs),             stat=istat)
     allocate(evec(ncfgs,ncfgs),       stat=istat)
     allocate(occu(ncfgs,ncfgs),       stat=istat)
     allocate(spin(ncfgs,ncfgs),       stat=istat)
     allocate(fmat(ncfgs,ncfgs,norbs), stat=istat)
     !
     allocate(hmat(ncfgs,ncfgs),       stat=istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_fock_eigen','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     eval = zero
     evec = zero
     occu = zero
     spin = zero
     fmat = zero

     hmat = czero

!! body]

     return
  end subroutine cat_alloc_fock_eigen

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!
!! @sub cat_free_fock_basis
!!
!! deallocate memory for the Fock basis
!!
  subroutine cat_free_fock_basis()
     implicit none

!! [body

     if ( allocated(dim_sub_n) ) deallocate(dim_sub_n)
     if ( allocated(bin_basis) ) deallocate(bin_basis)
     if ( allocated(dec_basis) ) deallocate(dec_basis)
     if ( allocated(ind_basis) ) deallocate(ind_basis)

!! body]

     return
  end subroutine cat_free_fock_basis

!!
!! @sub cat_free_fock_eigen
!!
!! deallocate memory for atomic eigensystem that
!! is initialized in the Fock basis
!!
  subroutine cat_free_fock_eigen()
     implicit none

!! [body

     if ( allocated(eval) ) deallocate(eval)
     if ( allocated(evec) ) deallocate(evec)
     if ( allocated(occu) ) deallocate(occu)
     if ( allocated(spin) ) deallocate(spin)
     if ( allocated(fmat) ) deallocate(fmat)

     if ( allocated(hmat) ) deallocate(hmat)

!! body]

     return
  end subroutine cat_free_fock_eigen

  end module m_fock

!!========================================================================
!!>>> module m_sector                                                  <<<
!!========================================================================

!!
!! @mod m_sector
!!
!! data structure for good quantum numbers (GQNs) algorithm
!!
  module m_sector
     use constants, only : dp
     use constants, only : zero, czero

     implicit none

!!
!! @struct t_fmat
!!
!! data structure for one F-matrix
!!
     private :: t_fmat
     type t_fmat

         ! the dimension, n x m
         integer :: n
         integer :: m

         ! the memory space for the matrix
         real(dp), allocatable :: val(:,:)

     end type t_fmat

!!
!! @struct t_sector
!!
!! data structure for one sector
!!
     public :: t_sector
     type t_sector

         ! the dimension of this sector
         integer :: ndim
 
         ! number of fermion operators
         integer :: nops

         ! the start index of this sector
         integer :: istart

         ! total number of electrons N
         integer :: nele

         ! z component of spin: Sz
         integer :: sz
 
         ! z component of spin-orbit momentum: Jz
         integer :: jz

         ! PS good quantum number
         integer :: ps

         ! the Fock basis index of this sector
         integer, allocatable  :: basis(:)

         ! the next sector after a fermion operator acts on this sector
         !
         ! next(nops,0) for annihilation
         ! next(nops,1) for creation operators
         !
         ! if it is -1, it means the next sector is null (outside of the
         ! Hilbert space). otherwise, it is the index of next sector
         integer, allocatable  :: next(:,:)

         ! the eigenvalues
         real(dp), allocatable :: eval(:)

         ! the eigenvectors
         ! since Hamiltonian must be real, then it is real as well
         real(dp), allocatable :: evec(:,:)

         ! the Hamiltonian of this sector
         complex(dp), allocatable :: hmat(:,:)

         ! the F-matrix between this sector and all other sectors
         !
         ! fmat(nops,0) for annihilation
         ! fmat(nops,1) for creation operators
         !
         ! if this sector doesn't point to some other sectors, then the
         ! pointer is null
         type (t_fmat), allocatable :: fmat(:,:)

     end type t_sector

!!
!! @var nsectors
!!
!! number of sectors
!!
     integer, public, save  :: nsectors

!!
!! @var max_dim_sect
!!
!! maximum dimension of sectors
!!
     integer, public, save  :: max_dim_sect

!!
!! @var ave_dim_sect
!!
!! average dimension of sectors
!!
     real(dp), public, save :: ave_dim_sect

!!
!! @var sectors
!!
!! all the sectors
!!
     type (t_sector), public, save, allocatable :: sectors(:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     ! declaration of module procedures: allocate memory
     public :: cat_alloc_fmat
     public :: cat_alloc_sector
     public :: cat_alloc_sectors

     ! declaration of module procedures: deallocate memory
     public :: cat_free_fmat
     public :: cat_free_sector
     public :: cat_free_sectors

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!
!! @sub cat_alloc_fmat
!!
!! allocate memory for one fmat
!!
  subroutine cat_alloc_fmat(one_fmat)
     implicit none

!! external arguments
     ! the fmat
     type (t_fmat), intent(inout) :: one_fmat

!! local variables
     ! the status flag
     integer :: istat

!! [body

     ! allocate memory
     allocate(one_fmat%val(one_fmat%n,one_fmat%m), stat=istat)

     ! check status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_fmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize it
     one_fmat%val = zero

!! body]

     return
  end subroutine cat_alloc_fmat

!!
!! @sub cat_alloc_sector
!!
!! allocate memory for one sector
!!
  subroutine cat_alloc_sector(one_sector)
     implicit none

!! external arguments
     ! the sector
     type (t_sector), intent(inout) :: one_sector

!! local variables
     ! loop index
     integer :: i
     integer :: j

     ! the status flag
     integer :: istat

!! [body

     ! allocate memory
     allocate(one_sector%basis(one_sector%ndim),                stat=istat)
     allocate(one_sector%next(one_sector%nops,0:1),             stat=istat)
     allocate(one_sector%eval(one_sector%ndim),                 stat=istat)
     allocate(one_sector%evec(one_sector%ndim,one_sector%ndim), stat=istat)
     allocate(one_sector%hmat(one_sector%ndim,one_sector%ndim), stat=istat)
     allocate(one_sector%fmat(one_sector%nops,0:1),             stat=istat)

     ! check status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_sector','can not allocate enough memory')
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
        enddo ! over j={0,1} loop
     enddo ! over i={1,one_sector%nops} loop

!! body]

     return
  end subroutine cat_alloc_sector

!!
!! @sub cat_alloc_sectors
!!
!! allocate memory for sectors
!!
  subroutine cat_alloc_sectors()
     implicit none

!! local variables
     ! the status flag
     integer :: istat

!! [body

     ! allocate memory
     allocate(sectors(nsectors), stat=istat)

     ! check status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_sectors','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

!! body]

     return
  end subroutine cat_alloc_sectors

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!
!! @sub cat_free_fmat
!!
!! deallocate memory for one fmat
!!
  subroutine cat_free_fmat(one_fmat)
     implicit none

!! external arguments
     ! the fmat
     type (t_fmat), intent(inout) :: one_fmat

!! [body

     if ( allocated(one_fmat%val) ) deallocate(one_fmat%val)

!! body]

     return
  end subroutine cat_free_fmat

!!
!! @sub cat_free_sector
!!
!! deallocate memory for one sector
!!
  subroutine cat_free_sector(one_sector)
     implicit none

!! external arguments
     ! the sector
     type (t_sector), intent(inout) :: one_sector

!! local variables
     ! loop index
     integer :: i
     integer :: j

!! [body

     if ( allocated(one_sector%basis) ) deallocate(one_sector%basis)
     if ( allocated(one_sector%next)  ) deallocate(one_sector%next )
     if ( allocated(one_sector%eval)  ) deallocate(one_sector%eval )
     if ( allocated(one_sector%evec)  ) deallocate(one_sector%evec )
     if ( allocated(one_sector%hmat)  ) deallocate(one_sector%hmat )

     ! deallocate fmat one by one
     if ( allocated(one_sector%fmat)  ) then
         do i=1,one_sector%nops
             do j=0,1
                 call cat_free_fmat(one_sector%fmat(i,j))
             enddo ! over j={0,1} loop
         enddo ! over i={1,one_sector%nops} loop
         deallocate(one_sector%fmat)
     endif ! back if ( allocated(one_sector%fmat)  ) block

!! body]

     return
  end subroutine cat_free_sector

!!
!! @sub cat_free_sectors
!!
!! deallocate memory of sectors
!!
  subroutine cat_free_sectors()
     implicit none

!! local variables
     ! loop index
     integer :: i

!! [body

     ! deallocate memory for arrays in T_sector
     ! before deallocating sectors to avoid memory leak
     if ( allocated(sectors) ) then
         do i=1,nsectors
             call cat_free_sector(sectors(i))
         enddo ! over i={1,nsectors} loop
         deallocate(sectors)
     endif ! back if ( allocated(sectors) ) block

!! body]

     return
  end subroutine cat_free_sectors

  end module m_sector

!!========================================================================
!!>>> module m_spmat                                                   <<<
!!========================================================================

!!
!! @mod m_spmat
!!
!! it contains single particle related matrices, including
!!     crystal field,
!!     spin-orbital coupling,
!!     Coulomb interaction U tensor, etc
!!
  module m_spmat
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs

     implicit none

!!
!! @var umat
!!
!! Coulomb interaction U tensor
!!
     complex(dp), public, allocatable, save :: umat(:,:,:,:)

!!
!! @var cmat
!!
!! crystal field (CF)
!!
     complex(dp), public, allocatable, save :: cmat(:,:)

!!
!! @var smat
!!
!! spin-orbital coupling (SOC)
!!
     complex(dp), public, allocatable, save :: smat(:,:)

!!
!! @var emat
!!
!! onsite energy (CF + SOC) of impurity
!!
     complex(dp), public, allocatable, save :: emat(:,:)

!!
!! @var tmat
!!
!! transformation matrix from original basis to natural basis
!!
     complex(dp), public, allocatable, save :: tmat(:,:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     ! declaration of module procedures: allocate memory
     public :: cat_alloc_spmat

     ! declaration of module procedures: deallocate memory
     public :: cat_free_spmat

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!
!! @sub cat_alloc_spmat
!!
!! allocate memory for single particle matrices
!!
  subroutine cat_alloc_spmat()
     implicit none

!! local variables
     ! the status flag
     integer :: istat

!! [body

     ! allocate memory
     allocate(umat(norbs,norbs,norbs,norbs), stat=istat)
     allocate(cmat(norbs,norbs),             stat=istat)
     allocate(smat(norbs,norbs),             stat=istat)
     allocate(emat(norbs,norbs),             stat=istat)
     allocate(tmat(norbs,norbs),             stat=istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_spmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     umat = czero
     cmat = czero
     smat = czero
     emat = czero
     tmat = czero

!! body]

     return
  end subroutine cat_alloc_spmat

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!
!! @sub cat_free_spmat
!!
!! deallocate memory for single particle matrices
!!
  subroutine cat_free_spmat()
     implicit none

!! [body

     if ( allocated(umat) ) deallocate(umat)
     if ( allocated(cmat) ) deallocate(cmat)
     if ( allocated(smat) ) deallocate(smat)
     if ( allocated(emat) ) deallocate(emat)
     if ( allocated(tmat) ) deallocate(tmat)

!! body]

     return
  end subroutine cat_free_spmat

  end module m_spmat
