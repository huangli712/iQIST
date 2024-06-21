!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : m_fock   module
!!!           m_sector module
!!!           m_spmat  module
!!! source  : atomic_context.f90
!!! type    : modules
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           06/11/2024 by li huang (last modified)
!!! purpose : define global data structures, arrays, and variables for
!!!           the atomic eigenvalue problem solver.
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
!! density matrix (occupancy) in the atomic eigenbasis
!!
     real(dp), public, save, allocatable :: occu(:,:)

!!
!! @var spin
!!
!! magnetic moment (Sz) in the atomic eigenbasis
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
     public :: cat_alloc_hmat_only

     ! declaration of module procedures: deallocate memory
     public :: cat_free_fock_basis
     public :: cat_free_fock_eigen
     public :: cat_free_hmat_only

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
         call s_print_error('cat_alloc_fock_basis', &
             & 'can not allocate enough memory')
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
!! allocate memory for atomic eigensystem
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
         call s_print_error('cat_alloc_fock_eigen', &
             & 'can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     eval = zero
     evec = zero
     occu = zero
     spin = zero
     fmat = zero
     !
     hmat = czero

!! body]

     return
  end subroutine cat_alloc_fock_eigen

!!
!! @sub cat_alloc_hmat_only
!!
!! allocate memory for atomic Hamiltonian only. this subroutine is for
!! the automatic partition algorithm.
!!
  subroutine cat_alloc_hmat_only()
     implicit none

!! local variables
     ! the status flag
     integer :: istat

!! [body

     ! allocate memory
     allocate(hmat(ncfgs,ncfgs),       stat=istat)

     ! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_hmat_only', &
             & 'can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     hmat = czero

!! body]

     return
  end subroutine cat_alloc_hmat_only

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
!! deallocate memory for atomic eigensystem
!!
  subroutine cat_free_fock_eigen()
     implicit none

!! [body

     if ( allocated(eval) ) deallocate(eval)
     if ( allocated(evec) ) deallocate(evec)
     if ( allocated(occu) ) deallocate(occu)
     if ( allocated(spin) ) deallocate(spin)
     if ( allocated(fmat) ) deallocate(fmat)
     !
     if ( allocated(hmat) ) deallocate(hmat)

!! body]

     return
  end subroutine cat_free_fock_eigen

!!
!! @sub cat_free_hmat_only
!!
!! deallocate memory for atomic Hamiltonian only. this subroutine is for
!! the automatic partition.
!!
  subroutine cat_free_hmat_only()
     implicit none

!! [body

     if ( allocated(eval) ) deallocate(eval)
     if ( allocated(evec) ) deallocate(evec)
     if ( allocated(occu) ) deallocate(occu)
     if ( allocated(spin) ) deallocate(spin)
     if ( allocated(fmat) ) deallocate(fmat)
     !
     if ( allocated(hmat) ) deallocate(hmat)

!! body]

     return
  end subroutine cat_free_hmat_only

  end module m_fock

!!========================================================================
!!>>> module m_sector                                                  <<<
!!========================================================================

!!
!! @mod m_sector
!!
!! define global structures and arrays for subspace of atomic Hamiltonian
!!
  module m_sector
     use constants, only : dp
     use constants, only : zero, czero

     implicit none

!!
!! @struct Tf
!!
!! data structure for annihilation operator f or creation operator f^+,
!!
!!     < alpha | f | beta > or < alpha | f^+ | beta >
!!
!! where | alpha > and | beta > are the atomic eigenstates in the
!! given subspace labelled by good quantum numbers
!!
     private :: Tf
     type Tf

         ! dimension of operator matrix, n x m
         integer :: n
         integer :: m

         ! memory space for operator matrix
         real(dp), allocatable :: val(:,:)

     end type Tf

!!
!! @struct Ts
!!
!! data structure for subspace of atomic Hamiltonian. sometimes we
!! call subspace as sector
!!
     public :: Ts
     type Ts

         ! global index of the first Fock state in this subspace
         integer :: istart

         ! dimension of this subspace
         ! how many Fock states are there in this subspace
         integer :: ndim

         ! number of fermion operators
         ! it is actually equal to norbs
         integer :: nops

         ! we just use N, Sz, Jz, and PS (AP) to label the subspaces.
         ! they are the so-called good quantum numbers. the good quantum
         ! number AP is only used in the automatic partition algorithm.

         ! total number of electrons: N
         integer :: nele

         ! z component of spin: Sz
         integer :: sz

         ! z component of spin-orbit momentum: Jz
         integer :: jz

         ! SU(2) good quantum number: PS (only for ictqmc = 4)
         ! see: Phys. Rev. B 86, 155158 (2012)
         !
         ! or
         !
         ! auxiliary good quantum number: AP (only for ictqmc = 6)
         ! see: Comput. Phys. Commun. 200, 274 (2016)
         integer :: ps

         ! collection of global indices of Fock states for this subspace
         integer, allocatable  :: basis(:)

         ! pointer to (or index of) the next subspace after a fermion
         ! operator acts on this subspace
         !
         !     next(nops,0) for annihilation
         !     next(nops,1) for creation operators
         !
         ! if it is -1, it means the next subspace is null (outside of
         ! the Hilbert space). otherwise, it is the index of next subspace
         integer, allocatable  :: next(:,:)

         ! eigenvalues of atomic Hamiltonian in this subspace
         real(dp), allocatable :: eval(:)

         ! eigenvectors of atomic Hamiltonian in this subspace
         ! since Hamiltonian must be real, then it is real as well
         real(dp), allocatable :: evec(:,:)

         ! atomic Hamiltonian in this subspace
         complex(dp), allocatable :: hmat(:,:)

         ! annihilation operator f or creation operator f^+ matrix
         !
         !     < alpha | f | beta > or < alpha | f^+ | beta >
         !
         ! where | alpha > and | beta > are the atomic eigenstates.
         !
         !     fmat(nops,0) for annihilation
         !     fmat(nops,1) for creation operators
         !
         ! if the corresponding matrix element of next(:,:) is -1, then
         ! the Tf is invalid (not allocated)
         type(Tf), allocatable :: fmat(:,:)

     end type Ts

!!
!! @var nsectors
!!
!! number of subspaces (sectors)
!!
     integer, public, save  :: nsectors

!!
!! @var max_dim_sect
!!
!! maximum dimension of subspaces (sectors)
!!
     integer, public, save  :: max_dim_sect

!!
!! @var ave_dim_sect
!!
!! average dimension of subspaces (sectors)
!!
     real(dp), public, save :: ave_dim_sect

!!
!! @var sectors
!!
!! an array of structs that stores all the subspaces (sectors)
!!
     type(Ts), public, save, allocatable :: sectors(:)

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
!! allocate memory for annihilation operator f or creation operator f^+
!!
  subroutine cat_alloc_fmat(one_fmat)
     implicit none

!! external arguments
     ! struct for annihilation operator f or creation operator f^+
     ! we have to make sure one_fmat%n and one_fmat%m are valid
     type(Tf), intent(inout) :: one_fmat

!! local variables
     ! the status flag
     integer :: istat

!! [body

     ! allocate memory
     allocate(one_fmat%val(one_fmat%n,one_fmat%m), stat=istat)

     ! check status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_fmat', &
             & 'can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize it
     one_fmat%val = zero

!! body]

     return
  end subroutine cat_alloc_fmat

!!
!! @sub cat_alloc_sector
!!
!! allocate memory for a subspace
!!
  subroutine cat_alloc_sector(one_sector)
     implicit none

!! external arguments
     ! this subspace
     type(Ts), intent(inout) :: one_sector

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
         call s_print_error('cat_alloc_sector', &
             & 'can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialize them
     one_sector%basis = 0
     one_sector%next  = 0
     one_sector%eval  = zero
     one_sector%evec  = zero
     one_sector%hmat  = czero

     ! initialize fmat one by one
     ! memory of one_sector%fmat%val should be allocated elsewhere
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
!! allocate memory for subspaces
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
         call s_print_error('cat_alloc_sectors', &
             & 'can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! initialization of each subspaces should be done elsewhere

!! body]

     return
  end subroutine cat_alloc_sectors

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!
!! @sub cat_free_fmat
!!
!! deallocate memory for annihilation operator f or creation operator f^+
!!
  subroutine cat_free_fmat(one_fmat)
     implicit none

!! external arguments
     ! annihilation operator f or creation operator f^+
     type(Tf), intent(inout) :: one_fmat

!! [body

     if ( allocated(one_fmat%val) ) deallocate(one_fmat%val)
     !
     one_fmat%n = 0
     one_fmat%m = 0

!! body]

     return
  end subroutine cat_free_fmat

!!
!! @sub cat_free_sector
!!
!! deallocate memory for a subspace
!!
  subroutine cat_free_sector(one_sector)
     implicit none

!! external arguments
     ! this subspace
     type(Ts), intent(inout) :: one_sector

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
         !
         deallocate(one_sector%fmat)
     endif ! back if ( allocated(one_sector%fmat)  ) block

!! body]

     return
  end subroutine cat_free_sector

!!
!! @sub cat_free_sectors
!!
!! deallocate memory for subspaces
!!
  subroutine cat_free_sectors()
     implicit none

!! local variables
     ! loop index
     integer :: i

!! [body

     ! deallocate memory for an array of Ts
     if ( allocated(sectors) ) then
         do i=1,nsectors
             call cat_free_sector( sectors(i) )
         enddo ! over i={1,nsectors} loop
         !
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
!! it contains some physical quantities that defined in single particle
!! basis, including:
!!     crystal field splitting,
!!     spin-orbit coupling,
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
!! Coulomb interaction U, a rank-4 tensor
!!
     complex(dp), public, allocatable, save :: umat(:,:,:,:)

!!
!! @var cmat
!!
!! crystal field splitting (CFS)
!!
     complex(dp), public, allocatable, save :: cmat(:,:)

!!
!! @var smat
!!
!! spin-orbit coupling (SOC)
!!
     complex(dp), public, allocatable, save :: smat(:,:)

!!
!! @var emat
!!
!! onsite energy (CFS + SOC) of impurity (emat = cmat + smat)
!!
     complex(dp), public, allocatable, save :: emat(:,:)

!!
!! @var tmat
!!
!! transformation matrix from original basis to natural eigenbasis
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
         call s_print_error('cat_alloc_spmat', &
             & 'can not allocate enough memory')
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
