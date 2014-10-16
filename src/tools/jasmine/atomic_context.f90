!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : m_full  module
!!!           m_sector  module
!!!           m_spmat            module
!!! source  : atomic_context.f90
!!! type    : modules
!!! authors : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : define data structure for good quantum numbers (GQNs) algorithm
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> Fock basis of full Hilbert space
  module m_full
     use constants, only : dp, zero, czero
     use control, only : norbs, ncfgs
  
     implicit none
  
! dimension of subspace of total electron N
     integer, public, allocatable, save :: dim_sub_n(:)
  
! binary form of Fock basis
     integer, public, allocatable, save :: bin_basis(:,:)
  
! decimal form of Fock basis
     integer, public, allocatable, save :: dec_basis(:)
  
! index of Fock basis, given their decimal number
     integer, public, allocatable, save :: index_basis(:)

! atomic Hamiltonian (CF + SOC + CU)
     complex(dp), public, allocatable, save :: hmat(:,:)
  
! eigenvalues of hmat
     real(dp), public, allocatable, save :: eigval(:)
  
! eigenvectors of hmat
     real(dp), public, allocatable, save :: eigvec(:, :)
  
! fmat for annihilation fermion operators
     real(dp), public, allocatable, save :: fmat(:,:,:)
  
! occupany number for atomic eigenstates
     real(dp), public, allocatable, save :: occu(:,:)

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: alloc_m_basis_fullspace
     public :: dealloc_m_basis_fullspace
     public :: alloc_m_glob_fullspace 
     public :: dealloc_m_glob_fullspace 

  contains
  
!!>>> alloc_m_basis_fullspace: allocate memory for these matrices
  subroutine alloc_m_basis_fullspace()
     implicit none
  
! allocate them
     allocate( dim_sub_n(0:norbs),        stat=istat )
     allocate( bin_basis(norbs, ncfgs),   stat=istat )
     allocate( dec_basis(ncfgs),          stat=istat )
     allocate( index_basis(0:ncfgs-1),    stat=istat )
  
! check the status 
     if ( istat /= 0 ) then
         call s_print_error('alloc_m_basis_fullspace', &
                            'can not allocate enough memory')
     endif

! initialize them
     dim_sub_n = 0
     bin_basis = 0
     dec_basis = 0
     index_basis = 0
  
     return
  end subroutine alloc_m_basis_fullspace
  
!!>>> dealloc_m_basis_fullspace: deallocate memory for these matrices
  subroutine dealloc_m_basis_fullspace()
     implicit none
  
! deallocate them
     if (allocated(dim_sub_n))   deallocate(dim_sub_n)
     if (allocated(bin_basis))   deallocate(bin_basis)
     if (allocated(dec_basis))   deallocate(dec_basis) 
     if (allocated(index_basis)) deallocate(index_basis) 
  
     return
  end subroutine dealloc_m_basis_fullspace

!!>>> alloc_m_glob_fullspace: allocate memory for m_glob_fullspace module
  subroutine alloc_m_glob_fullspace()
     implicit none

     allocate( hmat(ncfgs, ncfgs),              stat=istat )
     allocate( eigval(ncfgs),              stat=istat )
     allocate( eigvec(ncfgs, ncfgs),       stat=istat )
     allocate( fmat(ncfgs, ncfgs, norbs),  stat=istat )
     allocate( occu(ncfgs, ncfgs),          stat=istat )
  
! check status
     if ( istat /= 0 ) then
         call s_print_error('alloc_m_glob_fullspace', 'can not allocate enough memory')
     endif

! init them
     hmat = czero
     eigval = zero
     eigvec = zero
     fmat = zero
     occu = zero
  
     return
  end subroutine alloc_m_glob_fullspace
  
!!>>> dealloc_m_glob_fullspace: deallocate memory for m_glob_fullspace module
  subroutine dealloc_m_glob_fullspace()
     implicit none
  
     if(allocated(hmat))        deallocate(hmat)
     if(allocated(eigval)) deallocate(eigval)
     if(allocated(eigvec)) deallocate(eigvec)
     if(allocated(fmat))   deallocate(fmat)
     if(allocated(occu))    deallocate(occu)
  
     return
  end subroutine dealloc_m_glob_fullspace
 
  end module m_full

!!>>> data structure for good quantum numbers (GQNs) algorithm
  module m_sector
     use constants, only : dp, zero, czero
  
     implicit none
  
! the fmat between any two sectors, it is just a matrix
     type t_fmat
! the dimension
         integer :: n, m

! the items of the matrix
         real(dp), pointer :: item(:,:)
     end type t_fmat
  
! one sector
     type :: t_sector 
! the dimension of this sector
         integer :: ndim

! total number of electrons n
         integer :: nelectron 

! number of fermion operators
         integer :: nops

! the start index of this sector
         integer :: istart

! the Fock basis index of this sector
         integer, pointer :: basis(:)

! the Hamiltonian of this sector
         complex(dp), pointer :: ham(:,:)

! the eigenvalues
         real(dp), pointer :: eigval(:) 

! the eigenvectors, Hamiltonian must be real
         real(dp), pointer :: eigvec(:,:) 

! the next sector it points to when a fermion operator acts on this sector
! -1: outside of the Hilbert space, otherwise, it is the index of next sector
! next_sector(nops,0:1), 0 for annihilation and 1 for creation operators, respectively
         integer, pointer :: next(:,:)

! the fmat between this sector and all other sectors
! if this sector doesn't point to some other sectors, the pointer is null
! mymfat(nops, 0:1), 0 for annihilation and 1 for creation operators, respectively
         type(t_fmat), pointer :: fmat(:,:)
     end type t_sector

! number of sectors
     integer, public, save :: nsectors
  
! maximum dimension of sectors
     integer, public, save :: max_dim_sect
  
! average dimension of sectors
     real(dp), public, save :: ave_dim_sect
  
! all the sectors
     type(t_sector), public, allocatable, save :: sectors(:)

! status of allocating memory
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: nullify_one_fmat
     public :: alloc_one_fmat
     public :: dealloc_one_fmat
     public :: nullify_one_sector
     public :: alloc_one_sector
     public :: dealloc_one_sector

     public :: alloc_m_glob_sectors
     public :: dealloc_m_glob_sectors

  contains
  
!!>>> nullify_one_fmat: nullify one fmat
  subroutine nullify_one_fmat(one_fmat)
     implicit none
  
! external variables
     type(t_fmat), intent(inout) :: one_fmat
  
     nullify(one_fmat%item)
  
     return
  end subroutine nullify_one_fmat
  
!!>>> alloc_one_fmat: allocate one fmat
  subroutine alloc_one_fmat(one_fmat)
     implicit none
  
! external variables
     type(t_fmat), intent(inout) :: one_fmat
  
     allocate( one_fmat%item(one_fmat%n, one_fmat%m),  stat=istat )
  
! check status
     if ( istat /= 0 ) then
         call s_print_error('alloc_one_fmat', 'can not allocate enough memory')
     endif
! initialize it
     one_fmat%item = zero
  
     return
  end subroutine alloc_one_fmat
  
!>>> dealloc_one_fmat: deallocate one fmat
  subroutine dealloc_one_fmat(one_fmat)
     implicit none
  
! external variables
     type(t_fmat), intent(inout) :: one_fmat
  
     if ( associated(one_fmat%item) ) deallocate(one_fmat%item)
  
     return
  end subroutine dealloc_one_fmat
  
!!>>> nullify_one_sector: nullify one sector
  subroutine nullify_one_sector(one_sector)
     implicit none
  
! external variables
     type(t_sector), intent(inout) :: one_sector
  
     nullify( one_sector%basis )
     nullify( one_sector%ham )
     nullify( one_sector%eigval )
     nullify( one_sector%eigvec )
     nullify( one_sector%next )
     nullify( one_sector%fmat )
  
     return
  end subroutine nullify_one_sector
  
!>>> alloc_one_sector: allocate memory for one sector
  subroutine alloc_one_sector(one_sector)
     implicit none
  
! external variables
     type(t_sector), intent(inout) :: one_sector
  
! local variables
     integer :: i, j
  
     allocate( one_sector%basis(one_sector%ndim),                   stat=istat ) 
     allocate( one_sector%ham(one_sector%ndim, one_sector%ndim),    stat=istat ) 
     allocate( one_sector%eigval(one_sector%ndim),                  stat=istat )
     allocate( one_sector%eigvec(one_sector%ndim, one_sector%ndim), stat=istat ) 
     allocate( one_sector%next(one_sector%nops,0:1),           stat=istat )
     allocate( one_sector%fmat(one_sector%nops,0:1),                stat=istat )
  
! check status
     if ( istat /= 0 ) then
         call s_print_error('alloc_one_sector', 'can not allocate enough memory')
     endif

! initialize them
     one_sector%basis = 0
     one_sector%ham = czero
     one_sector%eigval = zero
     one_sector%eigvec = zero
     one_sector%next = 0
  
! initialize myfmat one by one
     do i=1, one_sector%nops 
        do j=0, 1
            one_sector%fmat(i,j)%n = 0
            one_sector%fmat(i,j)%m = 0
            call nullify_one_fmat(one_sector%fmat(i,j))
        enddo
     enddo
  
     return
  end subroutine alloc_one_sector
  
!!>>> dealloc_one_sector: deallocate memory for one sector
  subroutine dealloc_one_sector(one_sector)
     implicit none
  
! external variables
     type(t_sector), intent(inout) :: one_sector 
  
! local variables  
     integer :: i, j
  
     if (associated(one_sector%basis))      deallocate(one_sector%basis)
     if (associated(one_sector%ham))        deallocate(one_sector%ham)
     if (associated(one_sector%eigval))     deallocate(one_sector%eigval)
     if (associated(one_sector%eigvec))     deallocate(one_sector%eigvec)
     if (associated(one_sector%next))  deallocate(one_sector%next)
  
! deallocate myfmat one by one
     do i=1, one_sector%nops
         do j=0,1
             call dealloc_one_fmat(one_sector%fmat(i,j))
         enddo
     enddo 
  
     return
  end subroutine dealloc_one_sector

!!>>> alloc_m_glob_sectors: allocate memory for sectors
  subroutine alloc_m_glob_sectors()
     implicit none
  
! local variables
! loop index
     integer :: i
  
     allocate( sectors(nsectors),   stat=istat ) 
  
! check status
     if ( istat /= 0 ) then
         call s_print_error('alloc_m_glob_sectors', &
                            'can not allocate enough memory')
     endif
 
! nullify each sector one by one
     do i=1, nsectors
         call nullify_one_sector(sectors(i))
     enddo          
  
     return
  end subroutine alloc_m_glob_sectors
  
!!>>> dealloc_m_glob_sectors: deallocate memory of sectors
  subroutine dealloc_m_glob_sectors()
     implicit none
     
     integer :: i

! deallocate memory for pointers in t_sectors
! before deallocating sectors to avoid memory leak 
     do i=1, nsectors
         call dealloc_one_sector(sectors(i)) 
     enddo 
  
     if (allocated(sectors))  deallocate(sectors) 
     
     return
  end subroutine dealloc_m_glob_sectors

  end module m_sector

!!>>> single particle related matrices, including: 
!!>>> crystal field, spin-orbital coupling, Coulomb interaction U tensor
  module m_spmat
     use constants, only : dp, czero
     use control, only : norbs

     implicit none
  
! crystal field (CF)
     complex(dp), public, allocatable, save :: cmat(:,:) 
  
! spin-orbital coupling (SOC)
     complex(dp), public, allocatable, save :: smat(:,:)
  
! on-site energy (CF+SOC) of impurity
     complex(dp), public, allocatable, save :: emat(:,:)
  
! Coulomb interaction U tensor
     complex(dp), public, allocatable, save :: umat(:,:,:,:)
  
! the transformation matrix from origional basis to natural basis 
     complex(dp), public, allocatable, save :: tmat(:,:)
  
! the status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: alloc_m_spmat
     public :: dealloc_m_spmat

  contains
  
!!>>> alloc_m_spmat: allocate memory for these matrix
  subroutine alloc_m_spmat()
     implicit none
  
! allocate them
     allocate( cmat(norbs, norbs),               stat=istat )
     allocate( smat(norbs, norbs),              stat=istat )
     allocate( emat(norbs, norbs),             stat=istat )
     allocate( umat(norbs, norbs, norbs, norbs), stat=istat )
     allocate( tmat(norbs, norbs),           stat=istat )
  
! check the status
     if ( istat /= 0 ) then
         call s_print_error('alloc_m_spmat', 'can not allocate enough memory')
     endif

! initialize them
     cmat    = czero
     smat   = czero
     emat  = czero
     umat    = czero
     tmat= czero
  
     return
  end subroutine alloc_m_spmat
  
!!>>> dealloc_m_spmat: deallocate memory for these matrix
  subroutine dealloc_m_spmat()
     implicit none
  
! deallocate them
     if (allocated(cmat))      deallocate(cmat)
     if (allocated(smat))     deallocate(smat)
     if (allocated(emat))    deallocate(emat)
     if (allocated(umat))      deallocate(umat)
     if (allocated(tmat))  deallocate(tmat)
  
     return
  end subroutine dealloc_m_spmat
  
  end module m_spmat
