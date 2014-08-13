!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : m_sector
!!!           m_sector@nullify_one_fmat
!!!           m_sector@alloc_one_fmat
!!!           m_sector@dealloc_one_fmat
!!!           m_sector@nullify_one_sector
!!!           m_sector@alloc_one_sector
!!!           m_sector@dealloc_one_sector
!!! source  : mod_control.f90
!!! type    : modules
!!! authors : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014
!!! purpose : define data structure for good quantum number algorithm
!!! input   :
!!! output  :
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> data structure for good quantum number algorithm
  module m_sector
     use constants,  only: dp, zero, czero
  
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
         integer, pointer :: mybasis(:)

! the Hamiltonian of this sector
         complex(dp), pointer :: myham(:,:)

! the eigenvalues
         real(dp), pointer :: myeigval(:) 

! the eigenvectors, Hamiltonian must be real
         real(dp), pointer :: myeigvec(:,:) 

! the next sector it points to when a fermion operator acts on this sector
! -1: outside of the Hilbert space, otherwise, it is the index of next sector
! next_sector(nops,0:1), 0 for annihilation and 1 for creation operators, respectively
         integer, pointer :: next_sector(:,:)

! the fmat between this sector and all other sectors
! if this sector doesn't point to some other sectors, the pointer is null
! mymfat(nops, 0:1), 0 for annihilation and 1 for creation operators, respectively
         type(t_fmat), pointer :: myfmat(:,:)
     end type t_sector
     
! status of allocating memory
     integer, private :: istat

     contains
  
!!>>> nullify one fmat
     subroutine nullify_one_fmat(one_fmat)
        implicit none
  
! external variables
        type(t_fmat), intent(inout) :: one_fmat
  
        nullify(one_fmat%item)
  
        return
     end subroutine nullify_one_fmat
  
!!>>> allocate one fmat
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
  
     !>>> deallocate one fmat
     subroutine dealloc_one_fmat(one_fmat)
        implicit none
  
! external variables
        type(t_fmat), intent(inout) :: one_fmat
  
        if ( associated(one_fmat%item) ) deallocate(one_fmat%item)
  
        return
     end subroutine dealloc_one_fmat
  
!!>>> nullify one sector
     subroutine nullify_one_sector(one_sector)
        implicit none
  
! external variables
        type(t_sector), intent(inout) :: one_sector
  
        nullify( one_sector%mybasis )
        nullify( one_sector%myham )
        nullify( one_sector%myeigval )
        nullify( one_sector%myeigvec )
        nullify( one_sector%next_sector )
        nullify( one_sector%myfmat )
  
        return
     end subroutine nullify_one_sector
  
!>>> allocate memory for one sector
     subroutine alloc_one_sector(one_sector)
        implicit none
  
! external variables
        type(t_sector), intent(inout) :: one_sector
  
! local variables
        integer :: i, j
  
        allocate( one_sector%mybasis(one_sector%ndim),                   stat=istat ) 
        allocate( one_sector%myham(one_sector%ndim, one_sector%ndim),    stat=istat ) 
        allocate( one_sector%myeigval(one_sector%ndim),                  stat=istat )
        allocate( one_sector%myeigvec(one_sector%ndim, one_sector%ndim), stat=istat ) 
        allocate( one_sector%next_sector(one_sector%nops,0:1),           stat=istat )
        allocate( one_sector%myfmat(one_sector%nops,0:1),                stat=istat )
  
! check status
        if ( istat /= 0 ) then
            call s_print_error('alloc_one_sector', 'can not allocate enough memory')
        endif

! init them
        one_sector%mybasis = 0
        one_sector%myham = czero
        one_sector%myeigval = zero
        one_sector%myeigvec = zero
        one_sector%next_sector = 0
  
! init myfmat one by one
        do i=1, one_sector%nops 
           do j=0, 1
               one_sector%myfmat(i,j)%n = 0
               one_sector%myfmat(i,j)%m = 0
               call nullify_one_fmat(one_sector%myfmat(i,j))
           enddo
        enddo
  
        return
     end subroutine alloc_one_sector
  
!!>>> deallocate memory for onespace
     subroutine dealloc_one_sector(one_sector)
        implicit none
  
! external variables
        type(t_sector), intent(inout) :: one_sector 
  
! local variables  
        integer :: i, j
  
        if (associated(one_sector%mybasis))      deallocate(one_sector%mybasis)
        if (associated(one_sector%myham))        deallocate(one_sector%myham)
        if (associated(one_sector%myeigval))     deallocate(one_sector%myeigval)
        if (associated(one_sector%myeigvec))     deallocate(one_sector%myeigvec)
        if (associated(one_sector%next_sector))  deallocate(one_sector%next_sector)
  
! deallocate myfmat one by one
        do i=1, one_sector%nops
            do j=0,1
                call dealloc_one_fmat(one_sector%myfmat(i,j))
            enddo
        enddo 
  
        return
     end subroutine dealloc_one_sector

  end module m_sector
