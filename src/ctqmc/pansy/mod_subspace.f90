!>>> this module defines the data structure for the good quantum number algorithm
  module subspace
     use constants

     implicit none

! the fmat between any two subspaces, it is just a matrix
     type fmat
! the dimension 
         integer :: n

         integer :: m

! the items of the matrix
         real(dp), pointer :: item(:,:)

     end type fmat

! a sector defines a superstate consists of all the eigenstates labeled by same quantum number
! for example: {N, Sz}, {N, Jz}
     type sector
! the index 
         integer :: indx

! the dimension
         integer :: ndim

! total number of electrons for this sector
         integer :: nelectron

! number of operators
         integer :: nops

! the next sector it points to when a fermi operator acts on this sector 
! 0 for outside of the space, otherwise, it is the index of sector
! next_sector(nops,0:1), 0 for annihilation and 1 for creation operators, respectively
         integer, pointer :: next_sector(:,:)

! the eigenvalues of this sector
! eigenvalue(ndim)
         real(dp), pointer :: eigenvalue(:)
 
! the fmat between this sector and all other sectors 
! if this sector doesn't point to some other sectors, the pointer is null
! mymfat(nops, 0:1), 0 for annihilation and 1 for creation operators, respectively
         type(fmat), pointer :: myfmat(:,:)
 
     end type sector
     
     contains

!>>> nullify one fmat
     subroutine nullify_one_fmat(one_fmat)
         implicit none

! external variables
         type(fmat), intent(inout) :: one_fmat

         nullify(one_fmat%item)

         return
     end subroutine nullify_one_fmat

!>>> allocate one fmat
     subroutine alloc_one_fmat(one_fmat)
         implicit none

! external variables
         type(fmat), intent(inout) :: one_fmat 

         allocate(one_fmat%item(one_fmat%n, one_fmat%m))

! initialize it
         one_fmat%item = zero

         return
     end subroutine alloc_one_fmat

!>>> deallocate one fmat
     subroutine dealloc_one_fmat(one_fmat)
         implicit none

! external variables
         type(fmat), intent(inout) :: one_fmat

         if (associated(one_fmat%item)) deallocate(one_fmat%item) 

         return
     end subroutine dealloc_one_fmat
     
!>>> nullify one sector
     subroutine nullify_one_sector(one_sector)
         implicit none

! external variables
         type(sector), intent(inout) :: one_sector

         nullify( one_sector%next_sector ) 
         nullify( one_sector%eigenvalue )
         nullify( one_sector%myfmat )

         return
     end subroutine nullify_one_sector

!>>> allocate one sector
     subroutine alloc_one_sector(one_sector)
         implicit none

! external variables
         type(sector), intent(inout) :: one_sector

! local variables
         integer :: i, j

         allocate( one_sector%next_sector(one_sector%nops, 0:1) )
         allocate( one_sector%eigenvalue(one_sector%ndim) )
         allocate( one_sector%myfmat(one_sector%nops, 0:1) )

! initialize them
         one_sector%next_sector = 0
         one_sector%eigenvalue = zero

! initialize mymfat one by one
         do i=1, one_sector%nops
             do j=0, 1
                 one_sector%myfmat(i,j)%n = 0
                 one_sector%myfmat(i,j)%m = 0
                 call nullify_one_fmat(one_sector%myfmat(i,j))
             enddo
         enddo

         return
     end subroutine alloc_one_sector

!>>> deallocate one sector
     subroutine dealloc_one_sector(one_sector)
         implicit none

! external variables
         type(sector), intent(inout) :: one_sector

! local variables
         integer :: i, j

         if ( associated(one_sector%next_sector) ) deallocate(one_sector%next_sector)
         if ( associated(one_sector%eigenvalue) )  deallocate(one_sector%eigenvalue) 

! dellocate fmat one by one
         if ( associated(one_secotr%myfmat) ) then
             do i=0, 1
                 do j=1, one_sector%nops
                     call dealloc_one_fmat(one_sector%myfmat(j,i))                  
                 enddo
             enddo
         endif
  
         return
     end subroutine dealloc_one_sector

  end module subspace
