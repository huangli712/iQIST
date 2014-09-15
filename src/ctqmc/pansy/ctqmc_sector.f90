!!!-------------------------------------------------------------------------
!!! project : pansy
!!! program : m_sector module
!!! source  : mod_control.f90
!!! type    : module
!!! authors : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014
!!!           07/19/2014
!!!           08/18/2014
!!! purpose : define data structure for good quantum numbers (GQNs) algorithm
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> m_sector: define data structure for good quantum numbers (GQNs) algorithm
  module m_sector
     use constants, only : dp, zero
     use control, only : mkink, idoub, norbs
     use context, only : type_v, flvr_v
 
     implicit none
  
! the fmat between any two sectors, it is just a matrix
     type :: t_fmat
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

! the eigenvalues
         real(dp), pointer :: myeigval(:) 

! the next sector it points to when a fermion operator acts on this sector
! -1: outside of the Hilbert space, otherwise, it is the index of next sector
! next_sector(nops,0:1), 0 for annihilation and 1 for creation operators, respectively
! F|i> --> |j>
         integer, pointer :: next_sector(:,:)

! the fmat between this sector and all other sectors
! if this sector doesn't point to some other sectors, the pointer is null
! mymfat(nops, 0:1), 0 for annihilation and 1 for creation operators, respectively
         type(t_fmat), pointer :: myfmat(:,:)

! the final product matrices, which will be used to calculate the nmat
         real(dp), pointer :: final_product(:,:,:)

! matrices of occupancy operator c^{\dagger}c 
         real(dp), pointer :: occu(:,:,:)

! matrices of double occupancy operator c^{\dagger}cc^{\dagger}c
         real(dp), pointer :: double_occu(:,:,:,:)

     end type t_sector
     
! some global variables
! status flag
     integer, private :: istat

! the total number of sectors
     integer, public, save :: nsectors

! the max dimension of the sectors
     integer, public, save :: max_dim_sect

! the average dimension of the sectors
     real(dp), public, save :: ave_dim_sect

! the array contains all the sectors
     type(t_sector), public, save, allocatable :: sectors(:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================
    
     public :: nullify_one_fmat
     public :: alloc_one_fmat
     public :: dealloc_one_fmat
     public :: nullify_one_sector
     public :: alloc_one_sector
     public :: dealloc_one_sector
     public :: ctqmc_allocate_memory_sect
     public :: ctqmc_deallocate_memory_sect
     public :: ctqmc_make_string

  contains ! encapsulated functionality
  
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
  
     allocate(one_fmat%item(one_fmat%n, one_fmat%m), stat=istat)
  
! check the status
     if ( istat /= 0 ) then
         call s_print_error('alloc_one_fmat', 'can not allocate enough memory')
     endif

! initialize it
     one_fmat%item = zero
  
     return
  end subroutine alloc_one_fmat
  
!!>>> dealloc_one_fmat: deallocate one fmat
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
  
     nullify( one_sector%myeigval )
     nullify( one_sector%next_sector )
     nullify( one_sector%myfmat )
     nullify( one_sector%final_product )
     nullify( one_sector%occu )
     nullify( one_sector%double_occu )
  
     return
  end subroutine nullify_one_sector
  
!!>>> alloc_one_sector: allocate memory for one sector
  subroutine alloc_one_sector(one_sector)
     implicit none
  
! external variables
     type(t_sector), intent(inout) :: one_sector
  
! local variables
     integer :: i, j
  
     allocate( one_sector%myeigval(one_sector%ndim),           stat=istat )
     allocate( one_sector%next_sector(one_sector%nops,0:1),    stat=istat )
     allocate( one_sector%myfmat(one_sector%nops,0:1),         stat=istat )

     allocate( one_sector%final_product(one_sector%ndim, & 
                                        one_sector%ndim, 2),   stat=istat )

     allocate( one_sector%occu(one_sector%ndim, &
                               one_sector%ndim, & 
                               one_sector%nops   ),            stat=istat )

     if (idoub == 2) then
         allocate( one_sector%double_occu(one_sector%ndim, & 
                                          one_sector%ndim, &
                                          one_sector%nops, &
                                          one_sector%nops   ), stat=istat )

         one_sector%double_occu = zero
     endif

! check the status
     if ( istat /= 0 ) then
         call s_print_error('alloc_one_sector', 'can not allocate enough memory')
     endif
 
! initialize them
     one_sector%myeigval = zero
     one_sector%next_sector = 0
     one_sector%final_product = zero
     one_sector%occu = zero
  
! initialize myfmat one by one
     do i=1, one_sector%nops 
         do j=0, 1
             one_sector%myfmat(i,j)%n = 0
             one_sector%myfmat(i,j)%m = 0
             call nullify_one_fmat(one_sector%myfmat(i,j))
         enddo
     enddo
  
     return
  end subroutine alloc_one_sector
  
!!>>> dealloc_one_sector: deallocate memory for onespace
  subroutine dealloc_one_sector(one_sector)
     implicit none
  
! external variables
     type(t_sector), intent(inout) :: one_sector 
  
! local variables  
     integer :: i, j
  
     if ( associated(one_sector%myeigval) )        deallocate(one_sector%myeigval)
     if ( associated(one_sector%next_sector) )     deallocate(one_sector%next_sector)
     if ( associated(one_sector%final_product) )   deallocate(one_sector%final_product)
     if ( associated(one_sector%occu) )            deallocate(one_sector%occu)
     if ( associated(one_sector%double_occu) )     deallocate(one_sector%double_occu)
  
! deallocate myfmat one by one
     do i=1, one_sector%nops
         do j=0,1
             call dealloc_one_fmat(one_sector%myfmat(i,j))
         enddo
     enddo 
  
     return
  end subroutine dealloc_one_sector

!!>>> ctqmc_allocate_memory_sect: allocate memory for 
!!>>> sect-related variables
  subroutine ctqmc_allocate_memory_sect()
     implicit none

! local variables
     integer :: i

! allocate memory
     allocate( sectors(nsectors), stat=istat )

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_sect', &
                            'can not allocate enough memory')
     endif

! initialize them
     do i=1, nsectors
         sectors(i)%ndim = 0
         sectors(i)%nelectron = 0
         sectors(i)%nops = norbs
         sectors(i)%istart = 0
         call nullify_one_sector(sectors(i))
     enddo 

     return
  end subroutine ctqmc_allocate_memory_sect

!!>>> ctqmc_deallocate_memory_sect: deallocate memory for 
!!>>> sect-related variables
  subroutine ctqmc_deallocate_memory_sect()
     implicit none

! local variables
     integer :: i

     if ( allocated(sectors) ) then
! first, loop over all the sectors and deallocate their memory
         do i=1, nsectors
             call dealloc_one_sector(sectors(i))
         enddo
! then, deallocate memory of sect
         deallocate(sectors)
     endif

     return
  end subroutine ctqmc_deallocate_memory_sect

!!>>> ctqmc_make_string: subroutine used to build a string
  subroutine ctqmc_make_string(csize, index_t_loc, is_string, string)
     implicit none

! external variables
! the number of fermion operators
     integer, intent(in) :: csize

! the address index of fermion operators
     integer, intent(in) :: index_t_loc(mkink)

! whether it is a string
     logical, intent(out) :: is_string(nsectors)

! the string
     integer, intent(out) :: string(csize+1, nsectors)

! local variables
! sector index
     integer :: curr_sect_left
     integer :: next_sect_left
     integer :: curr_sect_right
     integer :: next_sect_right
     integer :: left
     integer :: right

! flvr and type of fermion operators
     integer :: vf
     integer :: vt

! loop index
     integer :: i,j

     is_string = .true.
     string = -1

! we build a string from right to left, that is,  beta <------- 0
! build the string from the beginning sector, that is:
! S_a1(q1)-->q2, S_a2(q2)-->q3, ... S_ai(qi)-->qi+1, ..., Sak(qk)-->q1
! if we find some qi==0, we cycle this sector immediately
     do i=1,nsectors
         curr_sect_left = i
         curr_sect_right = i
         next_sect_left = i
         next_sect_right = i
         left = 0
         right = csize + 1
         do j=1,csize
             if ( mod(j,2) == 1) then
                 left = left + 1
                 string(left,i) = curr_sect_left 
                 vt = type_v( index_t_loc(left) )
                 vf = flvr_v( index_t_loc(left) ) 
                 next_sect_left = sectors(curr_sect_left)%next_sector(vf,vt)
                 if (next_sect_left == -1 ) then
                     is_string(i) = .false. 
                     EXIT   ! finish check, exit
                 endif
                 curr_sect_left = next_sect_left
             else
                 right = right - 1
                 vt = type_v( index_t_loc(right) )
                 vf = flvr_v( index_t_loc(right) ) 
                 vt = mod(vt+1,2)
                 next_sect_right = sectors(curr_sect_right)%next_sector(vf,vt)
                 if (next_sect_right == -1 ) then
                     is_string(i) = .false. 
                     EXIT   ! finish check, exit
                 endif
                 string(right,i) = next_sect_right
                 curr_sect_right = next_sect_right
             endif
         enddo 

! if it doesn't form a string, we cycle it, go to the next sector
         if (is_string(i) .eqv. .false.) then
             cycle
         endif
! add the last sector to string, and check whether string(csize+1,i) == string(1,i)
! important for csize = 0
         string(csize+1,i) = i
! this case will generate a non-diagonal block, it will not contribute to trace 
         if ( next_sect_right /= next_sect_left ) then
             is_string(i) = .false.
         endif
     enddo ! over i={1,nsectors} loop

     return
  end subroutine ctqmc_make_string
 
  end module m_sector
