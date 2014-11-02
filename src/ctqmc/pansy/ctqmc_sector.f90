!!!-------------------------------------------------------------------------
!!! project : pansy
!!! program : m_sect module
!!!           m_sect@alloc_one_fmat
!!!           m_sect@dealloc_one_fmat
!!!           m_sect@alloc_one_sect
!!!           m_sect@dealloc_one_sect
!!!           m_sect@ctqmc_allocate_memory_sect
!!!           m_sect@ctqmc_deallocate_memory_sect
!!!           m_sect@ctqmc_make_string
!!! source  : mod_control.f90
!!! type    : module
!!! authors : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014  by yilin wang
!!!           07/19/2014  by yilin wang
!!!           08/18/2014  by yilin wang
!!!           11/01/2014  by yilin wang 
!!!           11/02/2014  by yilin wang 
!!! purpose : define data structure for good quantum numbers (GQNs) algorithm
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> m_sect: define the data structure for good quantum numbers (GQNs) algorithm
  module m_sect
     use constants, only : dp, zero
     use control, only : idoub, mkink, norbs
     use context, only : type_v, flvr_v
 
     implicit none
  
! F-matrix between any two sectors, it is a rectangle matrix
     type :: t_fmat

! dimensions
         integer :: n, m

! items
         real(dp), dimension(:,:), pointer :: item => null()

     end type t_fmat
  

! t_sect type contains all the information of a subspace of H_{loc} 
     type :: t_sect

! dimension
         integer :: ndim

! total number of electrons
         integer :: nelec 

! number of fermion operators
         integer :: nops

! start index of this sector
         integer :: istart

! eigenvalues
         real(dp), dimension(:), pointer :: eval => null() 

! next sector it points to when a fermion operator acts on this sector, F|i> --> |j>
! next_sector(nops,0:1), 0 for annihilation and 1 for creation operators, respectively
! it is -1 if goes outside of the Hilbert space, 
! otherwise, it is the index of next sector
         integer, dimension(:,:), pointer :: next_sect => null()

! F-matrix between this sector and all other sectors
! if this sector doesn't point to some other sectors, the pointer is null
! fmat(nops, 0:1), 0 for annihilation and 1 for creation operators, respectively
         type(t_fmat), dimension(:,:), pointer :: fmat => null()

! final products of matrices, which will be used to calculate nmat and nnmat
         real(dp), dimension(:,:,:), pointer :: fprod => null()

! matrix of occupancy operator c^{\dagger}c 
         real(dp), dimension(:,:,:), pointer :: occu => null()

! matrix of double occupancy operator c^{\dagger}cc^{\dagger}c
         real(dp), dimension(:,:,:,:), pointer :: doccu => null()

     end type t_sect
     
! some global variables
! status flag
     integer, private :: istat

! total number of sectors
     integer, public, save :: nsect

! maximal dimension of the sectors
     integer, public, save :: mdim_sect

! average dimension of the sectors
     real(dp), public, save :: adim_sect

! array of t_sect contains all the sectors
     type(t_sect), public, save, allocatable :: sectors(:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================
    
     public :: alloc_one_fmat
     public :: dealloc_one_fmat
     public :: alloc_one_sect
     public :: dealloc_one_sect
     public :: ctqmc_allocate_memory_sect
     public :: ctqmc_deallocate_memory_sect
     public :: ctqmc_make_string

  contains ! encapsulated functionality
  
!!>>> alloc_one_fmat: allocate one fmat
  subroutine alloc_one_fmat(fmat)
     implicit none
  
! external variables
     type(t_fmat), intent(inout) :: fmat
  
     allocate( fmat%item(fmat%n, fmat%m), stat=istat )
  
! check the status
     if ( istat /= 0 ) then
         call s_print_error('alloc_one_fmat', 'can not allocate enough memory')
     endif ! back if ( istat /=0 ) block

! initialize it
     fmat%item = zero
  
     return
  end subroutine alloc_one_fmat
  
!!>>> dealloc_one_fmat: deallocate one fmat
  subroutine dealloc_one_fmat(fmat)
     implicit none
  
! external variables
     type(t_fmat), intent(inout) :: fmat
  
     if ( associated(fmat%item) ) deallocate(fmat%item)
  
     return
  end subroutine dealloc_one_fmat
  
!!>>> alloc_one_sector: allocate memory for one sector
  subroutine alloc_one_sect(sect)
     implicit none
  
! external variables
     type(t_sect), intent(inout) :: sect
  
! local variables
     integer :: i, j
  
     allocate( sect%eval(sect%ndim),                       stat=istat )
     allocate( sect%next_sect(sect%nops, 0:1),             stat=istat )
     allocate( sect%fmat(sect%nops, 0:1),                  stat=istat )
     allocate( sect%fprod(sect%ndim, sect%ndim, 2),        stat=istat )
     allocate( sect%occu(sect%ndim, sect%ndim, sect%nops), stat=istat )

     if ( idoub == 2 ) then
         allocate( sect%doccu(sect%ndim, sect%ndim, sect%nops, sect%nops), stat=istat )
         sect%doccu = zero
     endif ! back if ( idoub == 2 ) block

! check the status
     if ( istat /= 0 ) then
         call s_print_error('alloc_one_sect', 'can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
 
! initialize them
     sect%eval = zero
     sect%next_sect = 0
     sect%fprod = zero
     sect%occu = zero
  
! initialize fmat one by one
     do i=1,sect%nops 
         do j=0,1
             sect%fmat(i,j)%n = 0
             sect%fmat(i,j)%m = 0
         enddo ! over j={0,1} loop
     enddo ! over i={1,sect%nops} loop
  
     return
  end subroutine alloc_one_sect
  
!!>>> dealloc_one_sector: deallocate memory for one sector
  subroutine dealloc_one_sect(sect)
     implicit none
  
! external variables
     type(t_sect), intent(inout) :: sect
  
! local variables  
     integer :: i, j
  
     if ( associated(sect%eval) )        deallocate(sect%eval)
     if ( associated(sect%next_sect) )   deallocate(sect%next_sect)
     if ( associated(sect%fprod) )       deallocate(sect%fprod)
     if ( associated(sect%occu) )        deallocate(sect%occu)
     if ( associated(sect%doccu) )       deallocate(sect%doccu)
  
! deallocate fmat one by one
     do i=1,sect%nops
         do j=0,1
             call dealloc_one_fmat(sect%fmat(i,j))
         enddo ! over j={0,1} loop
     enddo ! over i={1,sect%nops} loop
  
     return
  end subroutine dealloc_one_sect

!!>>> ctqmc_allocate_memory_sect: allocate memory for sector related variables
  subroutine ctqmc_allocate_memory_sect()
     implicit none

! local variables
     integer :: i

! allocate memory
     allocate( sectors(nsect), stat=istat )

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_sect', &
                          'can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     do i=1,nsect
         sectors(i)%ndim = 0
         sectors(i)%nelec = 0
         sectors(i)%nops = norbs
         sectors(i)%istart = 0
     enddo ! over i={1,nsect} loop

     return
  end subroutine ctqmc_allocate_memory_sect

!!>>> ctqmc_deallocate_memory_sect: deallocate memory for sector related variables
  subroutine ctqmc_deallocate_memory_sect()
     implicit none

! local variables
     integer :: i

     if ( allocated(sectors) ) then
! first, loop over all the sectors and deallocate their component's memory
         do i=1,nsect
             call dealloc_one_sect(sectors(i))
         enddo ! over i={1,nsect} loop
! then, deallocate memory of sectors itself
         deallocate(sectors)
     endif ! back if ( allocated(sectors) ) block

     return
  end subroutine ctqmc_deallocate_memory_sect

!!>>> ctqmc_make_string: subroutine used to build a string
  subroutine ctqmc_make_string(csize, index_t_loc, is_string, string)
     implicit none

! external variables
! number of fermion operators for current diagram
     integer, intent(in) :: csize

! address index of fermion operators
     integer, intent(in) :: index_t_loc(mkink)

! whether it is a string ?
     logical, intent(out) :: is_string(nsect)

! string index
     integer, intent(out) :: string(csize+1, nsect)

! local variables
! sector index: from left direction
     integer :: left
     integer :: curr_sect_l
     integer :: next_sect_l

! sector index: from right direction
     integer :: right
     integer :: curr_sect_r
     integer :: next_sect_r

! flavor and type of fermion operators
     integer :: vf
     integer :: vt

! loop index
     integer :: i, j

     is_string = .true.
     string = -1

! we build a string from right to left, that is, beta <------- 0
! begin with S1: F1(S1)-->S2, F2(S2)-->S3, ... ,Fk(Sk)-->S1
! if find some Si==-1, cycle this sector immediately
     do i=1,nsect
         curr_sect_l = i
         curr_sect_r = i
         next_sect_l = i
         next_sect_r = i
         left = 0
         right = csize + 1
         do j=1,csize
             if ( mod(j,2) == 1 ) then
                 left = left + 1
                 string(left,i) = curr_sect_l
                 vt = type_v( index_t_loc(left) )
                 vf = flvr_v( index_t_loc(left) ) 
                 next_sect_l = sectors(curr_sect_l)%next_sect(vf,vt)
                 if ( next_sect_l == -1 ) then
                     is_string(i) = .false. 
                     EXIT   ! finish check, exit
                 endif ! back if ( next_sect_l == - 1 ) block
                 curr_sect_l = next_sect_l
             else
                 right = right - 1
                 vt = type_v( index_t_loc(right) )
                 vf = flvr_v( index_t_loc(right) ) 
                 vt = mod(vt+1,2)
                 next_sect_r = sectors(curr_sect_r)%next_sect(vf,vt)
                 if ( next_sect_r == -1 ) then
                     is_string(i) = .false. 
                     EXIT   ! finish check, exit
                 endif ! back if ( next_sect_r == -1 ) block
                 string(right,i) = next_sect_r
                 curr_sect_r = next_sect_r
             endif ! back if ( mod(j,2) == 1 ) block
         enddo ! over j={1,csize} loop

! if it doesn't form a string, we cycle it, go to the next sector
         if ( .not. is_string(i) ) then
             cycle
         endif ! back if ( .not. is_string(i) ) block

! add the last sector to string, and check whether string(csize+1,i) == string(1,i)
! important for csize = 0
         string(csize+1,i) = i

! this case will generate a non-diagonal block, it will not contribute to the trace 
         if ( next_sect_r /= next_sect_l ) then
             is_string(i) = .false.
         endif ! back if ( next_sect_r /= next_sect_l ) block
     enddo ! over i={1,nsect} loop

     return
  end subroutine ctqmc_make_string
 
  end module m_sect
