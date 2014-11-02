!!!-------------------------------------------------------------------------
!!! project : manjushaka
!!! program : m_sector  module
!!!           m_sector@alloc_one_fmat
!!!           m_sector@dealloc_one_fmat
!!!           m_sector@alloc_one_sqrmat
!!!           m_sector@dealloc_one_sqrmat
!!!           m_sector@alloc_one_sect
!!!           m_sector@dealloc_one_sect
!!!           m_sector@ctqmc_allocate_memory_sect
!!!           m_sector@ctqmc_deallocate_memory_sect
!!!           m_sector@ctqmc_allocate_memory_occu
!!!           m_sector@ctqmc_deallocate_memory_occu
!!!           m_sector@ctqmc_make_trunc
!!!           m_sector@cat_trunc_sect
!!!           m_sector@ctqmc_make_occu
!!!           m_sector@ctqmc_make_string
!!! source  : ctqmc_sector.f90
!!! type    : module
!!! authors : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           07/19/2014 by yilin wang
!!!           08/09/2014 by yilin wang
!!!           08/13/2014 by yilin wang
!!!           08/20/2014 by yilin wang
!!!           11/02/2014 by yilin wang
!!! purpose : define data structure for good quantum numbers (GQNs) algorithm
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> m_sect: define the data structure for good quantum numbers (GQNs) algorithm
  module m_sect
     use constants, only : dp, zero, one, mystd, mytmp
     use control, only : idoub, itrun, myid, master
     use control, only : mkink, norbs, nmini, nmaxi
     use context, only : type_v, flvr_v

     use mmpi

     implicit none
  
! F-matrix between any two sectors, it is a rectangle matrix
     type :: t_fmat

! dimensions
         integer :: n, m

! items
         real(dp), dimension(:,:), pointer :: item => null()

     end type t_fmat

! square matrix
     type :: t_sqrmat

! dimension
         integer :: n

! items 
         real(dp), dimension(:,:), pointer :: item => null()

     end type t_sqrmat

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

! index of next sector, for truncating the Hilbert space of H_{loc}
         integer, dimension(:,:), pointer :: next_sect_t => null()
   
! F-matrix between this sector and all other sectors
! if this sector doesn't point to some other sectors, the pointer is null
! fmat(nops, 0:1), 0 for annihilation and 1 for creation operators, respectively
         type(t_fmat), dimension(:,:), pointer :: fmat => null()

     end type t_sect
    
! some global variables
! status flag
     integer, private :: istat

! total number of sectors
     integer, public, save :: nsect

! total number of sectors after truncating H_{loc}
     integer, public, save :: nsect_t
     
! maximal dimension of the sectors
     integer, public, save :: mdim_sect

! maximal dimension of the sectors after truncating H_{loc}
     integer, public, save :: mdim_sect_t

! average dimension of the sectors
     real(dp), public, save :: adim_sect

! average dimension of the sectors after truncating H_{loc}
     real(dp), public, save :: adim_sect_t

! array of t_sect contains all the sectors
     type(t_sect), public, save, allocatable :: sectors(:)

! probability of each sector, used to truncate high energy states
     real(dp), public, save, allocatable :: prob_sect(:)

! which sectors should be truncated ?
     logical, public, save, allocatable :: is_trunc(:)

! whether it forms a string ?
     logical, public, save, allocatable :: is_string(:,:)

! final product of matrices multiplications, which will be used to calculate nmat
     type(t_sqrmat), public, save, allocatable :: fprod(:,:)

! matrix of occupancy operator c^{\dagger}c 
     type(t_sqrmat), public, save, allocatable :: occu(:,:)

! matrix of double occupancy operator c^{\dagger}cc^{\dagger}c
     type(t_sqrmat), public, save, allocatable :: doccu(:,:,:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================
    
     public :: alloc_one_fmat
     public :: dealloc_one_fmat
     public :: alloc_one_sqrmat
     public :: dealloc_one_sqrmat
     public :: alloc_one_sect
     public :: dealloc_one_sect
     public :: ctqmc_allocate_memory_sect
     public :: ctqmc_deallocate_memory_sect
     public :: ctqmc_allocate_memory_occu
     public :: ctqmc_deallocate_memory_occu
     public :: ctqmc_make_trunc
     public :: cat_trunc_sect
     public :: ctqmc_make_string
     public :: ctqmc_make_occu 

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
     endif ! back if ( istat /= 0 ) block

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

!!>>> alloc_one_sqrmat: allocate one sqrmat
  subroutine alloc_one_sqrmat(sqrmat)
     implicit none
  
! external variables
     type(t_sqrmat), intent(inout) :: sqrmat
  
! allocate it
     allocate( sqrmat%item(sqrmat%n,sqrmat%n),  stat=istat )

! check the status
     if ( istat /= 0 ) then
         call s_print_error('alloc_one_sqrmat', 'can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block
 
! initialize it
     sqrmat%item = zero
  
     return
  end subroutine alloc_one_sqrmat
  
!!>>> dealloc_one_sqrmat: deallocate one sqrmat
  subroutine dealloc_one_sqrmat(sqrmat)
     implicit none
  
! external variables
     type(t_sqrmat), intent(inout) :: sqrmat
  
     if ( associated(sqrmat%item) ) deallocate(sqrmat%item)
  
     return
  end subroutine dealloc_one_sqrmat
 
!!>>> alloc_one_sector: allocate memory for one sector
  subroutine alloc_one_sect(sect)
     implicit none
  
! external variables
     type(t_sect), intent(inout) :: sect
  
! local variables
     integer :: i, j
  
! allocate them
     allocate( sect%eval(sect%ndim),            stat=istat )
     allocate( sect%next_sect(sect%nops,0:1),   stat=istat )
     allocate( sect%next_sect_t(sect%nops,0:1), stat=istat )
     allocate( sect%fmat(sect%nops,0:1),        stat=istat )

! check the status
     if ( istat /= 0 ) then
         call s_print_error('alloc_one_sector','can not allocate enough memory')
     endif ! back if ( istat /=0 ) block
  
! initialize them
     sect%eval = zero
     sect%next_sect = 0
     sect%next_sect_t = 0
  
! initialize fmat one by one
     do i=1,sect%nops 
         do j=0,1
             sect%fmat(i,j)%n = 0
             sect%fmat(i,j)%m = 0
             sect%fmat(i,j)%item => null()
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
  
     if ( associated(sect%eval) )          deallocate(sect%eval)
     if ( associated(sect%next_sect) )     deallocate(sect%next_sect)
     if ( associated(sect%next_sect_t) )   deallocate(sect%next_sect_t)
  
! deallocate fmat one by one
     if ( associated(sect%fmat) ) then
         do i=1,sect%nops
             do j=0,1
                 call dealloc_one_fmat(sect%fmat(i,j))
             enddo ! over j={0,1} loop
         enddo ! over i={1,sect%nops} loop
         deallocate(sect%fmat)
     endif ! back if ( associated(sect%fmat) ) block
  
     return
  end subroutine dealloc_one_sect

!!>>> ctqmc_allocate_memory_sect: allocate memory for sectors related variables
  subroutine ctqmc_allocate_memory_sect()
     implicit none

! local variables
     integer :: i

! allocate memory
     allocate( sectors(nsect),      stat=istat )
     allocate( is_trunc(nsect),     stat=istat )
     allocate( is_string(nsect,2),  stat=istat )
     allocate( prob_sect(nsect),    stat=istat )

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_sect', 'can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     do i=1,nsect
         sectors(i)%ndim = 0
         sectors(i)%nelec = 0
         sectors(i)%nops = norbs
         sectors(i)%istart = 0
         sectors(i)%eval => null()
         sectors(i)%next_sect => null()
         sectors(i)%next_sect_t => null()
         sectors(i)%fmat => null()
     enddo ! over i={1,nsect} loop

     is_trunc = .false.
     is_string = .false.
     prob_sect = zero

     return
  end subroutine ctqmc_allocate_memory_sect

!!>>> ctqmc_deallocate_memory_sect: deallocate memory for sectors related variables
  subroutine ctqmc_deallocate_memory_sect()
     implicit none

! local variables
     integer :: i

     if ( allocated(sectors) ) then
! first, loop over all the sectors and deallocate their memory
         do i=1,nsect
             call dealloc_one_sect(sectors(i))
         enddo ! over i={1,nsect} loop

! then, deallocate memory of sectors
         deallocate(sectors)
     endif ! back if ( allocated(sectors) ) block

     if ( allocated(is_trunc) )  deallocate(is_trunc)
     if ( allocated(is_string) ) deallocate(is_string)
     if ( allocated(prob_sect) ) deallocate(prob_sect)

     return
  end subroutine ctqmc_deallocate_memory_sect

!!>>> ctqmc_allocate_memory_occu: allocate memory for occu
  subroutine ctqmc_allocate_memory_occu()
     implicit none

     integer :: i,j,k

! allocate them
     allocate( fprod(nsect,2),     stat=istat )
     allocate( occu(norbs, nsect), stat=istat )

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_occu', 'can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     do i=1,nsect
         if ( is_trunc(i) ) cycle
         
         do j=1,2
             fprod(i,j)%n = sectors(i)%ndim
             call alloc_one_sqrmat(fprod(i,j))
         enddo ! over j={1,2} loop

         do j=1,norbs
             occu(j,i)%n = sectors(i)%ndim
             call alloc_one_sqrmat(occu(j,i))
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,nsect} loop
 
     if ( idoub == 2 ) then
         allocate( doccu(norbs,norbs,nsect),  stat = istat )

         if ( istat /= 0 ) then
             call s_print_error('ctqmc_allocate_memory_occu', 'can not allocate enough memory')
         endif ! back if ( istat /= 0 ) block

         do i=1,nsect
             if (is_trunc(i)) cycle

             do j=1,norbs
                 do k=1,norbs
                     doccu(k,j,i)%n = sectors(i)%ndim
                     call alloc_one_sqrmat(doccu(k,j,i))
                 enddo ! over k={1,norbs} loop
             enddo ! over j={1,norbs} loop
         enddo ! over i={1,nsect} loop
     endif ! back if ( idoub == 2 ) block
  
     return
  end subroutine ctqmc_allocate_memory_occu

!!>>> ctqmc_deallocate_memory_occu: deallocate memory for occu
  subroutine ctqmc_deallocate_memory_occu()
     implicit none

     integer :: i,j,k

     if ( allocated(fprod) ) then
         do i=1,nsect
             if ( is_trunc(i) ) cycle
             do j=1,2
                 call dealloc_one_sqrmat(fprod(i,j))
             enddo ! over j={1,2} loop
         enddo ! over i={1,nsect} loop
         deallocate(fprod)
     endif ! back if ( allocated(fprod) ) block

     if ( allocated(occu) ) then
         do i=1,nsect
             if ( is_trunc(i) ) cycle
             do j=1,norbs
                 call dealloc_one_sqrmat(occu(j,i))
             enddo ! over j={1,norbs} loop
         enddo ! over i={1,nsect} loop
         deallocate(occu)
     endif ! back if ( allocated(occu) ) block

     if ( allocated(doccu) ) then
         do i=1,nsect
             if ( is_trunc(i) ) cycle
             do j=1,norbs
                 do k=1,norbs
                     call dealloc_one_sqrmat(doccu(k,j,i)) 
                 enddo ! over k={1,norbs} loop
             enddo ! over j={1,norbs} loop
         enddo ! over i={1,nsect} loop
         deallocate(doccu)
     endif ! back if ( allocated(doccu) ) block

     return
  end subroutine ctqmc_deallocate_memory_occu

!!>>> ctqmc_make_trunc: subroutine used to truncate the Hilbert space
  subroutine ctqmc_make_trunc()
     implicit none

! local variables
! loop index
     integer :: i

! temp integer
     integer :: i1

! file status
     logical :: exists

! don't truncate the Hilbert space at all
     if ( itrun == 1 ) then
         do i=1,nsect
             sectors(i)%next_sect_t = sectors(i)%next_sect
         enddo ! over i={1,nsect} loop
         nsect_t = nsect
         mdim_sect_t = mdim_sect
         adim_sect_t = adim_sect

         if ( myid == master ) then
             write(mystd,*)
             write(mystd,'(4X,a)') 'use full Hilbert space' 
         endif ! back if ( myid == master ) block

! truncate the Hilbert space according to the total occupancy number
     elseif ( itrun == 2 ) then
         is_trunc = .false.
         do i=1,nsect
             if ( sectors(i)%nelec < nmini .or. sectors(i)%nelec > nmaxi ) then
                 is_trunc(i) = .true.
             endif
         enddo ! over i={1,nsect} loop

         call cat_trunc_sect()

         if ( myid == master ) then
             write(mystd,*)
             write(mystd,'(4X,a,i2,a,i2)') 'truncate occupancy number, just keep ', &
                                               nmini,' ~ ',  nmaxi 
         endif ! back if ( myid == master ) block

! truncate the Hilbert space according to the total occupancy number and 
! the probatility of atomic states
     elseif ( itrun == 3 ) then
         if ( myid == master ) then
             write(mystd,*)
             write(mystd,'(4X,a,i2,a,i2)') 'truncate occupancy number, just keep ', &
                                               nmini, ' ~ ',  nmaxi
             is_trunc = .false.
             prob_sect = zero
             inquire(file='solver.psect.dat', exist=exists)
             if ( exists ) then
                 write(mystd,'(4X,a)') 'truncate high energy atomic states according to their probability'

                 open(mytmp, file='solver.psect.dat', form='formatted', status='unknown')
                 read(mytmp, *) ! skip header
                 do i=1,nsect
                     read(mytmp,*) i1, prob_sect(i) 
                 enddo ! over i={1,nsect} loop
                 close(mytmp)

                 do i=1,nsect
                     if ( sectors(i)%nelec < nmini .or. &
                          sectors(i)%nelec > nmaxi .or. &
                          prob_sect(i) < 1e-4           ) then

                         is_trunc(i) = .true.
                     endif
                 enddo ! over i={1,nsect} loop
             else
                 do i=1,nsect
                     if ( sectors(i)%nelec < nmini .or. sectors(i)%nelec > nmaxi ) then
                         is_trunc(i) = .true.
                     endif
                 enddo ! over i={1,nsect} loop
             endif ! back if ( exists ) block
         endif ! back if ( myid == master ) block

# if defined (MPI)
! block until all processes have reached here
         call mp_barrier()
         call mp_bcast(prob_sect, master)
         call mp_bcast(is_trunc, master)
         call mp_barrier()
# endif  /* MPI */

         call cat_trunc_sect()
     endif ! back if ( itrun == 1 ) block

! print summary of sectors after truncating
     if ( myid == master ) then
         if ( itrun == 1 ) then 
             write(mystd,'(4X,a,i5)')    'tot_num_sect:', nsect_t
             write(mystd,'(4X,a,i5)')    'max_dim_sect:', mdim_sect_t
             write(mystd,'(4X,a,f10.2)') 'ave_dim_sect:', adim_sect_t
             write(mystd,*)
         elseif ( itrun == 2 .or. itrun == 3 ) then
              write(mystd,'(4X,a)') 'before truncation:'
              write(mystd,'(4X,a,i5)')    'tot_num_sect: ', nsect
              write(mystd,'(4X,a,i5)')    'max_dim_sect: ', mdim_sect
              write(mystd,'(4X,a,f10.2)') 'ave_dim_sect: ', adim_sect
              write(mystd,'(4X,a)') 'after truncation:'
              write(mystd,'(4X,a,i5)')    'tot_num_sect: ', nsect_t
              write(mystd,'(4X,a,i5)')    'max_dim_sect: ', mdim_sect_t
              write(mystd,'(4X,a,f10.2)') 'ave_dim_sect: ', adim_sect_t
              write(mystd,*)
         endif ! back if ( itrun == 1 ) block
     endif ! back if ( myid == master ) block

     return
  end subroutine ctqmc_make_trunc

!!>>> cat_trunc_sect: recalculate information of sectors after truncation
  subroutine cat_trunc_sect()
     implicit none

     integer :: i,j,k,ii
     integer :: sum_dim
 
     mdim_sect_t = -1
     nsect_t = 0
     sum_dim = 0
     do i=1,nsect
         sectors(i)%next_sect_t = -1
         if ( is_trunc(i) ) then
             cycle
         endif ! back if ( is_trunc(i) ) block
         if ( mdim_sect_t < sectors(i)%ndim ) then
             mdim_sect_t = sectors(i)%ndim
         endif ! back if ( mdim_sect_t < sectors(i)%ndim ) block 

         sum_dim = sum_dim + sectors(i)%ndim
         nsect_t = nsect_t + 1 
         do j=1,sectors(i)%nops
             do k=0,1
                 ii = sectors(i)%next_sect(j,k) 
                 if ( ii == -1 ) cycle
                 if ( .not. is_trunc(ii) ) then
                     sectors(i)%next_sect_t(j,k) = ii
                 endif ! back if ( .not. is_trunc(ii) ) block
             enddo ! over k={0,1} loop
         enddo ! over j={1,sectors(i)%nops} loop
     enddo ! over i={1,nsect} loop

     adim_sect_t = real(sum_dim) / real(nsect_t)

     return
  end subroutine cat_trunc_sect

!!>>> ctqmc_make_occu: calculate < c^{\dag} c>, < c^{\dag} c c^{\dag} c >
  subroutine ctqmc_make_occu()
     implicit none

     integer :: i,j,k,ii,jj
     real(dp) :: t_mat1(mdim_sect_t, mdim_sect_t) 
     real(dp) :: t_mat2(mdim_sect_t, mdim_sect_t) 

     do i=1,norbs
         do j=1,nsect
             if ( is_trunc(j) ) cycle
             k=sectors(j)%next_sect(i,0)
             if ( k == -1 ) then
                 occu(i,j)%item = zero
                 cycle
             endif ! back if ( k == -1 ) block

             call dgemm( 'N', 'N', sectors(j)%ndim, sectors(j)%ndim, sectors(k)%ndim, &
                         one,  sectors(k)%fmat(i,1)%item,            sectors(j)%ndim, &
                               sectors(j)%fmat(i,0)%item,            sectors(k)%ndim, & 
                         zero, occu(i,j)%item,                       sectors(j)%ndim  ) 
         enddo ! over j={1,nsect} loop
     enddo ! over i={1,norbs} loop

     if ( idoub == 2 ) then
         do i=1,norbs
             do j=1,norbs
                 do k=1,nsect
                     if ( is_trunc(k) ) cycle 
                     jj = sectors(k)%next_sect(j,0) 
                     ii = sectors(k)%next_sect(i,0)
                     if ( ii == -1 .or. jj == -1 ) then
                         doccu(i,j,k)%item = zero
                         cycle
                     endif ! back if ( ii == -1 .or. jj == -1 ) block

                     call dgemm( 'N', 'N', sectors(k)%ndim, sectors(k)%ndim, sectors(jj)%ndim, &
                                 one,  sectors(jj)%fmat(j,1)%item,           sectors(k)%ndim,  & 
                                       sectors(k)%fmat(j,0)%item,            sectors(jj)%ndim, & 
                                 zero, t_mat1,                               mdim_sect_t       ) 

                     call dgemm( 'N', 'N', sectors(k)%ndim, sectors(k)%ndim, sectors(ii)%ndim, &
                                 one,  sectors(ii)%fmat(i,1)%item,           sectors(k)%ndim,  &
                                       sectors(k)%fmat(i,0)%item,            sectors(ii)%ndim, & 
                                 zero, t_mat2,                               mdim_sect_t       ) 

                     call dgemm( 'N', 'N', sectors(k)%ndim, sectors(k)%ndim, sectors(k)%ndim,  &
                                 one,  t_mat2,                               mdim_sect_t,      & 
                                       t_mat1,                               mdim_sect_t,      & 
                                 zero, doccu(i,j,k)%item,                    sectors(k)%ndim   )

                 enddo ! over k={1,nsect} loop
             enddo ! over j={1,norbs} loop
         enddo ! over i={1,norbs} loop
     endif ! back if ( idoub == 2 ) block

     return
  end subroutine ctqmc_make_occu

!!>>> ctqmc_make_string: subroutine used to build a string
  subroutine ctqmc_make_string(csize, index_t_loc, string)
     implicit none

! external variables
! number of fermion operators for current diagram
     integer, intent(in) :: csize

! address index of fermion operators
     integer, intent(in) :: index_t_loc(mkink)

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

     is_string(:,1) = .true.
     string = -1

! we build a string from right to left, that is, beta <------- 0
! begin with S1: F1(S1)-->S2, F2(S2)-->S3, ... ,Fk(Sk)-->S1
! if find some Si==-1, cycle this sector immediately
     do i=1,nsect
         if ( is_trunc(i) ) then
             is_string(i,1) = .false.
         endif ! back if ( is_trunc(i) ) block

         curr_sect_l = i
         curr_sect_r = i
         next_sect_l = i
         next_sect_r = i
         left = 0
         right = csize + 1

         do j=1,csize
             if ( mod(j,2) == 1) then
                 left = left + 1
                 string(left,i) = curr_sect_l
                 vt = type_v( index_t_loc(left) )
                 vf = flvr_v( index_t_loc(left) ) 
                 next_sect_l = sectors(curr_sect_l)%next_sect_t(vf,vt)
                 if ( next_sect_l == -1 ) then
                     is_string(i,1) = .false. 
                     EXIT   ! finish check, exit
                 endif ! back if ( next_sect_l == -1 ) block
                 curr_sect_l = next_sect_l
             else
                 right = right - 1
                 vt = type_v( index_t_loc(right) )
                 vf = flvr_v( index_t_loc(right) ) 
                 vt = mod(vt+1,2)
                 next_sect_r = sectors(curr_sect_r)%next_sect_t(vf,vt)
                 if ( next_sect_r == -1 ) then
                     is_string(i,1) = .false. 
                     EXIT   ! finish check, exit
                 endif ! back if ( next_sect_r == -1 ) block
                 string(right,i) = next_sect_r
                 curr_sect_r = next_sect_r
             endif ! back if ( mod(j,2) == 1 ) block
         enddo ! over j={1,csize} loop

! if it doesn't form a string, we cycle it, go to the next sector
         if ( .not. is_string(i,1) ) then
             cycle
         endif ! back if ( .not. is_string(i,1) ) block

! add the last sector to string, and check whether string(csize+1,i) == string(1,i)
! important for csize = 0
         string(csize+1,i) = i

! this case will generate a non-diagonal block, it will not contribute to trace 
         if ( next_sect_r /= next_sect_l ) then
             is_string(i,1) = .false.
         endif ! back if ( next_sect_r /= next_sect_l ) block

     enddo ! over i={1,nsect} loop

     return
  end subroutine ctqmc_make_string

  end module m_sector
