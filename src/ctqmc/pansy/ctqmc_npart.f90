!!!-------------------------------------------------------------------------
!!! project : pansy
!!! program : m_npart module
!!!           m_npart@ctqmc_allocate_memory_part
!!!           m_npart@ctqmc_deallocate_memory_part
!!!           m_npart@ctqmc_make_npart
!!!           m_npart@cat_sector_ztrace
!!!           m_npart@ctqmc_save_npart
!!! source  : ctqmc_npart.f90
!!! type    : module
!!! authors : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014  by yilin wang
!!!           07/19/2014  by yilin wang
!!!           08/09/2014  by yilin wang
!!!           08/18/2014  by yilin wang
!!!           11/02/2014  by yilin wang
!!! purpose : define data structure for divide and conquer (npart) algorithm
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> m_npart: contains some key global variables and subroutines for divide
!!>>> and conquer (npart) algorithm to speed up the trace evaluation
  module m_npart
     use constants, only : dp, zero, one
     use control, only : npart, mkink, beta, ncfgs
     use context, only : time_v, expt_v, type_v, flvr_v

     use m_sect, only : nsect, sectors, mdim_sect
 
     implicit none

! status flag for allocating memory
     integer, private :: istat

! total number of matrices products
     real(dp), public, save :: nprod = zero

! the first filled part 
     integer, public, save :: ffpart = 0

! how to treat each part when calculating trace
     integer, public, save, allocatable :: isave(:,:,:)

! whether to copy this part ?
     logical, public, save, allocatable :: is_cp(:,:)

! number of columns to be copied, in order to save copy time 
     integer, public, save, allocatable :: ncol_cp(:,:)
 
! the start positions of fermion operators for each part
     integer, public, save, allocatable :: ops(:)

! the end positions of fermion operators for each part
     integer, public, save, allocatable :: ope(:)

! saved parts of matrices product, for previous accepted configuration 
     real(dp), public, save, allocatable :: saved_p(:,:,:,:)

! saved parts of matrices product, for new proposed configuration
     real(dp), public, save, allocatable :: saved_n(:,:,:,:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================
    
     public :: ctqmc_allocate_memory_part
     public :: ctqmc_deallocate_memory_part
     public :: ctqmc_make_npart
     public :: cat_sector_ztrace
     public :: ctqmc_save_npart

  contains ! encapsulated functionality

!!>>> ctqmc_allocate_memory_part: allocate memory for npart related variables
  subroutine ctqmc_allocate_memory_part()
     implicit none

! allocate memory
     allocate( isave(npart, nsect, 2), stat=istat )
     allocate( is_cp(npart, nsect),    stat=istat )
     allocate( ncol_cp(npart, nsect),  stat=istat )
     allocate( ops(npart),             stat=istat )
     allocate( ope(npart),             stat=istat )

     allocate( saved_p(mdim_sect, mdim_sect, npart, nsect), stat=istat )
     allocate( saved_n(mdim_sect, mdim_sect, npart, nsect), stat=istat )

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_sect', &
                          'can not allocate enough memory')
     endif

! initialize them
     isave = 1
     is_cp = .false.
     ncol_cp = 0
     ops = 0
     ope = 0
     saved_p = zero
     saved_n = zero

     return
  end subroutine ctqmc_allocate_memory_part

!!>>> ctqmc_deallocate_memory_part: deallocate memory for 
!!>>> npart-related variables
  subroutine ctqmc_deallocate_memory_part()
     implicit none

     if ( allocated(isave) )    deallocate(isave)
     if ( allocated(is_cp) )    deallocate(is_cp)
     if ( allocated(ncol_cp) )  deallocate(ncol_cp)
     if ( allocated(ops) )      deallocate(ops)
     if ( allocated(ope) )      deallocate(ope)
     if ( allocated(saved_p) )  deallocate(saved_p)
     if ( allocated(saved_n) )  deallocate(saved_n)
     
     return
  end subroutine ctqmc_deallocate_memory_part

!!>>> ctqmc_make_npart: subroutine used to determine isave 
  subroutine ctqmc_make_npart(cmode, csize, index_t_loc, tau_s, tau_e)
     implicit none
   
! external arguments
! mode for different Monte Carlo moves
     integer,  intent(in)  :: cmode
   
! total number of operators for current diagram
     integer,  intent(in)  :: csize
   
! local version of index_t
     integer, intent(in) :: index_t_loc(mkink)
   
! imaginary time value of operator A, only valid in cmode = 1 or 2
     real(dp), intent(in) :: tau_s
   
! imaginary time value of operator B, only valid in cmode = 1 or 2
     real(dp), intent(in) :: tau_e
        
! local variables
! length of imaginary time axis for each part
     real(dp) :: interval
   
! number of fermion operators for each part
     integer :: nop(npart)
   
! position of the operator A and operator B, index of part
     integer  :: tis
     integer  :: tie
     integer  :: tip
   
! loop index
     integer :: i, j
   
! init key arrays
     nop = 0
     ops = 0
     ope = 0
     ffpart = 0
 
! isave: how to treat each part for each alive string
! isave = 0: matrices product for this part has been calculated previously
! isave = 1: this part should be recalculated, and the result must be
!            stored in saved_p, if this Monte Caro move has been accepted.
! isave = 2: this part is empty, we don't need to do anything with them.

! copy isave 
     isave(:,:,1) = isave(:,:,2)
   
! case 1: recalculate all the matrices products
     if ( npart == 1 ) then
         nop(1) = csize
         ops(1) = 1
         ope(1) = csize
         ffpart = 1
         if (nop(1) <= 0) then
             isave(1,:,1) = 2
         else
             isave(1,:,1) = 1
         endif
! case 2: use npart alogithm
     elseif ( npart > 1) then
         interval = beta / real(npart)
! calculate number of operators for each part
         do i=1,csize
             j = ceiling( time_v( index_t_loc(i) ) / interval )
             nop(j) = nop(j) + 1
         enddo 
! if no operators in this part, ignore them
         do i=1, npart
             if (ffpart == 0 .and. nop(i) > 0) then
                 ffpart = i
             endif
             if (nop(i) <= 0) then
                 isave(i,:,1) = 2 
             endif
         enddo
! calculate the start and end index of operators for each part
         do i=1,npart
             if ( nop(i) > 0 ) then
                 ops(i) = 1
                 do j=1,i-1
                     ops(i) = ops(i) + nop(j)
                 enddo 
                 ope(i) = ops(i) + nop(i) - 1
             endif 
         enddo 
   
! case 2A: use some saved matrices products from previous accepted Monte Carlo move
         if (cmode == 1 .or. cmode == 2) then
! get the position of operator A and operator B
             tis = ceiling( tau_s / interval )
             tie = ceiling( tau_e / interval )
! operator A:
             if (nop(tis)>0) then
                 isave(tis,:,1) = 1
             endif
! special attention: if operator A is on the left or right boundary, then
! the neighbour part should be recalculated as well
             if ( nop(tis) > 0 ) then
                 if ( tau_s >= time_v( index_t_loc( ope(tis) ) ) ) then
                     tip = tis + 1
                     do while ( tip <= npart )
                         if ( nop(tip) > 0 ) then
                             isave(tip,:,1) = 1;  EXIT
                         endif
                         tip = tip + 1
                     enddo ! over do while loop
                 endif
! for remove an operator, nop(tis) may be zero
             else
                 tip = tis + 1
                 do while ( tip <= npart )
                     if ( nop(tip) > 0 ) then
                         isave(tip,:,1) = 1; EXIT
                     endif
                     tip = tip + 1
                 enddo ! over do while loop
             endif ! back if ( nop(tis) > 0 ) block
   
! operator B:
             if (nop(tie)>0) then
                 isave(tie,:,1) = 1
             endif
! special attention: if operator B is on the left or right boundary, then
! the neighbour part should be recalculated as well
             if ( nop(tie) > 0 ) then
                 if ( tau_e >= time_v( index_t_loc( ope(tie) ) ) ) then
                     tip = tie + 1
                     do while ( tip <= npart )
                         if ( nop(tip) > 0 ) then
                             isave(tip,:,1) = 1; EXIT
                         endif
                         tip = tip + 1
                     enddo ! over do while loop
                 endif
! for remove an operator, nop(tie) may be zero
             else
                 tip = tie + 1
                 do while ( tip <= npart )
                     if ( nop(tip) > 0 ) then
                         isave(tip,:,1) = 1; EXIT
                     endif
                     tip = tip + 1
                 enddo ! over do while loop
             endif ! back if ( nop(tie) > 0 ) block
   
! case 2B: recalculate all the matrices products 
         elseif (cmode == 3 .or. cmode == 4) then
             do i=1, nsect
                 do j=1, npart
                     if (isave(j,i,1) == 0) then
                         isave(j,i,1) = 1
                     endif  
                 enddo
             enddo
         endif ! back if (cmode == 1 .or. cmode == 2) block
   
! npart should be larger than zero
     else
         call s_print_error('ctqmc_make_ztrace', 'npart is small than 1, &
                                    it should be larger than zero')
     endif ! back if ( npart == 1 ) block
   
     return
  end subroutine ctqmc_make_npart
   
!!>>> cat_sector_ztrace: calculate the trace for one sector
  subroutine cat_sector_ztrace(csize, string, index_t_loc, expt_t_loc, trace)
     implicit none
   
! external variables
! number of total fermion operators
     integer, intent(in) :: csize
   
! string for this sector
     integer, intent(in) :: string(csize+1)
   
! address index of fermion operators
     integer, intent(in) :: index_t_loc(mkink)
   
! diagonal elements of last time-evolution matrices
     real(dp), intent(in) :: expt_t_loc(ncfgs)
   
! the calculated trace of this sector
     real(dp), intent(out) :: trace
   
! local variables
! temp matrices
     real(dp) :: r_mat(mdim_sect, mdim_sect)
     real(dp) :: t_mat(mdim_sect, mdim_sect)
   
! temp index
     integer :: dim1, dim2, dim3, dim4
     integer :: isect, sect1, sect2
     integer :: indx
     integer :: counter
     integer :: vt, vf
   
! loop index
     integer :: i, j, k, l
   
     isect = string(1)

!--------------------------------------------------------------------
! from right to left: beta <------ 0
     dim1 = sectors(string(1))%ndim
     r_mat = zero 
     t_mat = zero

! loop over all the parts
     do i=1, npart
   
! this part has been calculated previously, just use its results
         if ( isave(i, isect, 1) == 0 ) then
             sect1 = string(ope(i)+1)
             sect2 = string(ops(i))
             dim2 = sectors(sect1)%ndim
             dim3 = sectors(sect2)%ndim
             if ( i > ffpart ) then
                 call dgemm( 'N', 'N', dim2, dim1, dim3, one,     &
                              saved_p(:, :, i, isect), mdim_sect, &
                              r_mat,                   mdim_sect, &
                              zero, t_mat,             mdim_sect  )

                 r_mat(:, 1:dim1) = t_mat(:, 1:dim1)
                 nprod = nprod + one
             else
                 r_mat(:, 1:dim1) = saved_p(:, 1:dim1, i, isect)
             endif ! back if ( i > ffpart ) block
  
! this part should be recalcuated 
         elseif ( isave(i, isect, 1) == 1 ) then 
             sect1 = string(ope(i)+1)
             sect2 = string(ops(i))
             dim4 = sectors(sect2)%ndim
             saved_n(:, :, i, isect) = zero
   
! loop over all the fermion operators in this part
             counter = 0
             do j=ops(i), ope(i)
                 counter = counter + 1
                 indx = sectors(string(j)  )%istart
                 dim2 = sectors(string(j+1))%ndim
                 dim3 = sectors(string(j)  )%ndim
   
! multiply the diagonal matrix of time evolution operator 
                 if (counter > 1) then
                     do l=1,dim4
                         do k=1,dim3
                             t_mat(k, l) = saved_n(k, l, i, isect) * expt_v(indx+k-1, index_t_loc(j))
                         enddo
                     enddo
                     nprod = nprod + one
                 else
                     t_mat = zero
                     do k=1,dim3
                         t_mat(k,k) = expt_v(indx+k-1, index_t_loc(j))
                     enddo
                 endif
   
! multiply the matrix of fermion operator
                 vt = type_v( index_t_loc(j) )
                 vf = flvr_v( index_t_loc(j) ) 
                 call dgemm( 'N', 'N', dim2, dim4, dim3, one,            &
                             sectors(string(j))%fmat(vf, vt)%item, dim2, &
                             t_mat,                           mdim_sect, &
                             zero, saved_n(:, :, i, isect),   mdim_sect  ) 
   
                 nprod = nprod + one
             enddo  

! set its save status and copy status
             isave(i, isect, 1) = 0
             is_cp(i, isect) = .true.
             ncol_cp(i, isect) = dim4
   
! multiply this part with the rest parts
             if ( i > ffpart ) then
                 call dgemm( 'N', 'N', dim2, dim1, dim4, one,    &
                             saved_n(:, :, i, isect), mdim_sect, &
                             r_mat,                   mdim_sect, &
                             zero, t_mat,             mdim_sect  ) 

                 r_mat(:, 1:dim1) = t_mat(:, 1:dim1)
                 nprod = nprod + one
             else
                 r_mat(:, 1:dim1) = saved_n(:, 1:dim1, i, isect)
             endif

! no operators in this part, do nothing
         elseif ( isave(i, isect, 1) == 2 ) then
             cycle
         endif ! back if ( is_save(i, isect) ==0 )  block
   
! the start sector for next part
         isect = sect1
   
     enddo  ! over i={1, npart} loop   
   
! special treatment of the last time-evolution operator
     indx = sectors(string(1))%istart

! no fermion operators
     if ( csize == 0 ) then
         do k=1, dim1
             r_mat(k, k) = expt_t_loc(indx+k-1)
         enddo
! multiply the last time-evolution operator
     else
         do l=1,dim1
             do k=1,dim1
                 r_mat(k, l) = r_mat(k, l) * expt_t_loc(indx+k-1)
             enddo
         enddo
         nprod = nprod + one
     endif
   
! store final product
     sectors( string(1) )%fprod(:, :, 1) = r_mat(1:dim1, 1:dim1)
   
! calculate the trace
     trace  = zero
     do j=1, sectors(string(1))%ndim
         trace = trace + r_mat(j,j)
     enddo
   
     return
  end subroutine cat_sector_ztrace

!!>>> ctqmc_save_npart: copy data if propose has been accepted
  subroutine ctqmc_save_npart()
     implicit none

! loop index
     integer :: i,j

! copy save-state for all the parts 
     isave(:, :, 2) = isave(:, :, 1)

! when npart > 1, we used the npart algorithm, save the changed 
! matrices products when moves are accepted
     if ( npart > 1 ) then
         do i=1, nsect
             do j=1, npart
                 if ( is_cp(j, i) ) then
                     saved_p(:, 1:ncol_cp(j, i), j, i) = saved_n(:, 1:ncol_cp(j, i), j, i) 
                 endif
             enddo
         enddo
     endif

     return
  end subroutine ctqmc_save_npart

  end module m_npart
