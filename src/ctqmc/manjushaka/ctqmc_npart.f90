!!!-------------------------------------------------------------------------
!!! project : manjushaka
!!! program : m_npart  module
!!!           m_npart@ctqmc_allocate_memory_part
!!!           m_npart@ctqmc_deallocate_memory_part
!!!           m_npart@ctqmc_make_npart
!!!           m_npart@cat_sector_ztrace
!!!           m_npart@ctqmc_save_npart
!!! source  : ctqmc_npart.f90
!!! type    : module
!!! authors : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           07/19/2014 by yilin wang
!!!           08/09/2014 by yilin wang
!!!           08/13/2014 by yilin wang
!!!           08/20/2014 by yilin wang
!!!           11/02/2014 by yilin wang
!!! purpose : define data structure for divide conquer (npart) algorithm
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> containing the information for npart trace algorithm
  module m_npart
     use constants, only : dp, zero, one
     use control, only : npart, mkink, ncfgs, beta
     use context, only : time_v, type_v, flvr_v, expt_v

     use m_sector, only : nsect, sectors, is_trunc, t_matrix
     use m_sector, only : mdim_sect_t, fprod
     use m_sector, only : alloc_one_mat, dealloc_one_mat

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

! saved parts of matrices product, for previous configuration
     type(t_matrix), public, save, allocatable :: saved_p(:,:)

! saved parts of matrices product, for new configuration
     type(t_matrix), public, save, allocatable :: saved_n(:,:)

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: ctqmc_allocate_memory_part
     public :: ctqmc_deallocate_memory_part
     public :: ctqmc_make_npart
     public :: cat_sector_ztrace
     public :: ctqmc_save_npart

  contains ! encapsulated functionality

!!>>> ctqmc_allocate_memory_part: allocate memory for sectors related variables
  subroutine ctqmc_allocate_memory_part()
     implicit none

     integer :: i,j

! allocate memory
     allocate( isave(npart, nsect, 2), stat=istat )
     allocate( is_cp(npart, nsect),    stat=istat )
     allocate( ncol_cp(npart, nsect),  stat=istat )
     allocate( ops(npart),             stat=istat )
     allocate( ope(npart),             stat=istat )
     allocate( saved_p(npart, nsect),  stat=istat )
     allocate( saved_n(npart, nsect),  stat=istat )

! check the status
     if ( istat /= 0 ) then
         call s_print_error('ctqmc_allocate_memory_sect', 'can not allocate enough memory')
     endif

! initialize them
     do i=1,nsect
         if ( is_trunc(i) ) cycle
         do j=1,npart
             saved_p(j,i)%n = mdim_sect_t
             saved_p(j,i)%m = mdim_sect_t
             saved_n(j,i)%n = mdim_sect_t
             saved_n(j,i)%m = mdim_sect_t
             call alloc_one_mat(saved_p(j,i))
             call alloc_one_mat(saved_n(j,i))
         enddo ! over j={1,npart} loop
     enddo ! over i={1,nsect} loop

     isave = 1
     is_cp = .false.
     ncol_cp = 0
     ops = 0
     ope = 0

     return
  end subroutine ctqmc_allocate_memory_part

!!>>> ctqmc_deallocate_memory_part: deallocate memory for sectors related variables
  subroutine ctqmc_deallocate_memory_part()
     implicit none

     integer :: i,j

     if ( allocated(isave) )     deallocate(isave)
     if ( allocated(is_cp) )     deallocate(is_cp)
     if ( allocated(ncol_cp) )   deallocate(ncol_cp)
     if ( allocated(ops) )       deallocate(ops)
     if ( allocated(ope) )       deallocate(ope)

     if ( allocated(saved_p) ) then
         do i=1,nsect
             if (is_trunc(i)) cycle
             do j=1,npart
                 call dealloc_one_mat(saved_p(j,i))
             enddo ! over j={1,npart} loop
         enddo ! over i={1,nsect} loop
         deallocate(saved_p)
     endif ! back if ( allocated(saved_p) ) block

     if ( allocated(saved_n) ) then
         do i=1,nsect
             if ( is_trunc(i) ) cycle
             do j=1,npart
                 call dealloc_one_mat(saved_n(j,i))
             enddo ! over j={1,npart} loop
         enddo ! over i={1,nsect} loop
         deallocate(saved_n)
     endif ! back if ( allocated(saved_n) ) block

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
         if ( nop(1) <= 0 ) then
             isave(1,:,1) = 2
         else
             isave(1,:,1) = 1
         endif ! back if ( nop(1) <= 0 ) block
! case 2: use npart alogithm
     elseif ( npart > 1 ) then
         interval = beta / real(npart)
! calculate number of operators for each part
         do i=1,csize
             j = ceiling( time_v( index_t_loc(i) ) / interval )
             nop(j) = nop(j) + 1
         enddo  ! over i={1,csize} loop
! if no operators in this part, ignore them
         do i=1,npart
             if ( ffpart == 0 .and. nop(i) > 0 ) then
                 ffpart = i
             endif ! back if ( ffpart == 0 .and. nop(i) > 0 ) block
             if ( nop(i) <= 0 ) then
                 isave(i,:,1) = 2
             endif ! back if ( nop(i) <= 0 ) block
         enddo ! over i={1,npart} loop
! calculate the start and end index of operators for each part
         do i=1,npart
             if ( nop(i) > 0 ) then
                 ops(i) = 1
                 do j=1,i-1
                     ops(i) = ops(i) + nop(j)
                 enddo  ! over j={1,i-1} loop
                 ope(i) = ops(i) + nop(i) - 1
             endif  ! back if ( nop(i) > 0 ) block
         enddo  ! over i={1,npart} loop

! case 2A: use some saved matrices products from previous accepted Monte Carlo move
         if ( cmode == 1 .or. cmode == 2 ) then
! get the position of operator A and operator B
             tis = ceiling( tau_s / interval )
             tie = ceiling( tau_e / interval )
! operator A:
             if ( nop(tis) > 0 ) then
                 isave(tis,:,1) = 1
             endif ! back if ( nop(tis) > 0 ) block
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
                     enddo ! over do while ( tip <= npart ) loop
                 endif ! back if ( tau_s >= time_v( index_t_loc( ope(tis) ) ) ) block
! for remove an operator, nop(tis) may be zero
             else
                 tip = tis + 1
                 do while ( tip <= npart )
                     if ( nop(tip) > 0 ) then
                         isave(tip,:,1) = 1; EXIT
                     endif ! back if ( nop(tip) > 0 ) block
                     tip = tip + 1
                 enddo ! over do while ( tip <= npart ) loop
             endif ! back if ( nop(tis) > 0 ) block

! operator B:
             if ( nop(tie) > 0 ) then
                 isave(tie,:,1) = 1
             endif ! back if ( nop(tie) > 0 ) block
! special attention: if operator B is on the left or right boundary, then
! the neighbour part should be recalculated as well
             if ( nop(tie) > 0 ) then
                 if ( tau_e >= time_v( index_t_loc( ope(tie) ) ) ) then
                     tip = tie + 1
                     do while ( tip <= npart )
                         if ( nop(tip) > 0 ) then
                             isave(tip,:,1) = 1; EXIT
                         endif ! back if ( nop(tip) > 0 ) block
                         tip = tip + 1
                     enddo ! over do while ( tau_e >= time_v( index_t_loc( ope(tie) ) ) ) loop
                 endif ! back if ( nop(tie) > 0 ) block
! for remove an operator, nop(tie) may be zero
             else
                 tip = tie + 1
                 do while ( tip <= npart )
                     if ( nop(tip) > 0 ) then
                         isave(tip,:,1) = 1; EXIT
                     endif ! back if ( nop(tip) > 0 ) block
                     tip = tip + 1
                 enddo ! over do while ( tip <= npart ) loop
             endif ! back if ( nop(tie) > 0 ) block

! case 2B: recalculate all the matrices products
         elseif ( cmode == 3 .or. cmode == 4 ) then
             do i=1,nsect
                 do j=1,npart
                     if ( isave(j,i,1) == 0 ) then
                         isave(j,i,1) = 1
                     endif ! back if ( isave(j,i,1) == 0 ) block
                 enddo ! over j={1,npart} loop
             enddo ! over i={1,nsect} loop
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

! index of string for this sector
     integer, intent(in) :: string(csize+1)

! address index of fermion operators
     integer, intent(in) :: index_t_loc(mkink)

! diagonal elements of last time-evolution matrices
     real(dp), intent(in) :: expt_t_loc(ncfgs)

! the calculated trace of this sector
     real(dp), intent(out) :: trace

! local variables
! temp matrices
     real(dp) :: mat_r(mdim_sect_t, mdim_sect_t)
     real(dp) :: mat_t(mdim_sect_t, mdim_sect_t)

! temp index
     integer :: dim1
     integer :: dim2
     integer :: dim3
     integer :: dim4
     integer :: isect
     integer :: sect1
     integer :: sect2
     integer :: indx
     integer :: counter
     integer :: vt
     integer :: vf

! loop index
     integer :: i,j,k,l

     isect = string(1)

!--------------------------------------------------------------------
! from right to left: beta <------ 0
     dim1 = sectors(string(1))%ndim
     mat_r = zero
     mat_t = zero

! loop over all the parts
     do i=1,npart
! this part has been calculated previously, just use its results
         if ( isave(i,isect,1) == 0 ) then
             sect1 = string(ope(i)+1)
             sect2 = string(ops(i))
             dim2 = sectors(sect1)%ndim
             dim3 = sectors(sect2)%ndim

             if ( i > ffpart ) then
                 call dgemm( 'N', 'N', dim2, dim1, dim3,                &
                              one,  saved_p(i,isect)%item, mdim_sect_t, &
                                    mat_r,                 mdim_sect_t, &
                              zero, mat_t,                 mdim_sect_t  )

                 mat_r(:,1:dim1) = mat_t(:,1:dim1)
                 nprod = nprod + one
             else
                 mat_r(:,1:dim1) = saved_p(i,isect)%item(:,1:dim1)
             endif ! back if ( i > ffpart ) block

! this part should be recalcuated
         elseif ( isave(i,isect,1) == 1 ) then
             sect1 = string(ope(i)+1)
             sect2 = string(ops(i))
             dim4 = sectors(sect2)%ndim
             saved_n(i,isect)%item = zero

! loop over all the fermion operators in this part
             counter = 0
             do j=ops(i),ope(i)
                 counter = counter + 1
                 indx = sectors(string(j  ))%istart
                 dim2 = sectors(string(j+1))%ndim
                 dim3 = sectors(string(j  ))%ndim

                 if ( counter > 1 ) then
                     do l=1,dim4
                         do k=1,dim3
                             mat_t(k,l) = saved_n(i,isect)%item(k,l) * expt_v(indx+k-1,index_t_loc(j))
                         enddo ! over k={1,dim3} loop
                     enddo ! over l={1,dim4} loop
                     nprod = nprod + one
                 else
                     mat_t = zero
                     do k=1,dim3
                         mat_t(k,k) = expt_v(indx+k-1,index_t_loc(j))
                     enddo ! over k={1,dim3} loop
                 endif ! back if ( counter > 1 ) block

                 vt = type_v( index_t_loc(j) )
                 vf = flvr_v( index_t_loc(j) )
                 call dgemm( 'N', 'N', dim2, dim4, dim3,                       &
                             one,  sectors(string(j))%fmat(vf, vt)%item, dim2, &
                                   mat_t,                         mdim_sect_t, &
                             zero, saved_n(i,isect)%item,         mdim_sect_t  )

                 nprod = nprod + one
             enddo ! over j={ops(i),ope(i)} loop

             isave(i,isect,1) = 0
             is_cp(i,isect) = .true.
             ncol_cp(i,isect) = dim4

! multiply this part with the rest parts
             if ( i > ffpart ) then
                 call dgemm( 'N', 'N', dim2, dim1, dim4,               &
                             one,  saved_n(i,isect)%item, mdim_sect_t, &
                                   mat_r,                 mdim_sect_t, &
                             zero, mat_t,                 mdim_sect_t  )

                 mat_r(:,1:dim1) = mat_t(:,1:dim1)
                 nprod = nprod + one
             else
                 mat_r(:,1:dim1) = saved_n(i,isect)%item(:,1:dim1)
             endif ! back if ( i > ffpart ) block

         elseif ( isave(i,isect,1) == 2 ) then
             cycle
         endif ! back if ( isave(i,isect) == 0 )  block

! the start sector for next part
         isect = sect1

     enddo ! over i={1, npart} loop

! special treatment of the last time-evolution operator
     indx = sectors(string(1))%istart
     if ( csize == 0 ) then
         do k=1,dim1
             mat_r(k,k) = expt_t_loc(indx+k-1)
         enddo ! over k={1,dim1} loop
     else
         do l=1,dim1
             do k=1,dim1
                 mat_r(k,l) = mat_r(k,l) * expt_t_loc(indx+k-1)
             enddo ! over k={1,dim1} loop
         enddo ! over l={1,dim1} loop
         nprod = nprod + one
     endif ! back if ( csize == 0 ) block

! store final product
     fprod(string(1),1)%item = mat_r(1:dim1,1:dim1)

! calculate the trace
     trace  = zero
     do j=1,sectors(string(1))%ndim
         trace = trace + mat_r(j,j)
     enddo ! over j={1,sectors(string(1))%ndim} loop

     return
  end subroutine cat_sector_ztrace

!!>>> ctqmc_save_npart: copy data when propose has been accepted
  subroutine ctqmc_save_npart()
     implicit none

! loop index
     integer :: i,j

! copy save-state for all the parts
     isave(:,:,2) = isave(:,:,1)

! when npart > 1, we used the npart algorithm, save the changed
! matrices products f moves are accepted
     if ( npart > 1) then
         do i=1,nsect
             if ( is_trunc(i) ) cycle
             do j=1,npart
                 if ( is_cp(j,i) ) then
                     saved_p(j,i)%item(:,1:ncol_cp(j,i)) = saved_n(j,i)%item(:,1:ncol_cp(j,i))
                 endif ! back if ( is_cp(j,i) ) block
             enddo ! over j={1,npart}  loop
         enddo ! over i={1,nsect} loop
     endif ! back if ( npart > 1 ) block

     return
  end subroutine ctqmc_save_npart

  end module m_npart
