!>>> containing the information for npart trace algorithm
  module m_npart
     use constants
     use control
     use context
  
     use m_sector
 
     implicit none

! number of total matrices multiplication
     real(dp), public, save :: num_prod = zero

! the first filled part 
     integer, public, save :: first_fpart = 0

! how to treat each part when calculate trace
     integer, public, save, allocatable :: is_save(:,:,:)

! whether to copy this part
     logical, public, save, allocatable :: is_copy(:,:)

! the number of column to copied
     integer, public, save, allocatable :: col_copy(:,:)
 
! ops, ope
     integer, public, save, allocatable :: ops(:)
     integer, public, save, allocatable :: ope(:)

! saved parts of matrices product, previous configuration 
     type(t_sqrmat), public, save, allocatable :: saved_a(:,:)

! saved parts of matrices product, current configuration
     type(t_sqrmat), public, save, allocatable :: saved_b(:,:)

     contains

!>>> allocate memory for sect-related variables
     subroutine ctqmc_allocate_memory_part()

         implicit none

         integer :: istat
         integer :: i,j

! allocate memory
         
         allocate(is_save(npart, nsectors, 2),    stat=istat)
         allocate(is_copy(npart, nsectors),       stat=istat)
         allocate(col_copy(npart, nsectors),      stat=istat)
         allocate(ops(npart),                     stat=istat)
         allocate(ope(npart),                     stat=istat)
         allocate(saved_a(npart, nsectors),       stat=istat)
         allocate(saved_b(npart, nsectors),       stat=istat)

         do i=1, nsectors
             if (is_trunc(i)) cycle
             do j=1, npart
                 saved_a(j,i)%n = max_dim_sect_trunc
                 saved_b(j,i)%n = max_dim_sect_trunc
                 call alloc_one_sqrmat( saved_a(j,i) )
                 call alloc_one_sqrmat( saved_b(j,i) )
             enddo
         enddo    

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_sect','can not allocate enough memory')
         endif

         is_save = 1
         is_copy = .false.
         col_copy = 0
         ops = 0
         ope = 0

         return
     end subroutine ctqmc_allocate_memory_part

!>>> deallocate memory for sect-related variables
     subroutine ctqmc_deallocate_memory_part()
         implicit none
         
         integer :: i,j

         if ( allocated(is_save) )      deallocate(is_save)
         if ( allocated(is_copy) )      deallocate(is_copy)
         if ( allocated(col_copy) )     deallocate(col_copy)
         if ( allocated(ops) )          deallocate(ops)
         if ( allocated(ope) )          deallocate(ope)

         if ( allocated(saved_a) ) then
             do i=1, nsectors
                 if (is_trunc(i)) cycle
                 do j=1, npart
                     call dealloc_one_sqrmat(saved_a(j,i)) 
                 enddo
             enddo
             deallocate(saved_a)
         endif

         if ( allocated(saved_b) ) then
             do i=1, nsectors
                 if (is_trunc(i)) cycle
                 do j=1, npart
                     call dealloc_one_sqrmat(saved_b(j,i)) 
                 enddo
             enddo
             deallocate(saved_b)
         endif
        
         return
     end subroutine ctqmc_deallocate_memory_part

!>>> subroutine used to determine is_save 
     subroutine ctqmc_make_nparts(cmode, csize, index_t_loc, tau_s, tau_e)
        use constants 
        use control
        use context
    
        use m_sector
   
        implicit none
   
! external arguments
! the mode of how to calculating trace
        integer,  intent(in)  :: cmode
   
! the total number of operators for current diagram
        integer,  intent(in)  :: csize
   
! local version of index_t
        integer, intent(in) :: index_t_loc(mkink)
   
! imaginary time value of operator A, only valid in cmode = 1 or 2
        real(dp), intent(in) :: tau_s
   
! imaginary time value of operator B, only valid in cmode = 1 or 2
        real(dp), intent(in) :: tau_e
        
! local variables
! length in imaginary time axis for each part
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
        first_fpart = 0
 
! is_save: how to process each part for each success string
! is_save = 0: the matrices product for this part has been calculated
!              previoulsy, another part_indx is needed to label which
!              part we can fetch and use it directly.
! is_save = 1: this part should be recalculated, and the result must be
!              stored in saved_a, and saved_a_nm arrays if this Monte Caro
!              move has been accepted.
! is_save = 2: this part is empty, we don't need to do anything with them.
! first, set it to be 0
       is_save(:,:,1) = is_save(:,:,2)
   
!--------------------------------------------------------------------
! when npart > 1, we use npart alogithm
! otherwise, recalculate all the matrices products
        if ( npart == 1 ) then
            nop(1) = csize
            ops(1) = 1
            ope(1) = csize
            first_fpart = 1
            if (nop(1) <= 0) then
                is_save(1,:,1) = 2
            else
                is_save(1,:,1) = 1
            endif
        elseif ( npart > 1) then
   
            interval = beta / real(npart)
! calculate number of operators for each part
            do i=1,csize
                j = ceiling( time_v( index_t_loc(i) ) / interval )
                nop(j) = nop(j) + 1
            enddo 
! if no operators in this part, ignore them
            do i=1, npart
                if (first_fpart == 0 .and. nop(i) > 0) then
                    first_fpart = i
                endif
                if (nop(i) <= 0) then
                    is_save(i,:,1) = 2 
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
   
! when cmode == 1 or comde == 2, we can use some saved matrices products 
! by previous accepted Monte Carlo move
         if (cmode == 1 .or. cmode == 2) then
! get the position of operator A and operator B
             tis = ceiling( tau_s / interval )
             tie = ceiling( tau_e / interval )
! operator A:
                if (nop(tis)>0) then
                    is_save(tis,:,1) = 1
                endif
! special attention: if operator A is on the left or right boundary, then
! the neighbour part should be recalculated as well
                if ( nop(tis) > 0 ) then
                    if ( tau_s >= time_v( index_t_loc( ope(tis) ) ) ) then
                        tip = tis + 1
                        do while ( tip <= npart )
                            if ( nop(tip) > 0 ) then
                                is_save(tip,:,1) = 1;  EXIT
                            endif
                            tip = tip + 1
                        enddo ! over do while loop
                    endif
! for remove an operator, nop(tis) may be zero
                else
                    tip = tis + 1
                    do while ( tip <= npart )
                        if ( nop(tip) > 0 ) then
                            is_save(tip,:,1) = 1; EXIT
                        endif
                        tip = tip + 1
                    enddo ! over do while loop
                endif ! back if ( nop(tis) > 0 ) block
   
! operator B:
                if (nop(tie)>0) then
                    is_save(tie,:,1) = 1
                endif
! special attention: if operator B is on the left or right boundary, then
! the neighbour part should be recalculated as well
                if ( nop(tie) > 0 ) then
                    if ( tau_e >= time_v( index_t_loc( ope(tie) ) ) ) then
                        tip = tie + 1
                        do while ( tip <= npart )
                            if ( nop(tip) > 0 ) then
                                is_save(tip,:,1) = 1; EXIT
                            endif
                            tip = tip + 1
                        enddo ! over do while loop
                    endif
! for remove an operator, nop(tie) may be zero
                else
                    tip = tie + 1
                    do while ( tip <= npart )
                        if ( nop(tip) > 0 ) then
                            is_save(tip,:,1) = 1; EXIT
                        endif
                        tip = tip + 1
                    enddo ! over do while loop
                endif ! back if ( nop(tie) > 0 ) block
   
! when cmode == 3 or cmode == 4, recalculate all the matrices products 
            elseif (cmode == 3 .or. cmode == 4) then
                do i=1, nsectors
                    do j=1, npart
                        if (is_save(j,i,1) == 0) then
                            is_save(j,i,1) = 1
                        endif  
                    enddo
                enddo
            endif ! back if (cmode == 1 .or. cmode == 2) block
   
! npart should be larger than zero
        else
            call ctqmc_print_error('ctqmc_make_ztrace', 'npart is small than 1, &
                                    it should be larger than zero')
        endif ! back if (npart == 1) block
!--------------------------------------------------------------------
   
        return
     end subroutine ctqmc_make_nparts
   
!>>> calculate the trace for one sector
     subroutine cat_sector_ztrace(csize, string, index_t_loc, expt_t_loc, trace)
        use constants
        use control
        use context
   
        use m_sector
   
        implicit none
   
! external variables
! the number of total fermion operators
        integer, intent(in) :: csize
   
! the string for this sector
        integer, intent(in) :: string(csize+1)
   
! the address index of fermion operators
        integer, intent(in) :: index_t_loc(mkink)
   
! the diagonal elements of last time-evolution matrices
        real(dp), intent(in) :: expt_t_loc(ncfgs)
   
! the calculated trace of this sector
        real(dp), intent(out) :: trace
   
! local variables
! temp matrices
        real(dp) :: right_mat(max_dim_sect_trunc, max_dim_sect_trunc)
        real(dp) :: tmp_mat(max_dim_sect_trunc, max_dim_sect_trunc)
   
! temp index
        integer :: dim1, dim2, dim3, dim4
        integer :: isect, sect1, sect2
        integer :: indx
        integer :: counter
        integer :: vt, vf
   
! loop index
        integer :: i,j,k,l
   
! init isect
        isect = string(1)
!--------------------------------------------------------------------
! from right to left: beta <------ 0
        dim1 = sectors(string(1))%ndim
        right_mat = zero 
        tmp_mat = zero
! loop over all the parts
        do i=1, npart
   
! this part has been calculated previously, just use its results
            if (is_save(i,isect,1) == 0) then
                sect1 = string(ope(i)+1)
                sect2 = string(ops(i))
                dim2 = sectors(sect1)%ndim
                dim3 = sectors(sect2)%ndim
                if (i > first_fpart) then
                    call dgemm( 'N', 'N', dim2, dim1, dim3,  one,    &
                                 saved_a(i,isect)%item, max_dim_sect_trunc, &
                                 right_mat,             max_dim_sect_trunc, &
                                 zero, tmp_mat,         max_dim_sect_trunc   )
                    right_mat(:, 1:dim1) = tmp_mat(:, 1:dim1)
                    num_prod = num_prod + one
                else
                    right_mat(:, 1:dim1) = saved_a(i, isect)%item(:,1:dim1)
                endif
  
! this part should be recalcuated 
            elseif (is_save(i,isect,1) == 1) then 
                sect1 = string(ope(i)+1)
                sect2 = string(ops(i))
                dim4 = sectors(sect2)%ndim
                saved_b(i,isect)%item = zero
   
! loop over all the fermion operators in this part
                counter = 0
                do j=ops(i), ope(i)
                    counter = counter + 1
                    indx = sectors(string(j  ))%istart
                    dim2 = sectors(string(j+1))%ndim
                    dim3 = sectors(string(j  ))%ndim
   
                    if (counter > 1) then
                        do l=1,dim4
                            do k=1,dim3
                                tmp_mat(k,l) = saved_b(i,isect)%item(k,l) * expt_v(indx+k-1, index_t_loc(j))
                            enddo
                        enddo
                        num_prod = num_prod + one
                    else
                        tmp_mat = zero
                        do k=1,dim3
                            tmp_mat(k,k) = expt_v(indx+k-1, index_t_loc(j))
                        enddo
                    endif
   
                    vt = type_v( index_t_loc(j) )
                    vf = flvr_v( index_t_loc(j) ) 
                    call dgemm( 'N', 'N', dim2, dim4, dim3, one,              &
                                sectors(string(j))%myfmat(vf, vt)%item, dim2, &
                                tmp_mat,                         max_dim_sect_trunc, &
                                zero, saved_b(i,isect)%item,     max_dim_sect_trunc   ) 
   
                    num_prod = num_prod + one
                enddo  

                is_save(i, isect, 1) = 0
                is_copy(i, isect) = .true.
                col_copy(i,isect) = dim4
   
! multiply this part with the rest parts
                if (i > first_fpart) then
                    call dgemm( 'N', 'N', dim2, dim1, dim4, one,    &
                                saved_b(i,isect)%item, max_dim_sect_trunc, &
                                right_mat,            max_dim_sect_trunc, &
                                zero, tmp_mat,        max_dim_sect_trunc   ) 
                    right_mat(:,1:dim1) = tmp_mat(:,1:dim1)
                    num_prod = num_prod + one
                else
                    right_mat(:,1:dim1) = saved_b(i, isect)%item(:,1:dim1)
                endif

            elseif (is_save(i,isect,1) == 2) then
                cycle
            endif ! back if ( is_save(i,isect) ==0 )  block
   
! the start sector for next part
            isect = sect1
   
        enddo  ! over i={1, npart} loop   
   
! special treatment of the last time-evolution operator
        indx = sectors(string(1))%istart
        if (csize == 0) then
            do k=1, dim1
                right_mat(k,k) = expt_t_loc(indx+k-1)
            enddo
        else
            do l=1,dim1
                do k=1,dim1
                    right_mat(k,l) = right_mat(k,l) * expt_t_loc(indx+k-1)
                enddo
            enddo
            num_prod = num_prod + one
        endif
   
! store final product
        final_product(string(1),1)%item = right_mat(1:dim1, 1:dim1)
   
! calculate the trace
        trace  = zero
        do j=1, sectors(string(1))%ndim
            trace = trace + right_mat(j,j)
        enddo
   
        return
     end subroutine cat_sector_ztrace

!>>> copy data when propose has been accepted
     subroutine ctqmc_save_parts()
        implicit none

! loop index
        integer :: i,j

! copy save-state for all the parts 
        is_save(:,:,2) = is_save(:,:,1)

! when npart > 1, we used the npart algorithm, save the changed 
! matrices products when moves are accepted
        if ( npart > 1) then
            do i=1, nsectors
                if (is_trunc(i)) cycle
                do j=1, npart
                    if ( is_copy(j,i) ) then
                        saved_a(j, i)%item(:,1:col_copy(j,i)) = saved_b(j, i)%item(:,1:col_copy(j,i)) 
                    endif
                enddo
            enddo
        endif

        return
     end subroutine ctqmc_save_parts

  end module m_npart
