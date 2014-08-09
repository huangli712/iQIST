!-------------------------------------------------------------------------
! project : manjushaka
! program : m_sector
! source  : mod_control.f90
! type    : modules
! authors : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014
!           07/19/2014
! purpose : define data structure for good quantum number algorithm
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> data structure for good quantum number algorithm
  module m_sector
     use constants
     use control
     use context
 
     implicit none
  
! the fmat between any two sectors, it is just a matrix
     type :: t_fmat

! the dimension
         integer :: n, m

! the items of the matrix
         real(dp), dimension(:,:), pointer :: item => null()

     end type t_fmat
  
! a square matrix
     type :: t_sqrmat

! the dimension
         integer :: n

! the items of the matrix
         real(dp), dimension(:,:), pointer :: item => null()

     end type t_sqrmat

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
         real(dp), dimension(:), pointer :: myeigval => null()

! the next sector it points to when a fermion operator acts on this sector
! -1: outside of the Hilbert space, otherwise, it is the index of next sector
! next_sector(nops,0:1), 0 for annihilation and 1 for creation operators, respectively
! F|i> --> |j>
         integer, dimension(:,:), pointer :: next_sector => null()

! this is for truncating the Hilbert space
         integer, dimension(:,:), pointer :: next_sector_trunc => null()

! the fmat between this sector and all other sectors
! if this sector doesn't point to some other sectors, the pointer is null
! mymfat(nops, 0:1), 0 for annihilation and 1 for creation operators, respectively
         type(t_fmat), dimension(:,:), pointer :: myfmat => null()

     end type t_sector
     
! the total number of sectors
     integer, public, save :: nsectors

! total number of sectors after truncated
     integer, public, save :: nsectors_trunc 

! the max dimension of the sectors
     integer, public, save :: max_dim_sect

! max dimension of the sectors after truncated
     integer, public, save :: max_dim_sect_trunc

! the average dimension of the sectors
     real(dp), public, save :: ave_dim_sect

! average dimension of the sectors after truncated
     real(dp), public, save :: ave_dim_sect_trunc

! the array contains all the sectors
     type(t_sector), public, save, allocatable :: sectors(:)

! the probability of each sector
     real(dp), public, save, allocatable :: prob_sect(:)

! which sectors should be truncated ?
     logical, public, save, allocatable :: is_trunc(:)

! the final product matrices, which will be used to calculate the nmat
     type(t_sqrmat), public, save, allocatable :: final_product(:,:)

! matrices of occupancy operator c^{\dagger}c 
     type(t_sqrmat), public, save, allocatable :: occu(:,:)

! matrices of double occupancy operator c^{\dagger}cc^{\dagger}c
     type(t_sqrmat), public, save, allocatable :: double_occu(:,:,:)

     contains
  
!>>> allocate one fmat
     subroutine alloc_one_fmat(one_fmat)
        implicit none
  
! external variables
        type(t_fmat), intent(inout) :: one_fmat
  
        allocate(one_fmat%item(one_fmat%n, one_fmat%m))
  
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

!>>> allocate one fmat
     subroutine alloc_one_sqrmat(one_sqrmat)
        implicit none
  
! external variables
        type(t_sqrmat), intent(inout) :: one_sqrmat
  
        allocate(one_sqrmat%item(one_sqrmat%n, one_sqrmat%n))
  
! initialize it
        one_sqrmat%item = zero
  
        return
     end subroutine alloc_one_sqrmat
  
!>>> deallocate one fmat
     subroutine dealloc_one_sqrmat(one_sqrmat)
        implicit none
  
! external variables
        type(t_sqrmat), intent(inout) :: one_sqrmat
  
        if ( associated(one_sqrmat%item) ) deallocate(one_sqrmat%item)
  
        return
     end subroutine dealloc_one_sqrmat
 
!>>> allocate memory for one sector
     subroutine alloc_one_sector(one_sector)
        implicit none
  
! external variables
        type(t_sector), intent(inout) :: one_sector
  
! local variables
        integer :: i, j
  
        allocate(one_sector%myeigval(one_sector%ndim))
        allocate(one_sector%next_sector(one_sector%nops,0:1))
        allocate(one_sector%next_sector_trunc(one_sector%nops,0:1))
        allocate(one_sector%myfmat(one_sector%nops,0:1))
  
! init them
        one_sector%myeigval = zero
        one_sector%next_sector = 0
        one_sector%next_sector_trunc = 0
  
! init myfmat one by one
        do i=1, one_sector%nops 
           do j=0, 1
               one_sector%myfmat(i,j)%n = 0
               one_sector%myfmat(i,j)%m = 0
               one_sector%myfmat(i,j)%item => null()
           enddo
        enddo

        return
     end subroutine alloc_one_sector
  
!>>> deallocate memory for onespace
     subroutine dealloc_one_sector(one_sector)
        implicit none
  
! external variables
        type(t_sector), intent(inout) :: one_sector 
  
! local variables  
        integer :: i, j
  
        if (associated(one_sector%myeigval))            deallocate(one_sector%myeigval)
        if (associated(one_sector%next_sector))         deallocate(one_sector%next_sector)
        if (associated(one_sector%next_sector_trunc))   deallocate(one_sector%next_sector_trunc)
  
! deallocate myfmat one by one
        if ( associated(one_sector%myfmat) ) then
            do i=1, one_sector%nops
                do j=0,1
                    call dealloc_one_fmat(one_sector%myfmat(i,j))
                enddo
            enddo 
            deallocate(one_sector%myfmat)
        endif
  
        return
     end subroutine dealloc_one_sector

!>>> allocate memory for sect-related variables
     subroutine ctqmc_allocate_memory_sect()
         implicit none

! local variables
         integer :: i
         integer :: istat

! allocate memory
         allocate(sectors(nsectors),    stat=istat)
         allocate(is_trunc(nsectors),   stat=istat)
         allocate(prob_sect(nsectors),  stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_sect','can not allocate enough memory')
         endif

! initialize them
         do i=1, nsectors
             sectors(i)%ndim = 0
             sectors(i)%nelectron = 0
             sectors(i)%nops = norbs
             sectors(i)%istart = 0
             sectors(i)%myeigval => null()
             sectors(i)%next_sector => null()
             sectors(i)%next_sector_trunc => null()
             sectors(i)%myfmat => null()
         enddo 
         is_trunc = .false.
         prob_sect = zero

         return
     end subroutine ctqmc_allocate_memory_sect

!>>> deallocate memory for sect-related variables
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

        if ( allocated(is_trunc) ) deallocate(is_trunc)
        if ( allocated(prob_sect) ) deallocate(prob_sect)

        return
     end subroutine ctqmc_deallocate_memory_sect

!>>> allocate memory for occu
     subroutine ctqmc_allocate_memory_occu()
        implicit none

        integer :: i,j,k

        allocate(final_product(nsectors,2))
        allocate(occu(norbs, nsectors))

        do i=1, nsectors
            if (is_trunc(i)) cycle
            
            do j=1,2
                final_product(i,j)%n = sectors(i)%ndim
                call alloc_one_sqrmat( final_product(i,j) )
            enddo

            do j=1, norbs
                occu(j,i)%n = sectors(i)%ndim
                call alloc_one_sqrmat( occu(j,i) )
            enddo

        enddo
 
        if (idoub == 2) then
            allocate(double_occu(norbs, norbs, nsectors))
            do i=1, nsectors
                if (is_trunc(i)) cycle
                do j=1, norbs
                    do k=1, norbs
                        double_occu(k,j,i)%n = sectors(i)%ndim
                        call alloc_one_sqrmat( double_occu(k,j,i) )
                    enddo
                enddo
            enddo
        endif
   
        return
     end subroutine ctqmc_allocate_memory_occu

!>>> deallocate memory for occu
     subroutine ctqmc_deallocate_memory_occu()
        implicit none

        integer :: i,j,k

        if ( allocated(final_product) ) then
            do i=1, nsectors
                if (is_trunc(i)) cycle
                do j=1,2
                    call dealloc_one_sqrmat(final_product(i,j))
                enddo
            enddo
            deallocate(final_product)
        endif

        if ( allocated(occu) ) then
            do i=1, nsectors
                if (is_trunc(i)) cycle
                do j=1, norbs
                    call dealloc_one_sqrmat( occu(j,i) )
                enddo
            enddo
            deallocate(occu)
        endif

        if ( allocated(double_occu) ) then
            do i=1, nsectors
                if (is_trunc(i)) cycle
                do j=1, norbs
                    do k=1, norbs
                        call dealloc_one_sqrmat( double_occu(k,j,i) ) 
                    enddo
                enddo
            enddo
            deallocate(double_occu)
        endif

        return
     end subroutine ctqmc_deallocate_memory_occu

!>>> subroutine used to truncate the Hilbert space
     subroutine ctqmc_make_trunc()
        implicit none

! local variables
! loop index
        integer :: i,int1

! file status
        logical :: exists

! don't truncate the Hilbert space at all
        if (itrun == 1) then
            do i=1, nsectors
                sectors(i)%next_sector_trunc = sectors(i)%next_sector
            enddo
            nsectors_trunc     = nsectors
            max_dim_sect_trunc = max_dim_sect
            ave_dim_sect_trunc = ave_dim_sect
            if (myid == master) then
                write(mystd,*)
                write(mystd,'(4X,a)') 'use full Hilbert space' 
                write(mystd,'(4X,a,i5)')    'number of sectors:', nsectors_trunc 
                write(mystd,'(4X,a,i5)')    'maximum dimension of sectors:', max_dim_sect_trunc
                write(mystd,'(4X,a,f10.2)') 'averaged dimension of sectors:', ave_dim_sect_trunc
                write(mystd,*)
            endif

! truncate the Hilbert space according to the total occupancy number
        elseif (itrun == 2) then
            is_trunc = .false.
            do i=1, nsectors
                if (sectors(i)%nelectron < nmini .or. sectors(i)%nelectron > nmaxi) then
                    is_trunc(i) = .true.
                endif
            enddo
            call truncate_sectors()
            if (myid == master) then
                write(mystd, *)
                write(mystd,'(4X,a,i2,a,i2)') 'truncate occupancy number, just keep ', &
                                                  nmini,' ~ ',  nmaxi 
                write(mystd,'(4X,a)') 'before truncated:'
                write(mystd,'(4X,a,i5)')    'number of sectors: ', nsectors
                write(mystd,'(4X,a,i5)')    'maximum dimension of sectors: ', max_dim_sect
                write(mystd,'(4X,a,f10.2)') 'averaged dimension of sectors: ', ave_dim_sect
                write(mystd,*)
                write(mystd,'(4X,a)') 'after truncated:'
                write(mystd,'(4X,a,i5)')    'number of sectors: ', nsectors_trunc 
                write(mystd,'(4X,a,i5)')    'maximum dimension of sectors: ', max_dim_sect_trunc
                write(mystd,'(4X,a,f10.2)') 'averaged dimension of sectors: ', ave_dim_sect_trunc
                write(mystd,*)
            endif

! truncate the Hilbert space according to the total occupancy number and energy level
        elseif (itrun == 3) then
            if (myid == master) then
                write(mystd,*)
                write(mystd,'(4X,a,i2,a,i2)') 'truncate occupancy number, just keep ', &
                                                  nmini, ' ~ ',  nmaxi

            endif
            is_trunc = .false.
            inquire(file = 'solver.psect.dat', exist = exists)
            if (exists) then
                if (myid == master) then
                    write(mystd,'(4X,a)') 'truncate high energy states'
                endif
                open(mytmp, file='solver.psect.dat', form='formatted', status='unknown')
                read(mytmp, *) ! skip header
                do i=1, nsectors
                    read(mytmp, *) int1, prob_sect(i) 
                enddo
                do i=1, nsectors
                    if (sectors(i)%nelectron < nmini .or. sectors(i)%nelectron > nmaxi .or. prob_sect(i) < 1e-4) then
                        is_trunc(i) = .true.
                    endif
                enddo
            else
                do i=1, nsectors
                    if (sectors(i)%nelectron < nmini .or. sectors(i)%nelectron > nmaxi ) then
                        is_trunc(i) = .true.
                    endif
                enddo
            endif
            call truncate_sectors()
            if (myid == master) then
                write(mystd,'(4X,a)') 'before truncated:'
                write(mystd,'(4X,a,i5)')    'number of sectors: ', nsectors
                write(mystd,'(4X,a,i5)')    'maximum dimension of sectors: ', max_dim_sect
                write(mystd,'(4X,a,f10.2)') 'averaged dimension of sectors: ', ave_dim_sect
                write(mystd,*)
                write(mystd,'(4X,a)') 'after truncated:'
                write(mystd,'(4X,a,i5)')    'number of sectors: ', nsectors_trunc 
                write(mystd,'(4X,a,i5)')    'maximum dimension of sectors: ', max_dim_sect_trunc
                write(mystd,'(4X,a,f10.2)') 'averaged dimension of sectors: ', ave_dim_sect_trunc
                write(mystd,*)
            endif


        endif

        return
     end subroutine ctqmc_make_trunc

     subroutine truncate_sectors()
        implicit none

        integer :: i,j,k,ii
        integer :: sum_dim
 
        max_dim_sect_trunc = -1
        nsectors_trunc = 0
        sum_dim = 0
        do i=1, nsectors
            sectors(i)%next_sector_trunc = -1
            if (is_trunc(i)) then
                cycle
            endif
            if (max_dim_sect_trunc < sectors(i)%ndim) then
                max_dim_sect_trunc = sectors(i)%ndim
            endif 
            sum_dim = sum_dim + sectors(i)%ndim
            nsectors_trunc = nsectors_trunc + 1 
            do j=1, sectors(i)%nops
                do k=0,1
                    ii = sectors(i)%next_sector(j,k) 
                    if (ii == -1) cycle
                    if (.not. is_trunc(ii)) then
                        sectors(i)%next_sector_trunc(j,k) = ii
                    endif
                enddo
            enddo
        enddo
        ave_dim_sect_trunc = real(sum_dim) / real(nsectors_trunc)

        return
     end subroutine truncate_sectors


!>>> calculate < c^{\dag} c>, < c^{\dag} c c^{\dag} c >
     subroutine ctqmc_make_occu()
        implicit none

        integer :: i,j,k,ii,jj
        real(dp) :: tmp_mat1(max_dim_sect_trunc, max_dim_sect_trunc) 
        real(dp) :: tmp_mat2(max_dim_sect_trunc, max_dim_sect_trunc) 

        do i=1, norbs
            do j=1, nsectors
                if ( is_trunc(j) ) cycle
                k=sectors(j)%next_sector(i,0)
                if (k == -1) then
                    occu(i,j)%item = zero
                    cycle
                endif
                call dgemm( 'N', 'N', sectors(j)%ndim, sectors(j)%ndim, sectors(k)%ndim, one, &
                            sectors(k)%myfmat(i,1)%item,                     sectors(j)%ndim, &
                            sectors(j)%myfmat(i,0)%item,                     sectors(k)%ndim, & 
                            zero, occu(i,j)%item,                   sectors(j)%ndim   ) 
            enddo
        enddo 

        if (idoub == 2) then
            do i=1, norbs
                do j=1, norbs
                    do k=1, nsectors
                        if ( is_trunc(k) ) cycle 
                        jj = sectors(k)%next_sector(j,0) 
                        ii = sectors(k)%next_sector(i,0)
                        if (ii == -1 .or. jj == -1) then
                            double_occu(i,j,k)%item = zero
                            cycle
                        endif
                        call dgemm( 'N', 'N', sectors(k)%ndim, sectors(k)%ndim, sectors(jj)%ndim, one, &
                                    sectors(jj)%myfmat(j,1)%item,                    sectors(k)%ndim,  & 
                                    sectors(k)%myfmat(j,0)%item,                     sectors(jj)%ndim, & 
                                    zero, tmp_mat1,                                  max_dim_sect_trunc ) 

                        call dgemm( 'N', 'N', sectors(k)%ndim, sectors(k)%ndim, sectors(ii)%ndim, one, &
                                    sectors(ii)%myfmat(i,1)%item,                    sectors(k)%ndim,  &
                                    sectors(k)%myfmat(i,0)%item,                     sectors(ii)%ndim, & 
                                    zero, tmp_mat2,                                  max_dim_sect_trunc ) 

                        call dgemm( 'N', 'N', sectors(k)%ndim, sectors(k)%ndim, sectors(k)%ndim, one,   &
                                    tmp_mat2,                                        max_dim_sect_trunc,& 
                                    tmp_mat1,                                        max_dim_sect_trunc,& 
                                    zero, double_occu(i,j,k)%item,          sectors(k)%ndim     )

                    enddo
                enddo
            enddo
        endif

        return
     end subroutine ctqmc_make_occu

!>>> subroutine used to build a string
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

!--------------------------------------------------------------------
        is_string = .true.
        string = -1

! we build a string from right to left, that is,  beta <------- 0
! build the string from the beginning sector, that is:
! S_a1(q1)-->q2, S_a2(q2)-->q3, ... S_ai(qi)-->qi+1, ..., Sak(qk)-->q1
! if we find some qi==0, we cycle this sector immediately
        do i=1,nsectors
            if (is_trunc(i)) then
                is_string(i) = .false.
            endif
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
                    next_sect_left = sectors(curr_sect_left)%next_sector_trunc(vf,vt)
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
                    next_sect_right = sectors(curr_sect_right)%next_sector_trunc(vf,vt)
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
