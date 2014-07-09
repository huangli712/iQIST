!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_mksectors_n
!         : atomic_mksectors_nszps
!         : atomic_mksectors_njz
! source  : atomic_mksector.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : make sectors by using good quantum numbers
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> determine all the sectors for good quantum N case
! a sector consists of some many particle Fock states labeled by 
! good quantum number N 
subroutine atomic_mksectors_n()
    implicit none

    print *, "not implemented!"

    return
end subroutine atomic_mksectors_n

!>>> determine all the sectors for good quantum N, Sz, PS case
! a sector consists of some many particle Fock states labeled by 
! good quantum number N, Sz, PS
subroutine atomic_mksectors_nszps()
    implicit none

    print *, "not implemented!"

    return
end subroutine atomic_mksectors_nszps

!>>> determine all the sectors for good quantum N, Jz case
! a sector consists of some many particle Fock states labeled by 
! good quantum number N, Jz 
subroutine atomic_mksectors_njz()
    use control, only: norbs, ncfgs
    use m_basis_fullspace, only: dim_sub_n, bin_basis
    use m_sector
    use m_glob_sectors

    implicit none

    ! local variables
    ! the maximum number of sectors
    integer :: max_nsect
    ! the maximum dimension of each sector
    integer :: max_ndim
    ! the jz value for each |j2,jz> single particle orbital
    integer :: orb_good_jz(norbs)
    ! good quantum number N, Jz for each Fock state
    integer :: fock_good_ntot(ncfgs)
    integer :: fock_good_jz(ncfgs)
    ! good quantum number N, Jz for each sector
    integer, allocatable :: sect_good_ntot(:)
    integer, allocatable :: sect_good_jz(:)
    ! dimension of each sector
    integer, allocatable :: ndims(:)
    ! sector basis index
    integer, allocatable :: sector_basis(:,:)
    ! number of sectors
    integer :: nsect
    ! which sector point to
    integer :: which_sect
    ! tmp variables
    integer :: myntot
    integer :: myjz
    integer :: counter
    integer :: ibasis
    integer :: i,j,k,l
    logical :: can  


    max_nsect = ncfgs
    max_ndim = ncfgs
    ! allocate memory
    allocate(sect_good_ntot(max_nsect))
    allocate(sect_good_jz(max_nsect))
    allocate(ndims(max_nsect))
    allocate(sector_basis(max_ndim, max_nsect))

    ! make good_jz
    call make_good_jz(orb_good_jz)

    ! build good quantum numbers for each Fock state
    counter = 0
    do i=0, norbs
        do j=1, dim_sub_n(i)
            counter = counter + 1
            myjz = 0
            do k=1, norbs
                myjz = myjz + orb_good_jz(k) * bin_basis(k, counter) 
            enddo
            fock_good_ntot(counter) = i
            fock_good_jz(counter) = myjz
        enddo  
    enddo

    !----------------------------------------------------------------
    ! loop over all the Fock states to determine sectors
    nsect = 0
    ndims = 0
    sector_basis = 0
    do i=1, ncfgs    
        myntot = fock_good_ntot(i)
        myjz   = fock_good_jz(i)
        if (nsect==0) then
            sect_good_ntot(1) = myntot
            sect_good_jz(1)   = myjz
            nsect = nsect + 1
            ndims(1) = ndims(1) + 1 
            sector_basis(ndims(1),1) = i
        else
            ! loop over the exists sectors
            which_sect = -1
            do j=1, nsect
                ! compare two subspaces
                if ( sect_good_ntot(j) == myntot .and. sect_good_jz(j) == myjz) then
                    which_sect = j
                    EXIT
                endif
            enddo 
            ! new sector
            if( which_sect == -1 ) then
                nsect = nsect + 1
                sect_good_ntot(nsect) = myntot
                sect_good_jz(nsect)   = myjz
                ndims(nsect) = ndims(nsect) + 1
                sector_basis(ndims(nsect), nsect) = i
            ! old sector
            else
                ndims(which_sect) = ndims(which_sect) + 1 
                sector_basis(ndims(which_sect), which_sect) = i
            endif
        endif ! back to if (nsect == 0) then block 
    enddo 

    !----------------------------------------------------------------
    ! after we know how many sectors and the dimension of each sector,
    ! we can allocate memory for global variables for sectors
    nsectors = nsect
    call alloc_m_glob_sectors()
    ! now we will build each sector
    counter = 1
    do i=1, nsect
        sectors(i)%ndim = ndims(i)
        sectors(i)%nelectron = sect_good_ntot(i)
        sectors(i)%nops = norbs
        sectors(i)%istart = counter 
        counter = counter + ndims(i)
        ! allocate memory for pointers in each subspace
        ! WARNING: these memory should be deallocated before deallocating sectors
        call alloc_one_sector( sectors(i) )  
        ! set basis for each sector
        do j=1, ndims(i)
            sectors(i)%mybasis(j) = sector_basis(j,i) 
        enddo
    enddo

    !----------------------------------------------------------------
    ! make next_sector index
    ! for create operators
    do i=1, nsectors
        do j=1, norbs 
            do k=0,1 
                which_sect = -1
            ! we should lookup each basis in this subspace 
                can = .false.
                do l=1, sectors(i)%ndim
                    ibasis = sectors(j)%mybasis(l)
                    if (k==1 .and. bin_basis(j,ibasis) == 0 ) then
                        can = .true.
                        exit
                    elseif (k==0 .and. bin_basis(j, ibasis) == 1) then
                        can = .true. 
                        exit
                    endif 
                enddo 
                if (can == .true.) then
                    if (k==1) then
                        myntot = sect_good_ntot(i) + 1
                        myjz   = sect_good_jz(i) + orb_good_jz(j)
                    else
                        myntot = sect_good_ntot(i) - 1
                        myjz   = sect_good_jz(i) - orb_good_jz(j)
                    endif
                    ! loop over all sectors to see which sector it will point to 
                    do l=1, nsectors
                        if (sect_good_ntot(l) == myntot .and. sect_good_jz(l) == myjz) then
                            which_sect = l
                            exit 
                        endif 
                    enddo 
                endif  ! back to if (can == .true.) block
                ! set next_sector index
                sectors(i)%next_sector(j,k) = which_sect 
            enddo ! over k={0,1} loop
        enddo ! over j={1,norbs} loop
    enddo ! over i={1, nsectors} loop

    ! free memeory
    if (allocated(sect_good_ntot)) deallocate(sect_good_ntot)
    if (allocated(sect_good_jz))   deallocate(sect_good_jz)
    if (allocated(ndims))          deallocate(ndims)
    if (allocated(sector_basis))   deallocate(sector_basis)

    return
end subroutine atomic_mksectors_njz

subroutine make_good_jz(good_jz)
    use control, only: norbs

    implicit none

    ! external variables
    integer, intent(out) :: good_jz(norbs)

    if (norbs == 6) then
        ! j=1/2
        good_jz(1) = -1
        good_jz(2) =  1
        ! j=3/2
        good_jz(3) = -3
        good_jz(4) = -1
        good_jz(5) =  1
        good_jz(6) =  3
    elseif (norbs == 10) then
        ! j=3/2
        good_jz(1) = -3
        good_jz(2) = -1
        good_jz(3) =  1
        good_jz(4) =  3
        ! j=5/2
        good_jz(5) = -5
        good_jz(6) = -3
        good_jz(7) = -1
        good_jz(8) =  1
        good_jz(9) =  3
        good_jz(10)=  5
    elseif (norbs == 14) then
        ! j=5/2
        good_jz(1) = -5
        good_jz(2) = -3
        good_jz(3) = -1
        good_jz(4) =  1
        good_jz(5) =  3
        good_jz(6) =  5
        ! j=7/2
        good_jz(7) = -7
        good_jz(8) = -5
        good_jz(9) = -3
        good_jz(10)= -1
        good_jz(11)=  1
        good_jz(12)=  3
        good_jz(13)=  5
        good_jz(14)=  7
    else
        call atomic_print_error('make_good_jz', 'not implemented for this norbs value !')
    endif

    return
end subroutine make_good_jz


