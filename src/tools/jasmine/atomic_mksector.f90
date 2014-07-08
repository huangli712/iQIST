!>>> determine all the sectors for good quantum N, Jz case
subroutine atomic_mksectors_njz(good_jz)
    use constants
    use control
    use m_basis_fullspace
    use m_good_quantum
    use m_sector
    use m_sectors_njz

    implicit none

    ! external variables
    integer, intent(in) :: good_jz(norbs)

    ! local variables
    type(t_njz) :: subs(ncfgs) 
    type(t_njz) :: tmp
    integer :: ndims(ncfgs)
    integer :: sub_basis(ncfgs, ncfgs)
    integer :: nsubs
    integer :: ntot
    integer :: myjz 
    integer :: ibasis
    integer :: isubs
    logical :: eql
    integer :: i,j,k 
    integer :: counter
     
    ! loop over all of the configurations to determine subspaces
    ibasis = 0
    nsubs = 0
    ndims = 0
    sub_basis = 0

    ! loop over subspaces by good quantum number $N$
    do i=0, norbs
        ntot = i
        ! nstat(ntot): the dimension of ntot subspace 
        do j=1, dim_sub_n(ntot) 
            ibasis = ibasis + 1
            myjz = 0
            ! determine the jz
            do k=1, norbs
                myjz = myjz + good_jz(k) * bin_basis(k, ibasis) 
            enddo  
            tmp = t_good_njz(ntot, myjz)
            if (nsubs==0) then
                subs(1) = tmp
                nsubs = nsubs + 1
                ndims(1) = ndims(1) + 1 
                sub_basis(ndims(1),1) = ibasis
            else
                ! loop over the exists subspace
                isubs = -1
                SUBLOOP: do k=1, nsubs
                ! compare two subspaces
                    eql = .false.
                    call compare_good_njz(subs(k), tmp, eql)   
                    if (eql) then
                        isubs = k
                        EXIT SUBLOOP 
                    endif
                enddo SUBLOOP 
                ! new subsapce
                if( isubs == -1 ) then
                    nsubs = nsubs + 1
                    subs(nsubs) = tmp 
                    ndims(nsubs) = ndims(nsubs) + 1
                    sub_basis(ndims(nsubs), nsubs) = ibasis
                ! old subspace
                else
                    ndims(isubs) = ndims(isubs) + 1 
                    sub_basis(ndims(isubs), isubs) = ibasis
                endif
            endif ! back to {if (nsubs==0) then}
        enddo ! back to {do j=1, nstat(ntot)} 
    enddo ! back to {do i=0, norbs}

    ! after we know how many subspaces and the dimension of each subspace,
    ! we can allocate memory for global variables subspaces
    nsectors = nsubs
    call alloc_m_sectors_njz()

    ! now we will build each subspace 
    counter = 1
    do i=1, nsubs
        sectors(i)%ndim = ndims(i)
        sectors(i)%ndim = subs(i)%ntot
        sectors(i)%nops = norbs
        sectors(i)%istart = counter 
        counter = counter + ndims(i)
        ! allocate memory for pointers in each subspace
        ! WARNING: these memory should be deallocated before deallocating subspaces
        call alloc_one_sector( sectors(i) )  
        ! set basis for each subspace
        do j=1, ndims(i)
            sectors(i)%mybasis(j) = sub_basis(j,i) 
        enddo
    enddo

    return
end subroutine atomic_mksectors_njz

subroutine atomic_mknextsector_njz(good_jz)
    use m_basis_fullspace
    use m_good_quantum
    use m_sector
    use m_sectors_njz

    implicit none

    ! external variables
    integer, intent(in) :: good_jz(norbs)

    integer :: i, j, k, m
    integer :: isub, ibasis
    integer :: ntot_tmp, jz_tmp
    type(t_good_njz) :: good_tmp
    logical :: eql
    logical :: eql2
   
    ! for create operators
    do i=1, norbs
        do j=1, nsectors
            isub = -1
            ! we should lookup each basis in this subspace 
            eql2 = .false.
            do m=1, sectors(j)%ndim
                ibasis = sectors(j)%mybasis(m)
                if ( bin_basis(i,ibasis) == 0 ) then
                    eql2 = .true.
                    exit
                endif 
            enddo 
            if (eql2 == .true.) then
                ntot_tmp = good_njz(j)%ntot + 1
                jz_tmp   = good_njz(j)%jz + good_jz(i)
                good_njz_tmp = t_good_njz(ntot_tmp, jz_tmp)

                ! loop over all subspaces to see which subspace it will be 
                do k=1, nsectors
                    call compare_good_njz(good_njz_tmp, good_njz(k), eql) 
                    if (eql) then
                        isub = k
                        exit 
                    endif 
                enddo 
            endif
            sectors(j)%next_sector(i,1) = isub 
        enddo
    enddo

    ! for destroy operators
    do i=1, norbs
        do j=1, nsectors
            isub = -1
            ! we should lookup each basis in this subspace 
            eql2 = .false.
            do m=1, sectors(j)%ndim
                ibasis = sectors(j)%mybasis(m)
                if (bin_basis(i,ibasis)==1) then
                    eql2 = .true.
                    exit
                endif 
            enddo 
            if (eql2 == .true.) then
                ntot_tmp = good_njz(j)%ntot - 1
                jz_tmp   = good_njz(j)%jz - good_jz(i)
                good_njz_tmp = t_good_njz(ntot_tmp, jz_tmp)
                ! loop over all subspaces to see which subspace it will be 
                do k=1, nsectors
                    call compare_gqn(good_njz_tmp, good_njz(k), eql) 
                    if (eql) then
                        isub = k
                        exit 
                    endif 
                enddo 
            endif
            sectors(j)%next_sector(i,0) = isub 
        enddo
    enddo

    return
end subroutine atomic_mknextsector_njz
