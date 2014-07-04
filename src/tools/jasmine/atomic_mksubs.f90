subroutine atomic_make_subspaces()
    use constants
    use control
    use mod_global
    use mod_subspace
    implicit none

    ! local variables
    type(type_gqn) :: subs(maxsubs) 
    type(type_gqn) :: tmp
    integer :: ndims(maxsubs)
    integer :: sub_basis(maxdim, maxsubs)
    integer :: ntot
    integer :: myjz 
    integer :: ibasis
    integer :: isubs
    logical :: eql
    integer :: i,j,k 
     
    ! loop over all of the configurations to determine subspaces
    ibasis = 0
    nsubs = 0
    ndims = 0
    sub_basis = 0
    ! loop over subspaces by goog quantum number $N$
    do i=0, norbs
        ntot = i
        ! nstat(ntot): the dimension of ntot subspace 
        do j=1, nstat(ntot) 
            ibasis = ibasis + 1
            myjz = 0
            ! determine the jz
            do k=1, norbs
                myjz = myjz + good(k) * invcd(k, ibasis) 
            enddo  
            tmp = type_gqn(ntot, myjz)
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
                    call compare_gqn(subs(k), tmp, eql)   
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

    ! DEBUG 
    !-----------------------------------------------------------------
    open(mytmp, file="test_subs.dat") 
    do i=1, nsubs
        write(mytmp,"(4i5)") i, subs(i)%ntot, subs(i)%jz, ndims(i)
        do j=1, ndims(i)
            write(mytmp,"(i3,a,10i1)") j,": ", invcd(:, sub_basis(j, i)) 
        enddo
    enddo
    close(mytmp)
    !-----------------------------------------------------------------

    ! after we know how many subspaces and the dimension of each subspace,
    ! we can allocate memory for global variables subspaces
    call alloc_glob_subspaces()

    ! now we will build each subspaces 
    do i=1, nsubs
        subspaces(i)%indx = i
        subspaces(i)%ndim = ndims(i)
        subspaces(i)%mygqn = subs(i)
        ! allocate memory for pointers in each subspace
        ! WARNING: these memory should be deallocated before deallocating subspaces
        call alloc_one_subspace(subspaces(i))  
        ! set basis for each subspace
        do j=1, ndims(i)
            subspaces(i)%mybasis(j) = sub_basis(j,i) 
        enddo
    enddo

    return
end subroutine atomic_make_subspaces

subroutine atomic_make_towhich()
    use mod_subspace
    use mod_global
    implicit none

    integer :: i, j, k, m
    integer :: isub, ibasis
    integer :: ntot_tmp, jz_tmp
    type(type_gqn) :: gqn_tmp
    logical :: eql
    logical :: eql2
   
    ! for create operators
    do i=1, norbs
        do j=1, nsubs
            isub = -1
            ! we should lookup each basis in this subspace 
            eql2 = .false.
            do m=1, subspaces(j)%ndim
                ibasis = subspaces(j)%mybasis(m)
                if (invcd(i,ibasis)==0) then
                    eql2 = .true.
                    exit
                endif 
            enddo 
            if (eql2 == .true.) then
                ntot_tmp = subspaces(j)%mygqn%ntot + 1
                jz_tmp   = subspaces(j)%mygqn%jz + good(i)
                gqn_tmp = type_gqn(ntot_tmp, jz_tmp)
                ! loop over all subspaces to see which subspace it will be 
                do k=1, nsubs
                    call compare_gqn(gqn_tmp, subspaces(k)%mygqn, eql) 
                    if (eql) then
                        isub = k
                        exit 
                    endif 
                enddo 
            endif

            c_towhich(j,i) = isub 
        enddo
    enddo
    ! for destroy operators
    do i=1, norbs
        do j=1, nsubs
            isub = -1
            ! we should lookup each basis in this subspace 
            eql2 = .false.
            do m=1, subspaces(j)%ndim
                ibasis = subspaces(j)%mybasis(m)
                if (invcd(i,ibasis)==1) then
                    eql2 = .true.
                    exit
                endif 
            enddo 
            if (eql2 == .true.) then
                ntot_tmp = subspaces(j)%mygqn%ntot - 1
                jz_tmp   = subspaces(j)%mygqn%jz - good(i)
                gqn_tmp = type_gqn(ntot_tmp, jz_tmp)
                ! loop over all subspaces to see which subspace it will be 
                do k=1, nsubs
                    call compare_gqn(gqn_tmp, subspaces(k)%mygqn, eql) 
                    if (eql) then
                        isub = k
                        exit 
                    endif 
                enddo 
            endif

            d_towhich(j,i) = isub 
        enddo
    enddo

    ! DEBUG
    !-------------------------------------------------------------------
    open(mytmp, file="test_ctowhich.dat")
    do i=1, norbs
        do j=1, nsubs
            if (c_towhich(j,i) /= -1) then
            write(mytmp,"(8i3)") i, good(i), j, subspaces(j)%mygqn%ntot,subspaces(j)%mygqn%jz, &
                                c_towhich(j,i), &  
                            subspaces(c_towhich(j,i))%mygqn%ntot,subspaces(c_towhich(j,i))%mygqn%jz
            else
            write(mytmp,"(5i3)") i, good(i), j, subspaces(j)%mygqn%ntot,subspaces(j)%mygqn%jz
 
            endif
        enddo
    enddo
    close(mytmp)
    open(mytmp, file="test_dtowhich.dat")
    do i=1, norbs
        do j=1, nsubs
            if (d_towhich(j,i) /= -1) then
            write(mytmp,"(8i3)") i, good(i), j, subspaces(j)%mygqn%ntot,subspaces(j)%mygqn%jz, &
                                d_towhich(j,i), &  
                            subspaces(d_towhich(j,i))%mygqn%ntot,subspaces(d_towhich(j,i))%mygqn%jz
            else
            write(mytmp,"(5i3)") i, good(i), j, subspaces(j)%mygqn%ntot,subspaces(j)%mygqn%jz
 
            endif

        enddo
    enddo
    close(mytmp)
    !-------------------------------------------------------------------
    return
end subroutine atomic_make_towhich

subroutine atomic_trunk_space()
    use constants
    use control
    use mod_global
    implicit none

    integer :: i, j
    integer :: ntot1
    integer :: ntot2
    integer :: ntot3

    c_towhich_trunk = c_towhich
    d_towhich_trunk = d_towhich
    do i=1, nsubs
        ntot1 = subspaces(i)%mygqn%ntot 
        ntot2 = ntot1 + 1 
        ntot3 = ntot1 - 1 
        ! for create operator
        if ( .not. (ntot1 >= nmin .and. ntot2 <= nmax) ) then
            c_towhich_trunk(i,:) = -1
        endif
        ! for destroy operator
        if ( .not. (ntot3 >= nmin .and. ntot1 <= nmax) ) then
            d_towhich_trunk(i,:) = -1
        endif
    enddo
    ! DEBUG
    open(mytmp, file="test_ctowhich_trunk.dat")
    do i=1, norbs
        do j=1, nsubs
            if (c_towhich_trunk(j,i) /= -1) then
            write(mytmp,"(8i5)") i, good(i), j, subspaces(j)%mygqn%ntot,subspaces(j)%mygqn%jz, &
                                c_towhich_trunk(j,i), &  
                            subspaces(c_towhich_trunk(j,i))%mygqn%ntot,subspaces(c_towhich_trunk(j,i))%mygqn%jz
            else
            write(mytmp,"(5i5)") i, good(i), j, subspaces(j)%mygqn%ntot,subspaces(j)%mygqn%jz
 
            endif
        enddo
    enddo
    close(mytmp)
    ! DEBUG
end subroutine atomic_trunk_space
