subroutine atomic_write_eigval()
    use constants
    use control
    use mod_global
    implicit none

    integer :: i, j
    integer :: counter

    open(mytmp, file="atom.eigval.dat")
    counter = 0
    do i=1, nsubs
        do j=1, subspaces(i)%ndim
            counter = counter + 1
            write(mytmp, "(2i5, f16.8)") counter, j, subspaces(i)%myeigval(j)
        enddo
    enddo
    close(mytmp)

    return
end subroutine atomic_write_eigval

subroutine atomic_write_eigvec()
    use constants
    use control
    use mod_global
    implicit none

    integer :: i, j, k
    integer :: counter

    open(mytmp, file="atom.eigvec.dat")
    counter = 0
    do i=1, nsubs
        do j=1, subspaces(i)%ndim
            do k=1, subspaces(i)%ndim
                write(mytmp, "(3i5,f16.8)") i, k+counter, j+counter, subspaces(i)%myeigvec(k,j)
            enddo
        enddo
        counter = counter + subspaces(i)%ndim
    enddo
    close(mytmp)

    return
end subroutine atomic_write_eigvec

subroutine atomic_write_fmat()
    use constants
    use control
    use mod_global
    implicit none

    integer :: iorb, isub
    integer :: i, j

    open(mytmp, file="atom.fmat.dat")
    do iorb=1, norbs
        do isub=1, nsubs
            if (c_towhich(isub, iorb) == -1) cycle
            write(mytmp, "(5i5)") iorb, isub, c_towhich(isub, iorb), &
                    subspaces(isub)%ndim, subspaces(c_towhich(isub,iorb))%ndim 
            do i=1, c_fmat(isub, iorb)%ndimx 
                do j=1, c_fmat(isub, iorb)%ndimy
                    write(mytmp, "(2i5,f16.8)") j, i, c_fmat(isub, iorb)%elem(j,i)
                enddo
            enddo
        enddo
    enddo
    close(mytmp)

    return
end subroutine atomic_write_fmat

subroutine atomic_write_ctqmc()
    use constants
    use control
    use mod_global
    implicit none

    integer :: isub, iorb
    integer :: i,j 

    open(mytmp, file="atom.cix")
    ! first, dump the eigval of each subspaces
    write(mytmp,"(a)") "# eigenvalues of each subspaces"
    do isub=1, nsubs
        write(mytmp,"(2i5)") isub, subspaces(isub)%ndim 
        do i=1, subspaces(isub)%ndim
            write(mytmp, "(i5, f20.14)") i, subspaces(isub)%myeigval(i)
        enddo
        ! one blank line
        write(mytmp, *)
    enddo
    ! then, dump the c_towhich, c_towhich_trunk, d_towhich, d_towhich_trunk
    write(mytmp,"(a)") "# index from one subspace to another subspace, and fmat"
    do iorb=1, norbs
        do isub=1, nsubs
            ! index from one subspace to another subspace
            write(mytmp, "(6i5)") iorb, isub, c_towhich(isub, iorb), c_towhich_trunk(isub, iorb), &
                                              d_towhich(isub, iorb), d_towhich_trunk(isub, iorb)
            ! fmat of creation operator 
            if (c_towhich(isub, iorb) /= -1) then 
                write(mytmp, "(2i)") subspaces(c_towhich(isub,iorb))%ndim, subspaces(isub)%ndim
                do i=1, c_fmat(isub, iorb)%ndimx
                    do j=1, c_fmat(isub, iorb)%ndimy
                        write(mytmp, "(2i5, f20.14)") j, i, c_fmat(isub,iorb)%elem(j,i)
                    enddo
                enddo
            endif

            ! one blank line
            write(mytmp, *) 
        enddo
    enddo
    close(mytmp)
end subroutine atomic_write_ctqmc
