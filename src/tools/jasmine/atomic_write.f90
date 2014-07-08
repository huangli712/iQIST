!>>> write basis of fullspace to file 'atom.basis.dat'
subroutine atomic_write_basis()
    use constants
    use control
    use m_basis_fullspace

    implicit none

    ! local variables
    integer :: i

    ! open file 'atom.basis.dat' to write
    open(mytmp, file='atom.basis.dat')
    ! write the header
    do i=1, ncfgs
        write(mytmp, "(3I5,3X,14I2)") i, dec_basis(i), index_basis(i), bin_basis(:,i)   
    enddo 

    close(mytmp)

    return
end subroutine atomic_write_basis

!>>> write eigenvalue of fullspace to file 'atom.eigval.dat'
subroutine atomic_write_eigval_fullspace()
    use constants
    use control
    use m_glob_fullspace

    implicit none

    ! local variables
    integer :: i
    
    ! open file 'atom.eigval.dat' to write
    open(mytmp, file='atom.eigval.dat')
    do i=1, ncfgs
        write(mytmp, "(I5, F20.14)") i, hmat_eigval(i)
    enddo

    close(mytmp)

    return
end subroutine atomic_write_eigval_fullspace

!>>> write eigenvector of fullspace to file 'atom.eigvec.dat'
subroutine atomic_write_eigvec_fullspace()
    use constants
    use control
    use m_basis_fullspace
    use m_glob_fullspace

    implicit none

    ! local variables
    integer :: i, j
 
    ! open file 'atom.eigvec.dat' to write
    open(mytmp, file='atom.eigvec.dat')
    do i=1, ncfgs
        do j=1, ncfgs
            if ( abs(hmat_eigvec(j,i)) < eps6 ) cycle
            write(mytmp, "(2I5, F20.14, 4X, 14I1") j, i, hmat_eigvec(j,i), bin_basis(:,j) 
        enddo
    enddo 

    close(mytmp)

    return
end subroutine atomic_write_eigvec_fullspace

!>>> write atom.cix for CTQMC input
subroutine atomic_write_atomcix_fullspace()
    use constants
    use control
    use m_glob_fullspace

    implicit none

    ! local variables
    integer :: i, j, k
    integer :: s_order

    ! open file 'atom.cix' to write
    open(mytmp, file='atom.cix')

    ! write eigenvalues
    write(mytmp,'(a)') '# eigenvalues: index | energy | occupy | spin'
    do i=1,ncfgs
        write(mytmp,'(I5,3F20.14)') i, hmat_eigval(i), occu_mat(i,i), zero
    enddo 

    ! write fmat
    write(mytmp,'(a)') '# f matrix element: alpha | beta | orbital | fmat'
    do i=1,norbs
        ! for non-soc case, the spin order of CTQMC is like up, up, up, dn, dn, dn
        ! but the spin order of this program is up, dn, up, dn, up, dn
        ! we adjust it here
        ! for soc case, it doesn't matter
        if (isoc==0) then
            if (i <= nband) then
                s_order = 2*i-1  
            else
                s_order = 2*(i-nband)
            endif
        else 
            s_order = i
        endif 

        do j=1,ncfgs
            do k=1,ncfgs
                write(mytmp,'(3I5, F20.14)') k, j, i, anni_fmat(k, j, s_order)
            enddo 
        enddo 
    enddo 
    close(mytmp)

    return
end subroutine atomic_write_atomicx_fullspace

!>>> write eigenvalue of sectors to file 'atom.eigval.dat'
subroutine atomic_write_eigval_sectors(nsect, sectors)
    use constants
    use control
    use m_sector

    implicit none

    ! external variables
    integer, intent(in) :: nsect
    type(t_sector), intent(in) :: sectors(nsect)

    ! local variables
    integer :: i, j
    integer :: counter

    ! open file 'atom.eigval.dat' to write
    open(mytmp, file='atom.eigval.dat')
    counter = 0
    do i=1, nsect
        do j=1, sectors(i)%ndim
            counter = counter + 1
            write(mytmp, "(3I5, F20.14)") counter, i, j, sectors(i)%myeigval(j)
        enddo
    enddo
    close(mytmp)

    return
end subroutine atomic_write_eigval_sectors

!>>> write eigenvector of sectors to file 'atom.eigval.dat'
subroutine atomic_write_eigvec_sectors(nsect, sectors)
    use constants
    use control
    use m_sector

    implicit none

    ! external variables
    integer, intent(in) :: nsect
    type(t_sector), intent(in) :: sectors(nsect)

    ! local variables
    integer :: i, j
    integer :: counter

    open(mytmp, file="atom.eigvec.dat")
    counter = 0
    do i=1, nsect
        do j=1, sectors(i)%ndim
            do k=1, sectors(i)%ndim
                if ( abs(sectors(i)%myeigvec(k,j)) < eps6 ) cycle
                    write(mytmp, "(3I5, F20.14, 4X, 14I1)") i, k+counter, j+counter, &
                    sectors(i)%myeigvec(k,j), bin_basis(:,sectors(i)%mybasis(k))
            enddo
        enddo
        counter = counter + sectors(i)%ndim
    enddo
    close(mytmp)

    return
end subroutine atomic_write_eigvec_sectors

!>>> write atom.cix for CTQMC input, good quantum number algorithm
subroutine atomic_write_atomcix_sectors(nsect, sectors)
    use constants
    use control
    use m_sector

    implicit none

    ! external variables
    integer, intent(in) :: nsect
    type(t_sector), intent(in) :: sectors(nsect)

    ! local variables
    integer :: i, j, k, ii, row, col

    ! open 'atom.cix' to write
    open(mytmp, file='atom.cix')
    ! write number of sectors
    write(mytmp, "(a)") "#NUMBER OF SECTORS"
    write(mytmp, "(I5)") nsect

    ! write dimension, total electrons, next_sector, eigenvalue of each sector
    do i=1, nsect
        write(mytmp, "(a, I5)") "#SECT_INFO: ", i 
        write(mytmp, "(2I5, F20.14)") i, sectors(i)%ndim, sectors(i)%nelectrons

        ! write next_sector
        write(mytmp, "(2X,a)") "#NEXT_SECTOR"
        do j=1, sectors(i)%nops
            write(mytmp, "(2X, 3I5)") j, sectors(i)%next_sector(j,0), sectors(i)%next_sector(j,1)  
        enddo

        ! write eigeanvalue
        write(mytmp, "(2X,a)") "#EIGENVALUE"
        do j=1, sectors(i)%nops
            write(mytmp, "(2X,I5, F20.14)") j, sectors(i)%myeigval(j) 
        enddo
    enddo

    ! write fmat
    write(mytmp, "(a)") "#BEGIN_WRITE_FMAT"
    do i=1, nsect
        write(mytmp, "(2X,a, I5)") "#SECT_FMAT: ", i
        do j=1, sectors(i)%nops
            write(mytmp, "(4X,a, I5)") "#ORBIT: ", j
            do k=0,1
                write(mytmp, "(6X,a, I5)") "#SPIN: ", k
                ii = sectors(i)%myfmat(j,k)
                do col=1, sectors(i)%ndim
                    do row=1, sectors(ii)%ndim
                        write(mytmp, "(6X,2I5,F20.14)") row, col, sectors(i)%myfmat(j,k)%item(row, col)
                    enddo
                enddo
            enddo  ! over k={0,1} loop
        enddo ! over j={1, sectors(i)%nops} loop
    enddo  ! over i={1, nsect} loop

    close(mytmp)
 
    return
end subroutine atomic_write_atomcix_sectors
