!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_write_basis
!         : atomic_write_eigval_fullspace
!         : atomic_write_eigvec_fullspace
!         : atomic_write_atomcix_fullspace
!         : atomic_write_eigval_sectors
!         : atomic_write_eigvec_sectors
!         : atomic_write_atomcix_sectors
! source  : atomic_write.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : write output files
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> write basis of fullspace to file 'atom.basis.dat'
subroutine atomic_write_basis()
    use constants,         only: mytmp
    use control,           only: ncfgs
    use m_basis_fullspace, only: dec_basis, index_basis, bin_basis

    implicit none

    ! local variables
    integer :: i

    ! open file 'atom.basis.dat' to write
    open(mytmp, file='atom.basis.dat')
    ! write the header
    write(mytmp, '(a)') '#      i |  decimal(i) |    index(i) |      binary(i) |'
    do i=1, ncfgs
        write(mytmp, "(I10,4X,I10,4X,I10,8X,14I1)") i, dec_basis(i), index_basis(dec_basis(i)), bin_basis(:,i)   
    enddo 
    close(mytmp)

    return
end subroutine atomic_write_basis

!>>> write eigenvalue of fullspace to file 'atom.eigval.dat'
subroutine atomic_write_eigval_fullspace()
    use constants,         only: mytmp
    use control,           only: ncfgs
    use m_glob_fullspace,  only: hmat_eigval, occu_mat

    implicit none

    ! local variables
    integer :: i
    
    ! open file 'atom.eigval.dat' to write
    open(mytmp, file='atom.eigval.dat')
    write(mytmp, '(a)') '#       i |     eigenvalue(i) |      occupancy(i) |'
    do i=1, ncfgs
        write(mytmp, "(I10, 2F20.10)") i, hmat_eigval(i), occu_mat(i,i)
    enddo
    close(mytmp)

    return
end subroutine atomic_write_eigval_fullspace

!>>> write eigenvector of fullspace to file 'atom.eigvec.dat'
subroutine atomic_write_eigvec_fullspace()
    use constants,         only: mytmp, eps6
    use control,           only: ncfgs
    use m_basis_fullspace, only: bin_basis
    use m_glob_fullspace,  only: hmat_eigvec

    implicit none

    ! local variables
    integer :: i, j
 
    ! open file 'atom.eigvec.dat' to write
    open(mytmp, file='atom.eigvec.dat')
    write(mytmp, '(a)') '#      i |       j |  eigenvector(i,j) |       fock_basis(i) |'
    do i=1, ncfgs
        do j=1, ncfgs
            if ( abs(hmat_eigvec(j,i)) < eps6 ) cycle
            write(mytmp, "(2I10, F20.10, 8X, 14I1)") j, i, hmat_eigvec(j,i), bin_basis(:,j) 
        enddo
    enddo 
    close(mytmp)

    return
end subroutine atomic_write_eigvec_fullspace

!>>> write atom.cix for CTQMC input
subroutine atomic_write_atomcix_fullspace()
    use constants,        only: mytmp, zero
    use control,          only: nband, norbs, ncfgs, isoc
    use m_glob_fullspace, only: hmat_eigval, occu_mat, anni_fmat

    implicit none

    ! local variables
    integer :: i, j, k
    integer :: s_order

    ! open file 'atom.cix' to write
    open(mytmp, file='atom.cix')

    ! write eigenvalues
    write(mytmp,'(a)') '# eigenvalues: index | energy | occupy | spin'
    do i=1,ncfgs
        write(mytmp,'(I10,3F20.10)') i, hmat_eigval(i), occu_mat(i,i), zero
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
                write(mytmp,'(3I10, F20.10)') k, j, i, anni_fmat(k, j, s_order)
            enddo 
        enddo 
    enddo 
    close(mytmp)

    return
end subroutine atomic_write_atomcix_fullspace

!>>> write eigenvalue of sectors to file 'atom.eigval.dat'
subroutine atomic_write_eigval_sectors()
    use constants,      only: mytmp
    use m_glob_sectors, only: nsectors, sectors

    implicit none

    ! local variables
    integer :: i, j
    integer :: counter

    ! open file 'atom.eigval.dat' to write
    open(mytmp, file='atom.eigval.dat')
    counter = 0
    write(mytmp, '(a)') "#      i |     sect(i) | electron(i) |           j |   eigenvalue(j,i) |"
    do i=1, nsectors
        do j=1, sectors(i)%ndim
            counter = counter + 1
            write(mytmp, "(I10,4X,I10,4X,I10,4X,I10, F20.10)") counter, i, sectors(i)%nelectron, j, sectors(i)%myeigval(j)
        enddo
    enddo
    close(mytmp)

    return
end subroutine atomic_write_eigval_sectors

!>>> write eigenvector of sectors to file 'atom.eigval.dat'
subroutine atomic_write_eigvec_sectors()
    use constants,         only: mytmp, eps6
    use m_basis_fullspace, only: bin_basis
    use m_glob_sectors,    only: nsectors, sectors

    implicit none

    ! local variables
    integer :: i, j, k
    integer :: counter

    open(mytmp, file="atom.eigvec.dat")
    counter = 0
    write(mytmp, '(a)') '#      i |       j |       k |       eigvec(j,k) |       fock_basis(j) |'
    do i=1, nsectors
        do j=1, sectors(i)%ndim
            do k=1, sectors(i)%ndim
                if ( abs(sectors(i)%myeigvec(k,j)) < eps6 ) cycle
                    write(mytmp, "(3I10, F20.10, 8X, 14I1)") i, k+counter, j+counter, &
                    sectors(i)%myeigvec(k,j), bin_basis(:,sectors(i)%mybasis(k))
            enddo
        enddo
        counter = counter + sectors(i)%ndim
    enddo
    close(mytmp)

    return
end subroutine atomic_write_eigvec_sectors

!>>> write atom.cix for CTQMC input, good quantum number algorithm
subroutine atomic_write_atomcix_sectors()
    use constants,      only: mytmp
    use control,        only: isoc
    use m_glob_sectors, only: nsectors, sectors, max_dim_sect, ave_dim_sect

    implicit none

    ! local variables
    integer :: i, j, k, m, n, ii
    integer :: s_order
 
    ! open 'atom.cix' to write
    open(mytmp, file='atom.cix')
    ! write header
    write(mytmp, "(a)") "# THIS FILE IS GENERATED BY JASMINE CODE, CAN BE USED BY PANSY AND MANJUSHAKA CTQMC CODE"
    write(mytmp, *)
    ! write number of sectors
    write(mytmp, "(a)") "#NUMBER OF SECTORS | MAXIMUM DIMENSION OF SECTORS | AVERAGE DIMENSION OF SECTORS"
    write(mytmp, "(I10,15X,I10, 15X, F20.10)") nsectors, max_dim_sect, ave_dim_sect 

    ! write dimension, total electrons, next_sector, eigenvalue of each sector
    do i=1, nsectors
        write(mytmp, "(a)") "#SECT_INFO: INDEX  |  NDIM  |  NELEC  |   NOPS  |  ISTART"  
        write(mytmp, "(4X,5I10)") i, sectors(i)%ndim, sectors(i)%nelectron, sectors(i)%nops, sectors(i)%istart

        ! write next_sector
        write(mytmp, "(4X,a)") "#NEXT_SECTOR   ANNIH     CREAT"
        do j=1, sectors(i)%nops
            ! adjust the orbital order for CTQMC, up, up, up, dn, dn, dn
            if (isoc==0) then
                if (j <= sectors(i)%nops / 2) then
                    s_order = 2*j-1  
                else
                    s_order = 2*(j - sectors(i)%nops / 2)
                endif
            else 
                s_order = j
            endif 
            write(mytmp, "(2X, 3I10)") j, sectors(i)%next_sector(s_order,0), sectors(i)%next_sector(s_order,1)  
        enddo

        ! write eigeanvalue
        write(mytmp, "(4X,a)") "#EIGENVALUE"
        do j=1, sectors(i)%ndim
            write(mytmp, "(2X,I10, F20.10)") j, sectors(i)%myeigval(j) 
        enddo
    enddo

    ! write fmat
    ! write header
    write(mytmp, '(a)') '# BEGIN FMAT '
    do i=1, nsectors
        do j=1, sectors(i)%nops
            if (isoc==0) then
                if (j <= sectors(i)%nops / 2) then
                    s_order = 2*j-1  
                else
                    s_order = 2*(j-sectors(i)%nops / 2)
                endif
            else 
                s_order = j
            endif 
            do k=0,1
                ii = sectors(i)%next_sector(s_order,k)
                if (ii == -1) cycle 
                !write(mytmp)  sectors(i)%myfmat(s_order,k)%item(:,:)
                do m=1, sectors(i)%myfmat(s_order,k)%m
                do n=1, sectors(i)%myfmat(s_order,k)%n
                    if ( abs(sectors(i)%myfmat(s_order,k)%item(n,m)) < 1e-10 ) cycle
                    write(mytmp, '(5I6, F16.10)') n, m, k, j, i, sectors(i)%myfmat(s_order,k)%item(n,m)
                enddo
                enddo
            enddo  ! over k={0,1} loop
        enddo ! over j={1, sectors(i)%nops} loop
    enddo  ! over i={1, nsect} loop
    close(mytmp)

    return
end subroutine atomic_write_atomcix_sectors

!>>> write the transformation matrix from the original basis to natural basis
subroutine atomic_write_natural(info)
    use constants, only: dp, mytmp
    use control,   only: norbs
    use m_spmat,   only: tran_umat

    ! external variables
    character(len=*), intent(in) :: info

    ! local variables
    integer :: i,j

    open(mytmp, file='atom.natural.dat')
    write(mytmp,'(a)') info
    write(mytmp,'(a)') '#      i |       j |    umat_real(i,j) |    umat_imag(i,j) |'
    do i=1, norbs
        do j=1, norbs
            write(mytmp, '(2I10,2F20.10)') j, i, tran_umat(j,i)
        enddo
    enddo 
    close(mytmp)

    return
end subroutine atomic_write_natural
