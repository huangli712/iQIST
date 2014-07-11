!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_read_cf
!         : atomic_read_eimp
!         : atomic_read_umat
! source  : atomic_natural.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : print information
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> real crystal field from file 'atomic.cf.in'
subroutine atomic_read_cf()
    use constants, only: mytmp, dp, zero
    use m_spmat,   only: cfmat

    implicit none

    ! local variables
    ! file status
    logical :: exists
    ! iostat
    integer :: ierr
    ! dummy variables
    integer :: i, j
    real(dp) :: r1

    ! we read crystal field from file "atom.cf.in"
    ! inquire file
    inquire(file='atom.cf.in', exist=exists)

    if (exists .eqv. .true.) then
        open(mytmp, file='atom.cf.in')
        do while(.true.)
            read(mytmp, *, iostat=ierr) i, j, r1
            ! crystal field is actually real
            cfmat(i,j) = dcmplx(r1, zero)
            if (ierr /= 0) exit
        enddo
    else
        call atomic_print_error('atomic_read_cf', 'no file atomic.cf.in !')
    endif 

    return
end subroutine atomic_read_cf

!>>> read eimp from file 'atomic.eimp.in'
subroutine atomic_read_eimp()
    use constants,  only: mytmp, dp, zero
    use control,    only: norbs
    use m_spmat,    only: eimpmat

    implicit none

    ! local variables
    ! file status
    logical :: exists
    ! loop index
    integer :: i
    ! dummy variables
    integer :: i1, i2
    real(dp) :: r1

    ! we read eimp from file 'atomic.eimp.in'
    inquire(file='atom.eimp.in', exist=exists)

    if (exists .eqv. .true.) then
        open(mytmp, file='atom.eimp.in')
        do i=1, norbs
            read(mytmp, *) i1, i2, r1
            ! eimpmat is actually real in natural basis
            eimpmat(i,i) = dcmplx(r1, zero)
        enddo 
    else
        call atomic_print_error('atomic_read_eimp', 'no file atomic.eimp.in !')
    endif

    return
end subroutine atomic_read_eimp

!>>> read tran_umat from file 'atomic.umat.in'
subroutine atomic_read_umat()
    use constants, only: mytmp, dp, zero
    use control,   only: norbs
    use m_spmat,   only: tran_umat

    implicit none

    ! local variables
    ! file status
    logical :: exists
    ! loop index
    integer :: i, j
    ! dummy variables
    integer :: i1, i2
    real(dp) :: r1

    ! we read ran_umat from file 'atomic.umat.in'
    inquire(file='atom.umat.in', exist=exists)

    if (exists .eqv. .true.) then
        open(mytmp, file='atom.umat.in')
        do i=1, norbs
            do j=1, norbs
                read(mytmp, *) i1, i2, r1
                tran_umat(j,i) = dcmplx(r1, zero)
            enddo
        enddo
    else
        call atomic_print_error('atomic_read_umat', 'no file atomic.umat.in')
    endif

    return
end subroutine atomic_read_umat

!>>> read 'atom.sector.in' and 'atom.fmat.in', test purpose
subroutine test_read_sector_fmat()
    use constants
    use control
    use m_sector
    use m_glob_sectors
    implicit none

    ! local variables
    integer :: nsect
    type(t_sector), allocatable :: my_sectors(:)
    integer :: i,j,k,ii, row, col
    integer :: i1

    open(mytmp, file='atom.sector.in')
    read(mytmp,*) ! skip header
    read(mytmp,*) nsect
    ! allocate memory for my_sectors
    allocate(my_sectors(nsect))
    ! nullify each sector one by one
    do i=1, nsect
        call nullify_one_sector(my_sectors(i))
    enddo
    ! read information of each sector
    do i=1, nsect
        read(mytmp,*) ! skip header
        read(mytmp,*) i1, my_sectors(i)%ndim, my_sectors(i)%nelectron, my_sectors(i)%nops, my_sectors(i)%istart 
        call alloc_one_sector(my_sectors(i))
        ! read next_sector
        read(mytmp,*) ! skip header
        do j=1, my_sectors(i)%nops
            read(mytmp,*) i1, my_sectors(i)%next_sector(j,0), my_sectors(i)%next_sector(j,1) 
        enddo
        ! read eigenvalue
        read(mytmp,*) ! skip header
        do j=1, my_sectors(i)%ndim
            read(mytmp,*) i1, my_sectors(i)%myeigval(j)
        enddo
    enddo
    close(mytmp)
    ! read fmat
    open(mytmp, file='atom.fmat.in', form='unformatted')
    do i=1, nsect
        do j=1, my_sectors(i)%nops
            do k=0,1
                ii = my_sectors(i)%next_sector(j,k)
                if (ii == -1) cycle 
                my_sectors(i)%myfmat(j,k)%n = my_sectors(ii)%ndim
                my_sectors(i)%myfmat(j,k)%m = my_sectors(i)%ndim
                call alloc_one_fmat(my_sectors(i)%myfmat(j,k))
                read(mytmp)  my_sectors(i)%myfmat(j,k)%item(:,:)
            enddo  ! over k={0,1} loop
        enddo ! over j={1, sectors(i)%nops} loop
    enddo  ! over i={1, nsect} loop
    close(mytmp)

    ! check the result
    open(mytmp, file='test_fmat.dat')
    do i=1, nsect
        do j=1, my_sectors(i)%nops
        do k=0,1
            do row=1, my_sectors(i)%myfmat(j,k)%n 
            do col=1, my_sectors(i)%myfmat(j,k)%m 
                if( abs( my_sectors(i)%myfmat(j,k)%item(row,col) - sectors(i)%myfmat(j,k)%item(row,col) ) > 1e-10 )then
                    write(mytmp,'(5I5, 2F20.14)') i, j, k, row, col, my_sectors(i)%myfmat(j,k)%item(row,col), sectors(i)%myfmat(j,k)%item(row,col)
                endif 
            enddo
            enddo
        enddo 
        enddo
    enddo

    close(mytmp) 
end subroutine test_read_sector_fmat
