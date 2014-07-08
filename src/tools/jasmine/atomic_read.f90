!>>> real crystal field from file 'atomic.cf.in'
subroutine atomic_read_cf()
    use constants
    use m_spmat

    implicit none

    ! local variables
    ! file status
    logical :: exists
    ! iostat
    integer :: ierr
    ! dummy variables
    integer :: i, j
    real(dp) :: r1

    ! we read crystal field from file "atomic.cf.in"
    ! inquire file
    inquire(file='atomic.cf.in', exist=exists)

    if (exists .eqv. .true.) then
        open(mytmp, file='atomic.cf.in')
        do while(.true.)
            read(mytmp, iostat=ierr) i, j, r1
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
    use constants
    use control
    use m_spmat

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
    inquire(file='atomic.eimp.in', exist=exists)

    if (exists .eqv. .true.) then
        open(mytmp, file='atomic.eimp.in')
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
    use constants
    use control
    use m_spmat

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
    inquire(file='atomic.umat.in', exist=exists)

    if (exists .eqv. .true.) then
        open(mytmp, file='atomic.umat.in')
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


