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
        call s_print_error('atomic_read_cf', 'no file atomic.cf.in !')
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
        call s_print_error('atomic_read_eimp', 'no file atomic.eimp.in !')
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
        call s_print_error('atomic_read_umat', 'no file atomic.umat.in')
    endif

    return
end subroutine atomic_read_umat
