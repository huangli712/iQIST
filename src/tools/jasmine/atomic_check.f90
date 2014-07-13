!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_check_config
!         : atomic_check_hmat_real
! source  : atomic_check.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : do some check
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> check the validity of input config parameters
subroutine atomic_check_config()
    use constants,    only: mystd, zero
    use control

    ! local variables
    logical :: lpass

    lpass = .true.

    ! check nband
    if (nband <= 0) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: number of bands must be larger than zero !'
        write(mystd, *)
        lpass = .false.
    endif

    ! check itask
    if (itask /= 1 .and. itask /= 2) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: itask must be 1 or 2 !'
        write(mystd, *)
        lpass = .false.
    endif

    ! check icf
    if (icf /= 0 .and. icf /=1 .and. icf /= 2) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: icf must be one of 0, 1, 2 !'
        write(mystd, *)
        lpass = .false.
    endif

    ! check isoc
    if (isoc /= 0 .and. isoc /=1) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: isoc must be 0 or 1 !'
        write(mystd, *)
        lpass = .false.
    endif
    if (isoc == 1 .and. nband /=3 .and. nband /= 5 .and. nband /=7) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: only support spin-orbital coupling for nband=3, 5, 7 !'
        write(mystd, *)
        lpass = .false.
    endif

    ! check icu
    if (icu /= 1 .and. icu /= 2) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: icu must be 1 or 2 !'
        write(mystd, *)
        lpass = .false.
    endif
    if (icu == 2 .and. nband /= 5 .and. nband /= 7) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: only support Slater-Cordon type Coulomb interaction for nband=5, 7 !'
        write(mystd, *)
        lpass = .false.
    endif


    ! check Uc, Uv, Jz, Js, Jp, Ud, JH
    if (Uc < zero .or. Uv < zero .or. Jz < zero .or. Js < zero .or. Jp < zero .or. Ud < zero .or. JH < zero) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: Uc, Uv, Jz, Js, Jp, Ud, JH must be larger than zero !'
        write(mystd, *)
        lpass = .false.
    endif

    ! check ictqmc
    if (ictqmc /=1 .and. ictqmc /=2 .and. ictqmc /=3 .and. ictqmc /=4 .and. ictqmc /=5) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: ictqmc must be one of 1, 2, 3, 4, 5 !'
        write(mystd, *)
        lpass = .false.
    endif
    if (ictqmc == 3 .and. isoc == 1) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: CTQMC good quantum numbers (N,Sz) algorithm is NOT&
                                 supported for spin-orbital coupling case ! '
        write(mystd, *)
        lpass = .false.
    endif
    if (ictqmc == 4 .and. isoc == 1) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: CTQMC good quantum numbers (N,Sz,Ps) algorithm is NOT &
                               supported for spin-orbital coupling case ! '
        write(mystd, *)
        lpass = .false.
    endif
    if (ictqmc == 4 .and. icu == 2) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: CTQMC good quantum numbers (N,Sz,Ps) algorithm is NOT &
                                 supported for Slater-Cordon type Coulomb interaction U'
        write(mystd, *)
        lpass = .false.
    endif
    if (ictqmc == 5 .and. isoc == 0) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: CTQMC good quantum numbers (N,Jz) algorithm is ONLY &
                               supported for spin-orbital coupling case !'
        write(mystd, *)
        lpass = .false.
    endif
    if (ictqmc == 5 .and. isoc == 1 .and. icf /= 0) then
        write(mystd, '(2X,a)') 'jasmine >>> ERROR: CTQMC good quantum numbers (N,Jz) algorithm is NOT &
                               supported for spin-orbital coupling plus crystal field case !'
        write(mystd, *)
        lpass = .false.
    endif
   
    if (lpass == .true.) then
        write(mystd, '(2X,a)') 'jasmine >>> Good News: all control parameters are OK !'
    else
        call atomic_print_error('atomic_check_config', 'Found some wrong setting of parameters, please check the "atom.config.in" file !')
    endif

    return
end subroutine atomic_check_config

!>>> check whether Hamiltonian is real
subroutine atomic_check_hmat_real(ndim, hmat, lreal)
    use constants, only: dp, eps6

    implicit none

    ! external variables
    ! dimension of the Hamiltonian
    integer, intent(in) :: ndim
    ! the Hamiltonian matrix
    complex(dp), intent(in) :: hmat(ndim, ndim)
    ! whether Hamiltonian is real
    logical, intent(out) :: lreal

    ! local variables
    integer :: i, j

    do i=1, ndim
        do j=1, ndim
            if ( aimag(hmat(j,i)) > eps6 ) then
                lreal = .false. 
                return
            endif
        enddo
    enddo

    lreal = .true.

    return
end subroutine atomic_check_hmat_real
