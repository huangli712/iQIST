!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_make_natural
!           atomic_2natural_case1
!           atomic_2natural_case2
!           atomic_2natural_case3
!           atomic_2natural_case4
!           atomic_mat_2nospin
!           atomic_mat_2spin
! source  : atomic_natural.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : make natural basis 
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> make natural basis, the natural basis is the basis on which the 
! impurity energy matrix is diagonal
subroutine atomic_make_natural()
    use constants, only: dp, czero
    use control,   only: norbs, itask, icf, isoc, icu
    use m_spmat,   only: cumat, tran_umat

    implicit none

    ! local variables
    complex(dp) :: tmp_mat(norbs, norbs, norbs, norbs)
    ! umat from real orbital basis to complex orbital basis
    complex(dp) :: umat_r2c(norbs, norbs)
    ! umat from complex orbital basis to real orbital basis
    complex(dp) :: umat_c2r(norbs, norbs)

    tmp_mat  = czero
    umat_r2c = czero
    umat_c2r = czero

    ! make umat from origional basis to natural basis: tran_umat
    ! and set the eimp: eimpmat 
    if ( itask==1) then
        if ( isoc==0 .and. (icf==0 .or. icf==1) ) then
            ! for model calculation, no spin-orbital coupling, no crystal field or
            ! crystal field is diagonal, the real orbital basis is the natural basis
            call atomic_2natural_case1()
        elseif ( isoc==0 .and. icf==2 ) then
            call atomic_2natural_case2() 
        elseif ( isoc==1 .and. icf==0 ) then
            call atomic_2natural_case3()
        elseif ( isoc==1 .and. icf>0 ) then
            call atomic_2natural_case4()
        endif
    endif

    ! we need transform Coulomb interaction U
    ! for non-soc case, the tran_umat is defined as from real orbital basis to natural basis
    if (isoc==0) then
        ! for Slater-Cordon parameters Coulomb interaction U
        ! we first need to transfrom cumat from complex orbital basis to real orbital basis
        if ( icu == 2 ) then
            call atomic_make_umat_r2c( umat_r2c )
            umat_c2r = transpose( dconjg(umat_r2c) )
            call atomic_tran_cumat( umat_c2r, cumat, tmp_mat )
            cumat = tmp_mat
        endif
    ! for soc case, the tran_umat is defined as from complex orbital basis to natural basis
    else
        ! for Kanamori parameters Coulomb interaction U
        ! we first need to transfrom cumat from real orbital basis to complex orbital basis
        if ( icu == 2 ) then
            call atomic_make_umat_r2c( umat_r2c )
            call atomic_tran_cumat( umat_r2c, cumat, tmp_mat )
            cumat = tmp_mat
        endif
    endif

    ! finally, transform cumat to natural basis
    call atomic_tran_cumat(tran_umat, cumat, tmp_mat) 
    cumat = tmp_mat

    return
end subroutine atomic_make_natural

!>>> make natural basis for no crystal or diagonal crystal without spin-orbital coupling
subroutine atomic_2natural_case1()
    use constants, only: zero, cone
    use control,   only: norbs
    use m_spmat,   only: cfmat, eimpmat, tran_umat

    implicit none

    ! local variables
    integer :: i

    ! set eimpmat
    eimpmat = cfmat
    ! for this case, the natural basis is the real orbital basis
    ! so, the tran_umat is a unity matrix
    tran_umat = zero
    do i=1, norbs
        tran_umat(i,i) = cone
    enddo
    
    return
end subroutine atomic_2natural_case1

!>>> make natural basis for non-diagonal crystal field without spin-orbital coupling
subroutine atomic_2natural_case2()
    use constants, only: dp, czero
    use control,   only: norbs, nband
    use m_spmat,   only: cfmat, eimpmat, tran_umat

    implicit none

    ! local variables
    ! eimp matrix for no spin freedom
    complex(dp) :: eimp_nospin(nband, nband)
    real(dp) :: eimp_nospin_real(nband, nband)
    ! umat for no spin freedom
    complex(dp) :: umat_nospin(nband, nband)
    ! eigenvalue
    real(dp) :: eigval(nband)
    ! eigen vector
    real(dp) :: eigvec(nband, nband)
    ! index loop
    integer :: i

    ! set eimpmat
    eimpmat = cfmat   

    ! get eimp for nospin freedom
    call atomic_mat_2nospin(norbs, eimpmat, eimp_nospin)

    ! diagonalize eimp_nospin to get natural basis
    eimp_nospin_real = real(eimp_nospin)
    call dmat_dsyev(nband, nband, eimp_nospin_real, eigval, eigvec)

    eimp_nospin = czero
    do i=1, nband
        eimp_nospin(i,i) = eigval(i) 
    enddo
    umat_nospin = eigvec

    call atomic_mat_2spin(nband, eimp_nospin, eimpmat) 
    call atomic_mat_2spin(nband, umat_nospin, tran_umat)

    return
end subroutine atomic_2natural_case2

!>>> make natural basis for the case without crystal field and with spin-orbital coupling
! for this special case, the natural basis is |J^2,Jz>
subroutine atomic_2natural_case3()
    use constants, only: dp
    use control,   only: norbs
    use m_spmat,   only: eimpmat, socmat, tran_umat

    implicit none

    ! local variables
    ! umat from complex orbital basis to |j^2, jz> basis
    complex(dp) :: umat_c2j(norbs, norbs)

    ! set eimpmat
    eimpmat = socmat   

    call atomic_make_umat_c2j(umat_c2j)

    ! for soc case, the tran_umat is from complex orbital basis to natural basis
    tran_umat = umat_c2j

    ! transform sp_eimp_mat to natural basis
    call atomic_tran_represent(norbs, eimpmat, tran_umat)   

    return
end subroutine atomic_2natural_case3

!>>> make natural basis for the case with crystal field and with spin-orbital coupling
subroutine atomic_2natural_case4()
    use constants, only: dp
    use control,   only: norbs 
    use m_spmat,  only: cfmat, socmat, eimpmat, tran_umat

    implicit none

    ! local variables
    ! umat from real orbital basis to complex orbital basis
    complex(dp) :: umat_r2c(norbs, norbs)
    ! umat from complex orbital basis to natural basis
    complex(dp) :: umat_c2n(norbs, norbs)
    ! real version of eimp_mat
    real(dp) :: tmp_mat(norbs, norbs)
    ! eigen value 
    real(dp) :: eigval(norbs)
    ! eigen vector
    real(dp) :: eigvec(norbs, norbs)

    ! get umat_r2c
    call atomic_make_umat_r2c(umat_r2c)

    ! transfrom cfmat and to complex orbital basis
    call atomic_tran_represent(norbs, cfmat, umat_r2c)

    ! set eimpmat
    eimpmat = socmat + cfmat   

    tmp_mat = real(eimpmat)

    call dmat_dsyev(norbs, norbs, tmp_mat, eigval, eigvec)

    umat_c2n = eigvec
    tran_umat = umat_c2n

  !  transform eimpmat to natural basis
    call atomic_tran_represent(norbs, eimpmat, umat_c2n)   
    
    return
end subroutine atomic_2natural_case4

!>>> utility subroutines 
subroutine atomic_mat_2nospin(norbs, amat, bmat)
    use constants, only: dp
    
    implicit none

    ! external variables
    ! number of orbitals
    integer, intent(in) :: norbs
    ! matrix with spin
    complex(dp), intent(in) :: amat(norbs, norbs)
    ! matrix without spin
    complex(dp), intent(out) :: bmat(norbs/2, norbs/2)

    ! local variables
    ! loop index
    integer :: i, j

    do i=1, norbs/2
        do j=1, norbs/2
            bmat(j,i) = amat(2*j-1,2*i-1)
        enddo
    enddo

    return
end subroutine atomic_mat_2nospin 

subroutine atomic_mat_2spin(nband, amat, bmat)
    use constants, only: dp
    
    implicit none

    ! external variables
    ! number of orbitals
    integer, intent(in) :: nband
    ! matrix with spin
    complex(dp), intent(in) :: amat(nband, nband)
    ! matrix without spin
    complex(dp), intent(out) :: bmat(2*nband, 2*nband)

    ! local variables
    ! loop index
    integer :: i, j

    do i=1, nband
        do j=1, nband
            bmat(2*j-1,2*i-1) = amat(j,i)
            bmat(2*j,2*i)     = amat(j,i)
        enddo
    enddo

    return
end subroutine atomic_mat_2spin 
