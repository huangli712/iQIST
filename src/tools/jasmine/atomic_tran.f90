!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_make_umat_c2r
!         : atomic_make_umat_r2c
!         : atomic_make_umat_c2j
!         : atomic_tran_represent
!         : atomic_tran_represent_real
!         : atomic_tran_cumat
! source  : atomic_natural.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : make transformation from one space to another space 
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> make transformation matrix from complex orbital basis to real orbital 
subroutine atomic_make_umat_c2r( umat_c2r )
    use constants, only: dp, czero, cone, czi
    use control,   only: nband, norbs

    implicit none

    ! external variables
    ! the transformation matrix from real orbitals to complex orbitals
    complex(dp), intent(out) :: umat_c2r( norbs, norbs )

    ! sqrt(2)
    real(dp) :: sqrt2

    sqrt2 = sqrt(2.0_dp)
    umat_c2r = czero

    if ( nband == 3) then
    ! the real orbital order is:  
    ! dxzup, dxzdn, dyzup, dyzdn, dxyup, dxydn
    ! the complex orbital |Lz,Sz> order is
    ! -1up, -1dn, 0up, 0dn, 1up, 1dn 
        umat_c2r(1,1) =  czi/sqrt2
        umat_c2r(5,1) =  czi/sqrt2
        umat_c2r(2,2) =  czi/sqrt2
        umat_c2r(6,2) =  czi/sqrt2
        umat_c2r(1,3) =  cone/sqrt2
        umat_c2r(5,3) =  -cone/sqrt2
        umat_c2r(2,4) =  cone/sqrt2
        umat_c2r(6,4) =  -cone/sqrt2
        umat_c2r(3,5) =  cone
        umat_c2r(4,6) =  cone

    elseif ( nband == 5 ) then
    ! the real orbital order is: 
    ! dz2up, dz2dn, dxzup, dxzdn, dyzup, dyzdn, dx2-y2up, dx2-y2dn, dxyup, dxydn 
    ! the complex orbital |Lz,Sz> order is:
    ! -2up, -2dn, -1up, -1dn, 0up, 0dn, 1up, 1dn, 2up, 2dn
        umat_c2r(5,1) = cone
        umat_c2r(6,2) = cone
        umat_c2r(3,3) =  cone/sqrt2 
        umat_c2r(7,3) =  -cone/sqrt2
        umat_c2r(4,4) =  cone/sqrt2 
        umat_c2r(8,4) =  -cone/sqrt2
        umat_c2r(3,5) =    czi/sqrt2
        umat_c2r(7,5) =    czi/sqrt2
        umat_c2r(4,6) =    czi/sqrt2
        umat_c2r(8,6) =    czi/sqrt2
        umat_c2r(1,7) =  cone/sqrt2
        umat_c2r(9,7) =  cone/sqrt2
        umat_c2r(2,8) =  cone/sqrt2
        umat_c2r(10,8) =  cone/sqrt2
        umat_c2r(1,9) =    czi/sqrt2
        umat_c2r(9,9) =   -czi/sqrt2
        umat_c2r(2,10) =    czi/sqrt2
        umat_c2r(10,10) =  -czi/sqrt2
    elseif (nband == 7) then
    ! the real orbital order is:
    ! fz3up, fz3dn, fxz2up, fxz2dn, fyz2up, fyz2dn, fz(x2-y2)up, fz(x2-y2)dn, fxyzup, fxyzdn,
    ! fx(x2-3y2)up, fx(x2-3y2)dn, fy(3x2-y2)up, fy(3x2-y2)dn
    ! the complex orbital order is:
    ! -3up, -3dn, -2up, -2dn, -1up, -1dn, 0up, 0dn, 1up, 1dn, 2up, 2dn, 3up, 3dn    
        umat_c2r( 7, 1) = cone 
        umat_c2r( 8, 2) = cone 
        umat_c2r( 5, 3) = cone/sqrt2
        umat_c2r( 9, 3) = -cone/sqrt2
        umat_c2r( 6, 4) = cone/sqrt2
        umat_c2r(10, 4) = -cone/sqrt2
        umat_c2r( 5, 5) = czi/sqrt2
        umat_c2r( 9, 5) = czi/sqrt2
        umat_c2r( 6, 6) = czi/sqrt2
        umat_c2r(10, 6) = czi/sqrt2
        umat_c2r( 3, 7) = cone/sqrt2
        umat_c2r(11, 7) = cone/sqrt2
        umat_c2r( 4, 8) = cone/sqrt2
        umat_c2r(12, 8) = cone/sqrt2
        umat_c2r( 3, 9) = czi/sqrt2
        umat_c2r(11, 9) = -czi/sqrt2
        umat_c2r( 4,10) = czi/sqrt2
        umat_c2r(12,10) = -czi/sqrt2
        umat_c2r( 1,11) = cone/sqrt2
        umat_c2r(13,11) = -cone/sqrt2
        umat_c2r( 2,12) = cone/sqrt2
        umat_c2r(14,12) = -cone/sqrt2
        umat_c2r( 1,13) = czi/sqrt2
        umat_c2r(13,13) = czi/sqrt2
        umat_c2r( 2,14) = czi/sqrt2
        umat_c2r(14,14) = czi/sqrt2
    else
        call atomic_print_error('atomic_make_umat_c2r', 'not implemented for this nband!')
    endif

    return
end subroutine atomic_make_umat_c2r

!>>> make umat from real orbital basis to complex orbital basis
subroutine atomic_make_umat_r2c(umat_r2c)
    use constants, only: czero
    use control,   only: norbs

    ! external variables
    complex(dp), intent(out) :: umat_r2c(norbs, norbs)
 
    ! local variables
    complex(dp) :: umat_c2r(norbs, norbs)

    umat_c2r = czero
    call atomic_make_umat_c2r(umat_c2r)

    umat_r2c = transpose(dconjg(umat_c2r))

    return
end subroutine atomic_make_umat_r2c

!>>> make CG coefficients
subroutine atomic_make_umat_c2j( umat_c2j ) 
    use constants,  only: dp, czero
    use control,    only: nband, norbs

    implicit none
    
    ! the transformation matrix from complex orbitals |lz,sz> to |j2,jz>
    complex(dp), intent(out) :: umat_c2j( norbs, norbs )
    
    umat_c2j = czero

    if (nband == 3) then
    ! the |lz,sz> order is:
    ! |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>
    ! the |j2,jz> order is:
    ! |1/2,-1/2>, |1/2,1/2>, |3/2,-3/2>, |3/2, -1/2>, |3/2, 1/2>, |3/2,3/2>
        umat_c2j(1,1) = -sqrt(2.0/3.0) 
        umat_c2j(4,1) =  sqrt(1.0/3.0) 
        umat_c2j(3,2) = -sqrt(1.0/3.0) 
        umat_c2j(6,2) =  sqrt(2.0/3.0) 
        umat_c2j(2,3) =  1.0_dp
        umat_c2j(1,4) =  sqrt(1.0/3.0) 
        umat_c2j(4,4) =  sqrt(2.0/3.0) 
        umat_c2j(3,5) =  sqrt(2.0/3.0) 
        umat_c2j(6,5) =  sqrt(1.0/3.0) 
        umat_c2j(5,6) =  1.0_dp
    elseif ( nband == 5 ) then
    ! the |lz,sz> order is:
    ! |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>, |2,up>, |2,dn>
    ! the |j2,jz> order is:
    ! |3/2,-3/2>, |3/2,-1/2>, |3/2,1/2>, |3/2,3/2>
    ! |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
        umat_c2j(1,1) = -sqrt(4.0/5.0) 
        umat_c2j(4,1) =  sqrt(1.0/5.0) 
        umat_c2j(3,2) = -sqrt(3.0/5.0) 
        umat_c2j(6,2) =  sqrt(2.0/5.0) 
        umat_c2j(5,3) = -sqrt(2.0/5.0) 
        umat_c2j(8,3) =  sqrt(3.0/5.0) 
        umat_c2j(7,4) = -sqrt(1.0/5.0) 
        umat_c2j(10,4)=  sqrt(4.0/5.0) 
        umat_c2j(2,5) = 1.0_dp 
        umat_c2j(1,6) = sqrt(1.0/5.0)
        umat_c2j(4,6) = sqrt(4.0/5.0)
        umat_c2j(3,7) = sqrt(2.0/5.0)
        umat_c2j(6,7) = sqrt(3.0/5.0)
        umat_c2j(5,8) = sqrt(3.0/5.0)
        umat_c2j(8,8) = sqrt(2.0/5.0)
        umat_c2j(7,9) = sqrt(4.0/5.0)
        umat_c2j(10,9)= sqrt(1.0/5.0)
        umat_c2j(9,10)= 1.0_dp
    elseif ( nband == 7 ) then
    ! the |lz,sz> order is:
    ! |-3,up>, |-3,dn>, |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>, |2,up>, |2,dn>, |3,up>, |3,dn>
    ! the |j2,jz> order is:
    ! |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
    ! |7/2,-7/2>, |7/2,-5/2>, |7/2,-3/2>, |7/2,-1/2>, |7/2,1/2>, |7/2,3/2>, |7/2,5/2>, |7/2, 7/2>
        umat_c2j(1, 1) = -sqrt(6.0/7.0)
        umat_c2j(4, 1) =  sqrt(1.0/7.0)
        umat_c2j(3, 2) = -sqrt(5.0/7.0)
        umat_c2j(6, 2) =  sqrt(2.0/7.0)
        umat_c2j(5, 3) = -sqrt(4.0/7.0)
        umat_c2j(8, 3) =  sqrt(3.0/7.0)
        umat_c2j(7, 4) = -sqrt(3.0/7.0)
        umat_c2j(10,4) =  sqrt(4.0/7.0)
        umat_c2j(9, 5) = -sqrt(2.0/7.0)
        umat_c2j(12,5) =  sqrt(5.0/7.0)
        umat_c2j(11,6) = -sqrt(1.0/7.0)
        umat_c2j(14,6) =  sqrt(6.0/7.0)

        umat_c2j(2, 7)  = 1.0_dp
        umat_c2j(1, 8) =  sqrt(1.0/7.0)
        umat_c2j(4, 8) =  sqrt(6.0/7.0)
        umat_c2j(3, 9) =  sqrt(2.0/7.0)
        umat_c2j(6, 9) =  sqrt(5.0/7.0)
        umat_c2j(5,10) =  sqrt(3.0/7.0)
        umat_c2j(8,10) =  sqrt(4.0/7.0)
        umat_c2j(7,11) =  sqrt(4.0/7.0)
        umat_c2j(10,11)=  sqrt(3.0/7.0)
        umat_c2j(9,12) =  sqrt(5.0/7.0)
        umat_c2j(12,12)=  sqrt(2.0/7.0)
        umat_c2j(11,13)=  sqrt(6.0/7.0)
        umat_c2j(14,13)=  sqrt(1.0/7.0)
        umat_c2j(13,14)=  1.0_dp
    else
        call atomic_print_error('atomic_make_umat_c2j','not implemented !')
    endif

    return
end subroutine atomic_make_umat_c2j

!>>> transformation from one representation to another representation
subroutine atomic_tran_represent( ndim, amat, umat )
    use constants, only: dp

    implicit none

    ! external variables
    integer, intent(in) :: ndim
    complex(dp), intent(inout) :: amat(ndim, ndim)
    complex(dp), intent(in) :: umat(ndim, ndim)

    ! local variables
    complex(dp) :: tmp_mat(ndim, ndim)

    call zmat_zgemm0( ndim, amat, umat, tmp_mat )
    call zmat_zgemm2( ndim, umat, tmp_mat, amat )      

    return
end subroutine atomic_tran_represent

!>>> transformation from one representation to another representation, real version
subroutine atomic_tran_represent_real( ndim, amat, umat )
    use constants, only: dp

    implicit none

    ! external variables
    integer, intent(in) :: ndim
    real(dp), intent(inout) :: amat(ndim, ndim)
    real(dp), intent(in) :: umat(ndim, ndim)

    ! local variables
    real(dp) :: tmp_mat(ndim, ndim)

    call dmat_dgemm0( ndim, amat, umat, tmp_mat )
    call dmat_dgemm1( ndim, umat, tmp_mat, amat )      

    return
end subroutine atomic_tran_represent_real

!>>> transform Coulomb interaction U tensor 
subroutine atomic_tran_cumat(amtrx, cumat, cumat_t)
    use constants, only: dp, czero, epst
    use control,   only: norbs

    implicit none

    ! transformation matrix from orginal basis to natural basis
    complex(dp), intent(in) :: amtrx(norbs, norbs)
    ! coefficents matrix for generalized interaction U in orginal basis
    complex(dp), intent(in) :: cumat(norbs, norbs, norbs, norbs)
    ! coefficents matrix for generalized interaction U in natural basis
    complex(dp), intent(out) :: cumat_t(norbs, norbs, norbs, norbs)

    ! local varoables
    ! loop index over orbits in orginal single particle basis
    integer :: alpha1, alpha2
    integer :: alpha3, alpha4
    integer :: sigma1, sigma2
    integer :: sigma3, sigma4
    ! auxiliary complex(dp) variables
    complex(dp) :: ctmp

    ! initialize cumat_t to be zero
    cumat_t = czero 

    sigma1loop: do sigma1=1,norbs
    sigma2loop: do sigma2=1,norbs
    sigma3loop: do sigma3=1,norbs
    sigma4loop: do sigma4=1,norbs
        ctmp = czero

        alpha1loop: do alpha1=1,norbs
        alpha2loop: do alpha2=1,norbs
        alpha3loop: do alpha3=1,norbs
        alpha4loop: do alpha4=1,norbs
            if (abs(cumat(alpha1, alpha2, alpha3, alpha4)) .lt. epst) cycle
            ctmp = ctmp + cumat(alpha1, alpha2, alpha3, alpha4)          &
                 * conjg(amtrx(alpha1, sigma1)) * amtrx(alpha3, sigma3)  &
                 * conjg(amtrx(alpha2, sigma2)) * amtrx(alpha4, sigma4)
        enddo alpha4loop ! over alpha4={1,norbs} loop
        enddo alpha3loop ! over alpha3={1,norbs} loop
        enddo alpha2loop ! over alpha2={1,norbs} loop
        enddo alpha1loop ! over alpha1={1,norbs} loop

        cumat_t(sigma1, sigma2, sigma3, sigma4) = ctmp
    enddo sigma4loop ! over sigma4={1,norbs} loop
    enddo sigma3loop ! over sigma3={1,norbs} loop
    enddo sigma2loop ! over sigma2={1,norbs} loop
    enddo sigma1loop ! over sigma1={1,norbs} loop

    return
end subroutine atomic_tran_cumat
