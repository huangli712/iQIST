  subroutine atomic_make_umat_r2c( umat_r2c )
     use constants
     use control
 
     implicit none

! external variables
! the transformation matrix from real orbitals to complex orbitals
     complex(dp), intent(out) :: umat_r2c( norbs, norbs )

! local variables
     complex(dp) :: umat_c2r( norbs, norbs )

! sqrt(2)
     real(dp) :: sqrt2

     sqrt2 = sqrt(2.0_dp)

! case 1: if nband == 5, then t2g+eg orbitals will be included, the order is 
!         dz2up, dz2dn, dxzup, dxzdn, dyzup, dyzdn, dx2-y2up, dx2-y2dn, dxyup, dxydn 
! case 2: if nband == 3, then t2g orbitals will be included, the order is
!         dxzup, dxzdn, dyzup, dyzdn, dxyup, dxydn


     if ( nband == 5 ) then
! the complex orbital |Lz,Sz> order is:
! -2up, -2dn, -1up, -1dn, 0up, 0dn, 1up, 1dn, 2up, 2dn
         umat_r2c = czero
         umat_c2r = czero

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

         umat_r2c = transpose(dconjg(umat_c2r))  
     endif

     if ( nband == 3) then
! the complex orbital |Lz,Sz> order is
! -1up, -1dn, 0up, 0dn, 1up, 1dn 
         umat_r2c = czero
         umat_c2r = czero

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
  
         umat_r2c = transpose(dconjg(umat_c2r))  
     endif

     return
  end subroutine atomic_make_umat_r2c

!>>> make CG coefficients
  subroutine atomic_make_umat_c2j( umat_c2j ) 
     use constants
     use control
 
     implicit none

! external variables
! external variables
! symmetry index
! case 1: if isym == 5, then t2g+eg orbitals will be included, the order is 
!         dz2up, dz2dn, dxzup, dxzdn, dyzup, dyzdn, dx2-y2up, dx2-y2dn, dxyup, dxydn 
! case 2: if isym == 3, then t2g orbitals will be included, the order is
!         dxzup, dxzdn, dyzup, dyzdn, dxyup, dxydn

! the transformation matrix from complex orbitals |lz,sz> to |j2,jz>
     complex(dp), intent(out) :: umat_c2j( norbs, norbs )
      
     if ( nband == 5 ) then
!---------------------------------------------------------------------------
! |j2,jz> order is:
!---------------------------------------------------------------------------
! |3/2,-3/2>, |3/2,-1/2>, |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, 
! |3/2, 3/2>, |3/2, 1/2>, |5/2, 5/2>, |5/2, 3/2>, |5/2, 1/2>
!---------------------------------------------------------------------------
         umat_c2j(1,1) = dcmplx(-sqrt(four/five), zero)
         umat_c2j(4,1) = dcmplx( sqrt(one /five), zero)

         umat_c2j(3,2) = dcmplx(-sqrt(three/five), zero)
         umat_c2j(6,2) = dcmplx( sqrt(two/five), zero)

         umat_c2j(2,3) = cone

         umat_c2j(1,4) = dcmplx(sqrt(one/five), zero)
         umat_c2j(4,4) = dcmplx(sqrt(four/five), zero)

         umat_c2j(3,5) = dcmplx(sqrt(two/five), zero)
         umat_c2j(6,5) = dcmplx(sqrt(three/five), zero)

         umat_c2j(7,6) = dcmplx(-sqrt(one/five), zero)
         umat_c2j(10,6) = dcmplx(sqrt(four/five), zero)

         umat_c2j(5,7) = dcmplx(-sqrt(two/five), zero)
         umat_c2j(8,7) = dcmplx( sqrt(three/five), zero)

         umat_c2j(9,8) = cone

         umat_c2j(7,9) = dcmplx(sqrt(four/five), zero)
         umat_c2j(10,9) = dcmplx(sqrt(one/five), zero)

         umat_c2j(5,10) = dcmplx(sqrt(three/five), zero)
         umat_c2j(8,10) = dcmplx(sqrt(two/five), zero)
     endif
 
     return
  end subroutine atomic_make_umat_c2j

!>>> transformation from one representation to another representation
  subroutine atomic_tran_represent( ndim, amat, umat )
     use constant

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
