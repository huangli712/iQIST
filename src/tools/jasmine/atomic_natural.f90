!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_make_natural
!!!           atomic_2natural_case1
!!!           atomic_2natural_case2
!!!           atomic_2natural_case3
!!!           atomic_2natural_case4
!!!           atomic_mat_2nospin
!!!           atomic_mat_2spin
!!! source  : atomic_natural.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : make natural basis 
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> make natural basis, the natural basis is the basis on which the 
!!>>> impurity energy matrix is diagonal
  subroutine atomic_make_natural()
     use constants, only : dp, czero, mystd
     use control, only : norbs, itask, icf, isoc, icu
     use m_spmat, only : cumat, tran_umat
  
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
     else
         write(mystd, '(2X,a)') 'jasmine >>> natural basis is made outside !'
         write(mystd, *)
     endif
  
! dump eimpmat for reference
     call atomic_write_eimpmat()
  
! we need transform Coulomb interaction U
! for non-soc case, the tran_umat is defined as from real orbital basis to natural basis
     if (isoc==0) then
! for Slater-Cordon parameters Coulomb interaction U
! we first need to transfrom cumat from complex orbital basis to real orbital basis
         if ( icu == 2 ) then
             call atomic_mkumat_c2r( umat_c2r )
             call atomic_tran_cumat( umat_c2r, cumat, tmp_mat )
             cumat = tmp_mat
         endif
! for soc case, the tran_umat is defined as from complex orbital basis to natural basis
     else
! for Kanamori parameters Coulomb interaction U
! we first need to transfrom cumat from real orbital basis to complex orbital basis
         if ( icu == 1 ) then
             call atomic_mkumat_r2c( umat_r2c )
             call atomic_tran_cumat( umat_r2c, cumat, tmp_mat )
             cumat = tmp_mat
         endif
     endif
  
! finally, transform cumat to natural basis
     call atomic_tran_cumat(tran_umat, cumat, tmp_mat) 
     cumat = tmp_mat
  
     return
  end subroutine atomic_make_natural

!!>>> atomic_2natural_case1: make natural basis for no crystal or 
!!>>> diagonal crystal without spin-orbital coupling
  subroutine atomic_2natural_case1()
     use constants, only : mystd, zero, cone
     use control, only : norbs, mune
     use m_spmat, only : cfmat, eimpmat, tran_umat
  
     implicit none
  
! local variables
     integer :: i
  
! set eimpmat
     eimpmat = cfmat

! add chemical potential to eimpmat
     do i=1, norbs
         eimpmat(i,i) = eimpmat(i,i) + mune
     enddo
! for this case, the natural basis is the real orbital basis
! so, the tran_umat is a unity matrix
     tran_umat = zero
     do i=1, norbs
         tran_umat(i,i) = cone
     enddo
  
     write(mystd,'(2X,a)') 'jasmine >>> natural basis is: real orbital basis'
     write(mystd, *)
  
     call atomic_write_natural('# natural basis is real orbital, umat: real to natural')
   
     return
  end subroutine atomic_2natural_case1

!!>>> atomic_2natural_case2: make natural basis for non-diagonal 
!!>>> crystal field without spin-orbital coupling
  subroutine atomic_2natural_case2()
     use constants, only : mystd, dp, czero
     use control, only : norbs, nband, mune
     use m_spmat, only : cfmat, eimpmat, tran_umat
  
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
     call s_eig_sy(nband, nband, eimp_nospin_real, eigval, eigvec)
  
     eimp_nospin = czero
     do i=1, nband
         eimp_nospin(i,i) = eigval(i) 
     enddo
     umat_nospin = dcmplx(eigvec)
  
     call atomic_mat_2spin(nband, eimp_nospin, eimpmat) 
     call atomic_mat_2spin(nband, umat_nospin, tran_umat)
  
! add chemical potential to eimpmat
     do i=1, norbs
         eimpmat(i,i) = eimpmat(i,i) + mune
     enddo
  
     write(mystd, '(2X,a)') 'jasmine >>> natural basis is: linear combination of real orbitals '
     write(mystd, *)
  
     call atomic_write_natural('# natural basis is linear combination of real orbitals, umat: real to natural')
  
     return
  end subroutine atomic_2natural_case2

!!>>> atomic_2natural_case3: make natural basis for the case without 
!!>>> crystal field and with spin-orbital coupling
!!>>> for this special case, the natural basis is |J^2,Jz>
  subroutine atomic_2natural_case3()
     use constants, only : dp, mystd
     use control, only : norbs, mune
     use m_spmat, only : eimpmat, socmat, tran_umat
  
     implicit none
  
! local variables
! umat from complex orbital basis to |j^2, jz> basis
     complex(dp) :: umat_c2j(norbs, norbs)

! loop inex
     integer :: i
  
! set eimpmat
     eimpmat = socmat   
  
     call atomic_mkumat_c2j(umat_c2j)
  
! for soc case, the tran_umat is from complex orbital basis to natural basis
     tran_umat = umat_c2j
  
! transform sp_eimp_mat to natural basis
     call atomic_tran_repr_cmpl(norbs, eimpmat, tran_umat)   
  
! add chemical potential to eimpmat
     do i=1, norbs
         eimpmat(i,i) = eimpmat(i,i) + mune
     enddo
  
     write(mystd, '(2X,a)') 'jasmine >>> natural basis is: |j2,jz> '
     write(mystd, *)
  
     call atomic_write_natural('# natural basis is |j2,jz>, umat: complex to natural')
  
     return
  end subroutine atomic_2natural_case3

!!>>> atomic_2natural_case4: make natural basis for the case with 
!!>>> crystal field and with spin-orbital coupling
  subroutine atomic_2natural_case4()
     use constants, only : dp, mystd
     use control, only : norbs, mune

     use m_spmat, only : cfmat, socmat, eimpmat, tran_umat
  
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

! loop index
     integer :: i

! whether cfmat is real on complex orbital basis
     logical :: lreal
  
! get umat_r2c
     call atomic_mkumat_r2c(umat_r2c)
  
! transfrom cfmat to complex orbital basis
     call atomic_tran_repr_cmpl(norbs, cfmat, umat_r2c)
  
! check whether cfmat is real, if not, we cann't make natural basis
     call atomic_check_realmat(norbs, cfmat, lreal)
     if (lreal .eqv. .false.) then
         call s_print_error('atomic_2natural_case4', 'crystal field on &
             complex orbital basis is not real, cannot make natural basis !')
     endif
  
! set eimpmat
     eimpmat = socmat + cfmat   
  
     tmp_mat = real(eimpmat)
  
     call s_eig_sy(norbs, norbs, tmp_mat, eigval, eigvec)
  
     umat_c2n = eigvec
     tran_umat = umat_c2n
  
! transform eimpmat to natural basis
     call atomic_tran_repr_cmpl(norbs, eimpmat, umat_c2n)   
  
! add chemical poential to eimpmat
     do i=1, norbs
         eimpmat(i,i) = eimpmat(i,i) + mune
     enddo
  
     write(mystd, '(2X,a)') 'jasmine >>> natural basis is: linear combination of complex orbitals '
     write(mystd, *)
  
     call atomic_write_natural('# natural basis is linear combination of complex orbitals, umat: complex to natural')
  
     return
  end subroutine atomic_2natural_case4

!!>>> atomic_mat_2nospin: convert matrix with spin to no-spin
  subroutine atomic_mat_2nospin(norbs, amat, bmat)
     use constants, only : dp
     
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
  
!!>>> atomic_mat_2spin: convert matrix without spin to with spin
  subroutine atomic_mat_2spin(nband, amat, bmat)
      use constants, only : dp
      
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
!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_mkumat_c2r
!!!           atomic_mkumat_r2c
!!!           atomic_mkumat_c2j
!!!           atomic_tran_cumat
!!!           atomic_tran_repr_cmpl
!!!           atomic_tran_repr_real
!!! source  : atomic_natural.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : make transformation from one representation to another 
!!!           representation 
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_mkumat_c2r: make transformation matrix from 
!!>>> complex orbital basis to real orbital 
  subroutine atomic_mkumat_c2r( umat_c2r )
     use constants, only : dp, czero, cone, czi
     use control, only : nband, norbs
  
     implicit none
  
! external variables
! the transformation matrix from real orbitals to complex orbitals
     complex(dp), intent(out) :: umat_c2r( norbs, norbs )
  
! sqrt(2)
     real(dp) :: sqrt2
  
     sqrt2 = sqrt(2.0_dp)
     umat_c2r = czero
  
     if ( nband == 3 ) then
! the real orbital order (t2g) is:  
! dxzup, dxzdn, dyzup, dyzdn, dxyup, dxydn
! the corresponding p orbital order is:
! pyup, pydn, pxup, pxdn, pzup, pzdn
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
     elseif ( nband == 7 ) then
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
         call s_print_error('atomic_make_umat_c2r', &
                    'not implemented for this nband!')
     endif
  
     return
  end subroutine atomic_mkumat_c2r

!!>>> atomic_mkumat_r2c: make umat from real orbital 
!!>>> basis to complex orbital basis
  subroutine atomic_mkumat_r2c(umat_r2c)
     use constants, only : dp, czero
     use control, only : norbs
  
! external variables
     complex(dp), intent(out) :: umat_r2c(norbs, norbs)
   
! local variables
     complex(dp) :: umat_c2r(norbs, norbs)
  
     umat_c2r = czero
     call atomic_mkumat_c2r(umat_c2r)
  
     umat_r2c = transpose(dconjg(umat_c2r))
  
     return
  end subroutine atomic_mkumat_r2c

!!>>> atomic_mkumat_c2j: make CG coefficients
  subroutine atomic_mkumat_c2j( umat_c2j ) 
     use constants, only : dp, czero
     use control, only : nband, norbs
  
     implicit none
     
! the transformation matrix from complex orbitals |lz,sz> to |j2,jz>
     complex(dp), intent(out) :: umat_c2j( norbs, norbs )
     
     umat_c2j = czero
  
     if (nband == 3) then
! the |lz,sz> order is:
! |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>
! the |j2,jz> order is:
! |1/2,-1/2>, |1/2,1/2>, |3/2,-3/2>, |3/2, -1/2>, |3/2, 1/2>, |3/2,3/2>
         umat_c2j(1,1) = -sqrt(2.0_dp/3.0_dp) 
         umat_c2j(4,1) =  sqrt(1.0_dp/3.0_dp) 
         umat_c2j(3,2) = -sqrt(1.0_dp/3.0_dp) 
         umat_c2j(6,2) =  sqrt(2.0_dp/3.0_dp) 
         umat_c2j(2,3) =  1.0_dp
         umat_c2j(1,4) =  sqrt(1.0_dp/3.0_dp) 
         umat_c2j(4,4) =  sqrt(2.0_dp/3.0_dp) 
         umat_c2j(3,5) =  sqrt(2.0_dp/3.0_dp) 
         umat_c2j(6,5) =  sqrt(1.0_dp/3.0_dp) 
         umat_c2j(5,6) =  1.0_dp
     elseif ( nband == 5 ) then
! the |lz,sz> order is:
! |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>, |2,up>, |2,dn>
! the |j2,jz> order is:
! |3/2,-3/2>, |3/2,-1/2>, |3/2,1/2>, |3/2,3/2>
! |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
         umat_c2j(1,1) = -sqrt(4.0_dp/5.0_dp) 
         umat_c2j(4,1) =  sqrt(1.0_dp/5.0_dp) 
         umat_c2j(3,2) = -sqrt(3.0_dp/5.0_dp) 
         umat_c2j(6,2) =  sqrt(2.0_dp/5.0_dp) 
         umat_c2j(5,3) = -sqrt(2.0_dp/5.0_dp) 
         umat_c2j(8,3) =  sqrt(3.0_dp/5.0_dp) 
         umat_c2j(7,4) = -sqrt(1.0_dp/5.0_dp) 
         umat_c2j(10,4)=  sqrt(4.0_dp/5.0_dp) 
         umat_c2j(2,5) = 1.0_dp 
         umat_c2j(1,6) =  sqrt(1.0_dp/5.0_dp)
         umat_c2j(4,6) =  sqrt(4.0_dp/5.0_dp)
         umat_c2j(3,7) =  sqrt(2.0_dp/5.0_dp)
         umat_c2j(6,7) =  sqrt(3.0_dp/5.0_dp)
         umat_c2j(5,8) =  sqrt(3.0_dp/5.0_dp)
         umat_c2j(8,8) =  sqrt(2.0_dp/5.0_dp)
         umat_c2j(7,9) =  sqrt(4.0_dp/5.0_dp)
         umat_c2j(10,9)=  sqrt(1.0_dp/5.0_dp)
         umat_c2j(9,10)= 1.0_dp
     elseif ( nband == 7 ) then
! the |lz,sz> order is:
! |-3,up>, |-3,dn>, |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>, 
! | 0,dn>, | 1,up>, | 1,dn>, | 2,up>, | 2,dn>, | 3,up>, |3,dn>
! the |j2,jz> order is:
! |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2>, |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
! |7/2,-7/2>, |7/2,-5/2>, |7/2,-3/2>, |7/2,-1/2>, |7/2,1/2>, |7/2,3/2>, 
! |7/2,5/2>, |7/2, 7/2>
         umat_c2j(1, 1) = -sqrt(6.0_dp/7.0_dp)
         umat_c2j(4, 1) =  sqrt(1.0_dp/7.0_dp)
         umat_c2j(3, 2) = -sqrt(5.0_dp/7.0_dp)
         umat_c2j(6, 2) =  sqrt(2.0_dp/7.0_dp)
         umat_c2j(5, 3) = -sqrt(4.0_dp/7.0_dp)
         umat_c2j(8, 3) =  sqrt(3.0_dp/7.0_dp)
         umat_c2j(7, 4) = -sqrt(3.0_dp/7.0_dp)
         umat_c2j(10,4) =  sqrt(4.0_dp/7.0_dp)
         umat_c2j(9, 5) = -sqrt(2.0_dp/7.0_dp)
         umat_c2j(12,5) =  sqrt(5.0_dp/7.0_dp)
         umat_c2j(11,6) = -sqrt(1.0_dp/7.0_dp)
         umat_c2j(14,6) =  sqrt(6.0_dp/7.0_dp)
  
         umat_c2j(2, 7)  = 1.0_dp
         umat_c2j(1, 8) =  sqrt(1.0_dp/7.0_dp)
         umat_c2j(4, 8) =  sqrt(6.0_dp/7.0_dp)
         umat_c2j(3, 9) =  sqrt(2.0_dp/7.0_dp)
         umat_c2j(6, 9) =  sqrt(5.0_dp/7.0_dp)
         umat_c2j(5,10) =  sqrt(3.0_dp/7.0_dp)
         umat_c2j(8,10) =  sqrt(4.0_dp/7.0_dp)
         umat_c2j(7,11) =  sqrt(4.0_dp/7.0_dp)
         umat_c2j(10,11)=  sqrt(3.0_dp/7.0_dp)
         umat_c2j(9,12) =  sqrt(5.0_dp/7.0_dp)
         umat_c2j(12,12)=  sqrt(2.0_dp/7.0_dp)
         umat_c2j(11,13)=  sqrt(6.0_dp/7.0_dp)
         umat_c2j(14,13)=  sqrt(1.0_dp/7.0_dp)
         umat_c2j(13,14)=  1.0_dp
     else
         call s_print_error('atomic_make_umat_c2j','not implemented !')
     endif
  
     return
  end subroutine atomic_mkumat_c2j
 
!!>>> atomic_tran_cumat: transform Coulomb interaction U tensor 
!!>>> from one representation to another representation
  subroutine atomic_tran_cumat(amtrx, cumat, cumat_t)
     use constants, only : dp, czero, epst
     use control, only : norbs
  
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

!!>>> atomic_tran_repr_cmpl: transformation from one representation 
!!>>> to another representation, complex version
  subroutine atomic_tran_repr_cmpl( ndim, amat, umat )
     use constants, only : dp, cone, czero
  
     implicit none
  
! external variables
     integer, intent(in) :: ndim
     complex(dp), intent(inout) :: amat(ndim, ndim)
     complex(dp), intent(in) :: umat(ndim, ndim)
  
! local variables
     complex(dp) :: tmp_mat(ndim, ndim)
     complex(dp) :: alpha
     complex(dp) :: betta
  
     alpha = cone; betta = czero
     call zgemm('N', 'N', ndim, ndim, ndim, & 
                         alpha, amat, ndim, &
                                umat, ndim, &
                      betta, tmp_mat, ndim  )
  
     alpha = cone; betta = czero
     call zgemm('C', 'N', ndim, ndim, ndim, &
                         alpha, umat, ndim, &
                             tmp_mat, ndim, &
                         betta, amat, ndim  )
  
  
     return
  end subroutine atomic_tran_repr_cmpl

!>>> atomic_tran_repr_real: transformation from one representation to 
!!>>> another representation, real version
  subroutine atomic_tran_repr_real( ndim, amat, umat )
     use constants, only : dp, zero, one
  
     implicit none
  
! external variables
     integer, intent(in) :: ndim
     real(dp), intent(inout) :: amat(ndim, ndim)
     real(dp), intent(in) :: umat(ndim, ndim)
  
! local variables
     real(dp) :: tmp_mat(ndim, ndim)
     real(dp) :: alpha
     real(dp) :: betta
  
     alpha = one; betta = zero
     call dgemm('N', 'N', ndim, ndim, ndim, &
                         alpha, amat, ndim, &
                                umat, ndim, &
                      betta, tmp_mat, ndim  )
  
     alpha = one; betta = zero 
     call dgemm('T', 'N', ndim, ndim, ndim, &
                         alpha, umat, ndim, &
                             tmp_mat, ndim, &
                         betta, amat, ndim  )
  
     return
  end subroutine atomic_tran_repr_real
 
