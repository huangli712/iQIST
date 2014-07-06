!>>> make natural basis
  subroutine atomic_make_natural()
     use constants
     use control
     use m_sp_mat

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

! make umat from origional basis to natural basis: sp_tran_umat
! and set the eimp: sp_eimp_mat 
     if ( itask==1) then
         if ( isoc==0 .and. (icf==0 .or. icf==1) ) then
! for model calculation, no spin-orbital coupling, no crystal field or
! crystal field is diagonal, the real orbital basis is the natural basis
! so, we do nothing here, just return
             write(mystd,"(4X, a)") "atomic >>> The natural basis is real orbital basis !"
             return
         elseif ( isoc==0 .and. icf==2 ) then
             call atomic_2natural_case1() 
         elseif ( isoc==1 .and. icf==0 ) then
             call atomic_2natural_case2()
         elseif ( isoc==1 .and. icf>0 ) then
             call atomic_2natural_case3()
         endif
     endif

! we need transform Coulomb interaction U
! for Slater-Cordon parameters Coulomb interaction U
! we first need to transfrom sp_cu_mat from complex orbital basis to real orbital basis
     if ( icu == 2 ) then
         call atomic_make_umat_r2c( umat_r2c )
         umat_r2c = transpose( umat_c2r )
         call atomic_tran_cumat( umat_r2c, sp_cu_mat, tmp_mat )
         sp_cu_mat = tmp_mat
     endif

     call atomic_tran_cumat(sp_tran_umat, sp_cu_mat, tmp_mat) 
     sp_cu_mat = tmp_mat

     return
  end subroutine atomic_make_natural

!>>> make natural basis for non-diagonal crystal field without spin-orbital coupling
  subroutine atomic_2natural_case1()
     use constants
     use control
     use m_sp_mat

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

! set sp_eimp_mat
     sp_eimp_mat = sp_cf_mat   

! get eimp for nospin freedom
     call atomic_mat_2nospin(sp_eimp_mat, eimp_nospin)

! diagonalize eimp_nospin to get natural basis
     eimp_nospin_real = real(eimp_nospin)
     call dmat_dsyev(nband, nband, eimp_nospin_real, eigval, eigvec)
 
     eimp_nospin = czero
     do i=1, nband
         eimp_nospin(i,i) = eigval(i) 
     enddo
     umat_nospin = eigvec
 
     call atomic_mat_2spin(eimp_nospin, sp_eimp_mat) 
     call atomic_mat_2spin(umat_nospin, sp_tran_mat)

     return
  end subroutine atomic_2natural_case1

!>>> make natural basis for the case without crystal field and with spin-orbital coupling
!    for this special case, the natural basis is |J^2,Jz>
  subroutine atomic_2natural_case2()
     use constants
     use control
     use m_sp_mat

     implicit none

! local variables
! umat from real orbital basis to complex orbital basis
     complex(dp) :: umat_r2c(norbs, norbs)

! umat from complex orbital basis to |j^2, jz> basis
     complex(dp) :: umat_c2j(norbs, norbs)

! umat from real orbital basis to |j^2, jz> basis
     complex(dp) :: umat_r2j(norbs, norbs)

! set sp_eimp_mat
     sp_eimp_mat = sp_soc_mat   

! get umat_r2c
     call atomic_make_umat_r2c(umat_r2c)

! get umat_c2j
     call atomic_make_umat_c2j(umat_c2j)

     call zmat_zgemm0(norbs, umat_r2c, umat_c2j, umat_r2j)
 
     sp_tran_umat = umat_r2j
 
! transform sp_eimp_mat to natural basis
     call atomic_tran_represent(norbs, sp_eimp_mat, sp_tran_umat)   

     return
  end subroutine atomic_2natural_case2

!>>> make natural basis for the case with crystal field and with spin-orbital coupling
  subroutine atomic_2natural_case3()
     use constants
     use control
     use m_sp_mat

     implicit none

! local variables
! umat from real orbital basis to complex orbital basis
     complex(dp) :: umat_r2c(norbs, norbs)

! umat from complex orbital basis to natural basis
     complex(dp) :: umat_c2n(norbs, norbs)

! real version of eimp_mat
     real(dp) :: eimp_mat(norbs, norbs)

! eigen value 
     real(dp) :: eigval(norbs)

! eigen vector
     real(dp) :: eigvec(norbs, norbs)

! get umat_r2c
     call atomic_make_umat_r2c(umat_r2c)

! transfrom sp_cf_mat and sp_soc_mat to complex orbital basis
     call atomic_tran_represent(norbs, sp_cf_mat, umat_r2c)
     call atomic_tran_represent(norbs, sp_soc_mat, umat_r2c)

! set sp_eimp_mat
     sp_eimp_mat = sp_soc_mat + sp_cf_mat   

     eimp_mat = real(sp_eimp_mat)

     call dmat_dsyev(norbs, norbs, eimp_mat, eigval, eigvec)

     umat_c2n = eigvec
     call zmat_zgemm0(norbs, umat_r2c, umat_c2n, sp_tran_umat)
 
! transform sp_eimp_mat to natural basis
     call atomic_tran_represent(norbs, sp_eimp_mat, umat_c2n)   
     
     return
  end subroutine atomic_2natural_case3
