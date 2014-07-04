! make single particle related matrices, including:
! crystal field, spin-orbital coupling, Coulomb inteartion U
  subroutine atomic_make_spmat()
     use m_sp_mat

     implicit none

! first, allocate memory for these single particle matrices
     call alloc_m_sp_mat()

! second, make crystal field 
     if (icf .eqv. .true.) then
         call atomic_make_cf()
     else
         sp_cf_mat = czero
     endif

! third, make spin-orbital coupling 
     if (isoc .eqv. .true.) then
         call atomic_make_soc()
     else
         sp_soc_mat = czero
     endif

! fourth, make Coulomb interaction U
     if (icu == 1) then
! Kanamori type
         call atomic_make_cumat_kanamori()
     else
! Slater Integral type
         call atomic_make_cumat_slater()
     endif

     return
  end subroutine atomic_make_spmat 

!>>> make crystal field
  subroutine atomic_make_cf()
     use constants
     use m_sp_mat

     implicit none

! local variables
     logical :: exists

! loop index
     integer :: i, j

! dummy variables
     integer :: int1, int2
     real(dp) :: real1, real2

! we read crystal field from file "atomic.cf.in"
! inquire file
     inquire(file='atomic.cf.in', exist=exists)

     if (exists .eqv. .true.) then
         open(mytmp, file='atomic.cf.in')
         do i=1, norbs
             do j=1, norbs
                 read(mytmp, *) int1, int2, real1, real2
                 sp_cf_mat(j,i) = dcmplx(real1, real2)
             enddo
         enddo 
     else
         call atomic_print_error('atomic_make_cf', 'no file atomic.cf.in !')
     endif 

     return
  end subroutine atomic_make_cf

!>>> make spin-orbital coupling 
  subroutine atomic_make_soc()
     use control
     use m_sp_mat

     implicit none

     if (nband == 3) then
         call atomic_make_soc_3band()
     elseif(nband == 5) then
         call atomic_make_soc_5band()
     elseif(nband == 7) then
         call atomic_make_soc_7band()
     else
         call atomic_print_error('atomic_make_soc', 'not implementd!')
     endif

     return 
  end subroutine atomic_make_soc

!>>> make spin-orbital coupling matrix for 3 band
   subroutine atomic_make_soc_3band()
      use constants
      use control
      use m_sp_mat

      implicit none

! we use the same default orbital order as in WANNIER90 package 
! the t2g orbital order of d(l=2) is 
! |dxz,up>, |dxz,dn>, |dyz,up>, |dyz,dn>, |dxy,up>, |dxy,dn>
! the coressponding p(l=1) orbital order is
! |py,up>,  |py,dn>,  |px,up>,  |px,dn>,  |pz,up>,  |pz,dn>
      sp_soc_mat = czero
      sp_soc_mat(1,3) =  czi;     sp_soc_mat(1,6) =  -czi
      sp_soc_mat(2,4) = -czi;     sp_soc_mat(2,5) =  -czi
      sp_soc_mat(3,1) = -czi;     sp_soc_mat(3,6) =  cone
      sp_soc_mat(4,2) =  czi;     sp_soc_mat(4,5) = -cone
      sp_soc_mat(5,2) =  czi;     sp_soc_mat(5,4) = -cone
      sp_soc_mat(6,1) =  czi;     sp_soc_mat(6,3) =  cone
! please note: minus sign for spin-orbital coupling strength         
      sp_soc_mat = -sp_soc_mat * lambda / two 

   end subroutine atomic_make_soc_3band

!>>> make spin-orbital coupling matrix for 5 band
   subroutine atomic_make_soc_5band()
      use constants
      use control
      use m_sp_mat

      implicit none

! local variables
! sqrt(3)
      real(dp) :: sqrt3

      sqrt3 = sqrt(3.0_dp)

! we use the same default orbital order as in WANNIER90 package 
! the orbital order is: 
! |dz2,up>, |dz2,dn>, |dxz,up>, |dxz,dn>, |dyz,up>, |dyz,dn>, |dx2-y2,up>, |dx2-y2,dn>, |dxy,up>, |dxy,dn>
      sp_soc_mat = czero
      sp_soc_mat(1, 4) = -sqrt3*cone;      sp_soc_mat(1, 6) = sqrt3*czi
      sp_soc_mat(2, 3) =  sqrt3*cone;      sp_soc_mat(2, 5) = sqrt3*czi
      sp_soc_mat(3, 2) =  sqrt3*cone;      sp_soc_mat(3, 5) = -czi
      sp_soc_mat(3, 8) =       -cone;      sp_soc_mat(3, 10)= czi
      sp_soc_mat(4, 1) = -sqrt3*cone;      sp_soc_mat(4, 6) = czi
      sp_soc_mat(4, 7) =        cone;      sp_soc_mat(4, 9) = czi
      sp_soc_mat(5, 2) =  -sqrt3*czi;      sp_soc_mat(5, 3) = czi
      sp_soc_mat(5, 8) =        -czi;      sp_soc_mat(5, 10)=-cone
      sp_soc_mat(6, 1) =  -sqrt3*czi;      sp_soc_mat(6, 4) = -czi
      sp_soc_mat(6, 7) =        -czi;      sp_soc_mat(6, 9) = cone
      sp_soc_mat(7, 4) =        cone;      sp_soc_mat(7, 6) = czi
      sp_soc_mat(7, 9) =     -two*czi
      sp_soc_mat(8, 3) =       -cone;      sp_soc_mat(8, 5) = czi
      sp_soc_mat(8,10) =      two*czi
      sp_soc_mat(9, 4) =        -czi;      sp_soc_mat(9, 6) = cone
      sp_soc_mat(9, 7) =     two*czi;    
      sp_soc_mat(10,3) =        -czi;      sp_soc_mat(10,5) = -cone
      sp_soc_mat(10,8) =    -two*czi; 
  
! scale the SOC strength lambda
      sp_soc_mat = sp_soc_mat * lambda / two

      return 
   end subroutine atomic_make_soc_5band

!>>> make spin-orbital coupling matrix for 7 band
   subroutine atomic_make_soc_7band()
      use constants
      use control
      use m_sp_mat

      implicit none
      
      write(mystd,*) "not implemented now!"

      return
   end subroutine atomic_make_soc_7band

!>>> make Coulomb interaction U, this subroutine is taken from 
!>>> duliang's atomic software
!>>>  $U_{\alpha\betta\delta\gamma}$ coefficent matrix for Uij  <<<!
!--------------------three band for example------------------------!
!> $f_{\alpha}^{\dagger}f_{\betta}^{\dagger}f_{\delta}f_{\gamma}$ <!
! norbs   bandindex(?band)    spinindex(?spin)    representation   !
!   1          1               1->\uparrow             1\up        !
!   2          1               0->\doarrow             1\do        !
!   3          2               1->\uparrow             2\up        !
!   4          2               0->\doarrow             2\do        !
!   5          3               1->\uparrow             3\up        !
!   6          3               0->\doarrow             3\do        !
!------------------------------------------------------------------!
  subroutine atomic_make_cumat_kanamori()
     use constants
     use control
     use m_sp_mat

     implicit none

! local varibales
! loop index over orbits
     integer :: alpha, betta
     integer :: delta, gamma

! band index and spin index
! band index of alpha and betta
     integer :: aband, bband 

! band index of delta and gamma
     integer :: dband, gband 

! spin index of alpha and betta
     integer :: aspin, bspin 

! spin index of delta and gamma
     integer :: dspin, gspin 

! auxiliary variables
     real(dp) :: dtmp

! initialize cumat to zero
     sp_cu_mat = czero

! loop for creation operator
     alphaloop: do alpha=1,norbs-1
     bettaloop: do betta=alpha+1,norbs

! loop for destroy operator
        gammaloop: do gamma=1,norbs-1
        deltaloop: do delta=gamma+1,norbs
            aband = (alpha+1)/2; aspin = mod(alpha,2)
            bband = (betta+1)/2; bspin = mod(betta,2)
            gband = (gamma+1)/2; gspin = mod(gamma,2)
            dband = (delta+1)/2; dspin = mod(delta,2)

! here we use "res" due to overlap between "Uv and Jz"
            dtmp = zero

! intraorbital Coulomb interaction
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if ((aband.eq.bband) .and. (aspin.ne.bspin)) then
                    dtmp = dtmp + Uc
                endif
            endif

! interorbital Coulomb interaction
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if (aband .ne. bband) then
                    dtmp = dtmp + Uv
                endif
            endif

! Hund's exchange interaction 
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if ((aband.ne.bband) .and. (aspin.eq.bspin)) then
                    dtmp = dtmp - Jz
                endif
            endif
           
! spin flip term
            if ((aband.eq.gband) .and. (bband.eq.dband)) then
                if ((aspin.ne.gspin) .and. (bspin.ne.dspin) .and. (aspin.ne.bspin)) then
                    dtmp = dtmp - Js
                endif
            endif
          
! pair hopping term
            if ((aband.eq.bband) .and. (dband.eq.gband) .and. (aband.ne.dband)) then
                if ((aspin.ne.bspin) .and. (dspin.ne.gspin) .and. (aspin.eq.gspin)) then
                    dtmp = dtmp + Jp
                endif
            endif
                 
            sp_cu_mat(alpha, betta, delta, gamma) = dtmp

        enddo deltaloop ! over delta={gamma+1,norbs} loop
        enddo gammaloop ! over gamma={1,norbs-1} loop
     enddo bettaloop ! over betta={alpha+1,norbs} loop
     enddo alphaloop ! over alpha={1,norbs-1} loop

     return
  end subroutine atomic_make_cumat_kanamori
