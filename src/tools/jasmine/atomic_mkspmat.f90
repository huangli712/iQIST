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
