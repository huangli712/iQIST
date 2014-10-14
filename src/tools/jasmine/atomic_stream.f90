!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_config
!!! source  : atomic_config.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : set control parameters 
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_config: read config parameters from file 'atom.config.in'
  subroutine atomic_config()
     use constants, only : dp, mytmp
     use control

     use parser, only : p_create, p_destroy, p_parse, p_get
  
     implicit none
  
! local variables
! file status
     logical :: exists
     
!----------------------------------------------------------------
     itask  = 1           ! type of task
     ictqmc = 1           ! type of CTQMC algorithm
     icf    = 0           ! type of crystal field
     isoc   = 0           ! type of spin-orbital coupling (SOC)
     icu    = 1           ! type of Coulomb interaction
     
!---------------------------------------------------------------- 
     nband = 1            ! number of bands
     nspin = 2            ! number of spins
     norbs = nband*nspin  ! number of orbits
     ncfgs = 2**norbs     ! number of many-body configurations

!----------------------------------------------------------------
     Uc = 2.00_dp         ! intraorbital Coulomb interaction
     Uv = 2.00_dp         ! interorbital Coulomb interaction
     Jz = 0.00_dp         ! Hund's exchange interaction
     Js = 0.00_dp         ! spin-flip interaction
     Jp = 0.00_dp         ! pair-hopping interaction
  
!----------------------------------------------------------------
     Ud = 2.00_dp         ! Ud
     JH = 0.00_dp         ! JH
     F0 = 0.00_dp         ! F0
     F2 = 0.00_dp         ! F2
     F4 = 0.00_dp         ! F4
     F6 = 0.00_dp         ! F6
  
!----------------------------------------------------------------
     lambda = 0.00_dp     ! spin-orbit coupling parameter
     mune   = 0.00_dp     ! chemical potential
  
!----------------------------------------------------------------
! file status
     exists = .false.
  
! inquire the input file status
     inquire( file="atom.config.in", exist=exists )
  
! read parameters from atom.config.in
     if ( exists .eqv. .true. ) then
!----------------------------------------------------------------
         call p_create()
         call p_parse('atom.config.in')
!----------------------------------------------------------------
         call p_get('nband', nband)
!----------------------------------------------------------------
         call p_get('itask',  itask)
         call p_get('ictqmc', ictqmc)
         call p_get('icf',    icf)
         call p_get('isoc',   isoc)
         call p_get('icu',    icu)
!----------------------------------------------------------------
         call p_get('Uc',     Uc) 
         call p_get('Uv',     Uv) 
         call p_get('Jz',     Jz) 
         call p_get('Js',     Js) 
         call p_get('Jp',     Jp) 
!----------------------------------------------------------------
         call p_get('Ud',     Ud) 
         call p_get('JH',     JH) 
!----------------------------------------------------------------
         call p_get('lambda', lambda)
         call p_get('mune',   mune)
!----------------------------------------------------------------
         call p_destroy()
!----------------------------------------------------------------
! calculate the norbs and ncfgs
         norbs = nband * nspin
         ncfgs = 2 ** norbs 

! calculate F0, F2, F4, F6 here
         F0 = Ud
         if (nband == 5) then
             F2 = JH * 14.0_dp / 1.625_dp 
             F4 = 0.625_dp * F2
         elseif(nband == 7) then
             F2 = JH * 6435.0_dp / (286.0_dp + (195.0_dp * 451.0_dp / 675.0_dp) &
                                            + (250.0_dp * 1001.0_dp / 2025.0_dp))
             F4 = 451.0_dp / 675.0_dp * F2
             F6 = 1001.0_dp / 2025.0_dp * F2
         endif
!----------------------------------------------------------------
     else
         call s_print_error('atomic_config', 'no file atom.config.in !')
     endif
  
     return
  end subroutine atomic_config
!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_read_cf
!!!           atomic_read_eimp
!!!           atomic_read_umat
!!! source  : atomic_natural.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : read data from files
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_read_cf: read crystal field from file 'atomic.cf.in'
  subroutine atomic_read_cf()
     use constants, only : mytmp, dp, zero
     use m_spmat, only : cfmat
  
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

!!>>> atomic_read_eimp: read on-site impurity energy 
!!>>> from file 'atomic.eimp.in'
  subroutine atomic_read_eimp()
     use constants, only : mytmp, dp, zero
     use control, only : norbs

     use m_spmat, only : eimpmat
  
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
  
!!>>> atomic_read_umat: read the transformation matrix 
!!>>> tran_umat from file 'atomic.umat.in'
  subroutine atomic_read_umat()
     use constants, only : mytmp, dp, zero
     use control, only : norbs

     use m_spmat, only : tran_umat
  
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
!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_check_config
!!!           atomic_check_realmat
!!! source  : atomic_check.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : do some checks
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_check_config: check the validity of input config parameters
  subroutine atomic_check_config()
     use constants, only : mystd, zero
     use control
  
! local variables
     logical :: lpass
  
     lpass = .true.
  
! check nband
     if (nband <= 0) then
         write(mystd, '(2X,a)') 'jasmine >>> ERROR: number of bands must &
                                                    be larger than zero !'
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
         write(mystd, '(2X,a)') 'jasmine >>> ERROR: only support spin-orbital &
                                                  coupling for nband=3, 5, 7 !'
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
         write(mystd, '(2X,a)') 'jasmine >>> ERROR: only support Slater-Cordon type &
                                                Coulomb interaction for nband=5, 7 !'
         write(mystd, *)
         lpass = .false.
     endif
  
  
! check Uc, Uv, Jz, Js, Jp, Ud, JH
     if (Uc < zero .or. Uv < zero .or. Jz < zero .or. Js < zero .or. &
                                  Jp < zero .or. Ud < zero .or. JH < zero) then
         write(mystd, '(2X,a)') 'jasmine >>> ERROR: Uc, Uv, Jz, Js, Jp, Ud, JH &
                                                        must be larger than zero !'
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
         write(mystd, '(2X,a)') 'jasmine >>> ERROR: CTQMC good quantum numbers (N,Sz) &
                           algorithm is NOT supported for spin-orbital coupling case ! '
         write(mystd, *)
         lpass = .false.
     endif
     if (ictqmc == 4 .and. isoc == 1) then
         write(mystd, '(2X,a)') 'jasmine >>> ERROR: CTQMC good quantum numbers (N,Sz,Ps) &
                              algorithm is NOT supported for spin-orbital coupling case ! '
         write(mystd, *)
         lpass = .false.
     endif
     if (ictqmc == 4 .and. icu == 2) then
         write(mystd, '(2X,a)') 'jasmine >>> ERROR: CTQMC good quantum numbers (N,Sz,Ps) &
                   algorithm is NOT supported for Slater-Cordon type Coulomb interaction U'
         write(mystd, *)
         lpass = .false.
     endif
     if (ictqmc == 5 .and. isoc == 0) then
         write(mystd, '(2X,a)') 'jasmine >>> ERROR: CTQMC good quantum numbers (N,Jz) &
                             algorithm is ONLY supported for spin-orbital coupling case !'
         write(mystd, *)
         lpass = .false.
     endif
     if (ictqmc == 5 .and. isoc == 1 .and. icf /= 0) then
         write(mystd, '(2X,a)') 'jasmine >>> ERROR: CTQMC good quantum numbers (N,Jz) &
          algorithm is NOT supported for spin-orbital coupling plus crystal field case !'
         write(mystd, *)
         lpass = .false.
     endif
     
     if (lpass == .true.) then
         write(mystd, '(2X,a)') 'jasmine >>> Good News: all control parameters are OK !'
     else
         call s_print_error('atomic_check_config', 'Found some wrong setting of parameters, &
                                                   please check the "atom.config.in" file !')
     endif
  
     return
  end subroutine atomic_check_config

!!>>> atomic_check_realmat: check whether a matrix is real
  subroutine atomic_check_realmat(ndim, mat, lreal)
     use constants, only : dp, eps6
  
     implicit none
  
! external variables
! dimension of the matrix
     integer, intent(in) :: ndim

! the matrix to be checked
     complex(dp), intent(in) :: mat(ndim, ndim)

! whether Hamiltonian is real
     logical, intent(out) :: lreal
  
! local variables
     integer :: i, j
  
     do i=1, ndim
         do j=1, ndim
             if ( aimag(mat(j,i)) > eps6 ) then
                 lreal = .false. 
                 return
             endif
         enddo
     enddo
  
     lreal = .true.
  
     return
  end subroutine atomic_check_realmat
!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_make_spmat
!!!           atomic_make_soc
!!!           atomic_mksoc_3band
!!!           atomic_mksoc_5band
!!!           atomic_mksoc_7band
!!!           atomic_mkcumat_kanamori
!!!           atomic_mkcumat_slater
!!! source  : atomic_mkspmat.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : make single particle related matrices, including crystal field,
!!!           spin-orbital coupling, and Coulomb interaction U tensor
!!! status  : unstable
!!! comment : subroutine atomic_mkcumat_kanamori and atomic_mkcumat_slater are 
!!!           modified from Dr. LiangDu's (duleung@gmail.com) atomic program
!!!-------------------------------------------------------------------------

!!>>> atomic_make_spmat: make single particle related matrices, including
!!>>> crystal field (CF), spin-orbit coupling (SOC), and Coulomb interaction
!!>>> U tensor
  subroutine atomic_make_spmat()
     use constants, only : czero
     use control, only : itask, icf, isoc, icu

     use m_spmat, only : cfmat, socmat
  
     implicit none
  
! make crystal field and spin-orbital coupling 
! make natural basis inside
     if (itask == 1) then 
! crysal field
         if (icf > 0) then
! we read the non-zero elements of 
! crystal field from a file "atom.cf.in".
! the crystal field is defined on real orbital basis
! at present, we only support real crystal field, 
! so, the elements in this file provided by user must be real
             call atomic_read_cf()
         else
             cfmat = czero
         endif
! spin-orbit coupling
         if (isoc > 0) then
! make an atomic on-site SOC, $\lambda * L * S$
! it is defined on the complex orbital basis
             call atomic_make_soc()
         else
             socmat = czero
         endif
! make natural basis outside
     else 
! read the eimp (CF+SOC) matrices on natural basis
! this matrix should be a diagonal matrix, and the elements must be real
         call atomic_read_eimp()
! read the transformation matrices used to transfer eimp 
! from original basis to natural basis 
! without SOC, the original basis is the real orbital basis
! with SOC, the original basis is the complex orbital basis
! at present, we just only support real numbers of this umat
         call atomic_read_umat()
     endif
  
! make Coulomb interaction U
     if (icu == 1) then
! Kanamori parameters type
! it is defined on real orbital basis 
         call atomic_mkcumat_kanamori()
     else
! Slater-Cordon parameters type
! it is defined on complex orbital basis
         call atomic_mkcumat_slater()
     endif
  
     return
  end subroutine atomic_make_spmat 

!!>>> atomic_make_soc: make spin-orbital coupling (SOC)
  subroutine atomic_make_soc()
     use constants, only : two
     use control, only : nband, lambda
     use m_spmat, only : socmat
  
     implicit none
  
     if (nband == 3) then
         call atomic_mksoc_3band(socmat)
! for 3 bands system, there is a minus sign
         socmat = -socmat * lambda / two
     elseif(nband == 5) then
         call atomic_mksoc_5band(socmat)
         socmat = socmat * lambda / two
     elseif(nband == 7) then
         call atomic_mksoc_7band(socmat)
         socmat = socmat * lambda / two
     else
         call s_print_error('atomic_make_soc', 'not implementd!')
     endif
  
     return 
  end subroutine atomic_make_soc

!>>> atomic_mksoc_3band: make spin-orbit coupling matrix for 3 bands
  subroutine atomic_mksoc_3band(socmat)
     use constants, only : dp

     implicit none
  
! external variables
     complex(dp), intent(out) :: socmat(6,6)
  
! local variables
     real(dp) :: sqrt2
  
     sqrt2 = sqrt(2.0_dp)
     
! make SOC on complex orbital basis, the orbital order is:
! |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>      
     socmat = dcmplx(0.0_dp, 0.0_dp)
  
     socmat(1,1) = -1.0_dp
     socmat(4,1) = sqrt2 
     socmat(2,2) =  1.0_dp
     socmat(6,3) = sqrt2
     socmat(1,4) = sqrt2
     socmat(5,5) = 1.0_dp
     socmat(3,6) = sqrt2
     socmat(6,6) = -1.0_dp
  
     return
  end subroutine atomic_mksoc_3band

!!>>> atomic_mksoc_5band: make spin-orbit coupling matrix for 5 bands
  subroutine atomic_mksoc_5band(socmat)
     use constants, only : dp

     implicit none
  
! external variables
     complex(dp), intent(out) :: socmat(10,10)
  
! local variables
     real(dp) :: sqrt6
  
     sqrt6 = sqrt(6.0_dp)
! make SOC on complex orbital basis, the orbital order is:
! |-2,up>, |-2,dn>, |-1,up>, |-1,dn>, |0,up>, |0,dn>, |1,up>, |1,dn>, |2,up>, |2,dn>      
  
     socmat = dcmplx(0.0_dp, 0.0_dp)
  
     socmat(1,1) = -2.0_dp 
     socmat(4,1) =  2.0_dp
     socmat(2,2) =  2.0_dp
     socmat(3,3) = -1.0_dp
     socmat(6,3) =  sqrt6
     socmat(1,4) =  2.0_dp
     socmat(4,4) =  1.0_dp
     socmat(8,5) =  sqrt6
     socmat(3,6) =  sqrt6
     socmat(7,7) =  1.0_dp
     socmat(10,7)=  2.0_dp
     socmat(5,8) =  sqrt6  
     socmat(8,8) = -1.0_dp
     socmat(9,9) =  2.0_dp
     socmat(7,10)=  2.0_dp 
     socmat(10,10)= -2.0_dp
  
     return 
  end subroutine atomic_mksoc_5band

!!>>> atomic_mksoc_7band: make spin-orbit coupling matrix for 7 bands
  subroutine atomic_mksoc_7band(socmat)
     use constants, only : dp

     implicit none
  
! local variables
     real(dp) :: sqrt6
     real(dp) :: sqrt10
     real(dp) :: sqrt12
  
! external variables
     complex(dp), intent(out) :: socmat(14,14)    
  
     sqrt6  = sqrt(6.0_dp)
     sqrt10 = sqrt(10.0_dp)
     sqrt12 = sqrt(12.0_dp)
  
     socmat = dcmplx(0.0_dp, 0.0_dp)
  
     socmat(1,1) = -3.0_dp
     socmat(4,1) = sqrt6
     socmat(2,2) = 3.0_dp
     socmat(3,3) = -2.0_dp
     socmat(6,3) = sqrt10 
     socmat(1,4) = sqrt6
     socmat(4,4) = 2.0_dp
     socmat(5,5) = -1.0_dp
     socmat(8,5) = sqrt12
     socmat(3,6) = sqrt10 
     socmat(6,6) = 1.0_dp
     socmat(10,7) = sqrt12
     socmat(5,8) = sqrt12
     socmat(9,9) = 1.0_dp
     socmat(12,9) = sqrt10
     socmat(7,10) = sqrt12
     socmat(10,10) = -1.0_dp
     socmat(11,11) =  2.0_dp
     socmat(14,11) =  sqrt6
     socmat(9,12) = sqrt10
     socmat(12,12) =  -2.0_dp
     socmat(13,13) =  3.0_dp
     socmat(11,14) =  sqrt6
     socmat(14,14) =  -3.0_dp
  
     return
  end subroutine atomic_mksoc_7band

!!>>> atomic_mkcumat_kanamori: make Coulomb interaction U according to
!!>>> Kanamori parameterized Hamiltonian
  subroutine atomic_mkcumat_kanamori()
     use constants, only : dp, czero, zero
     use control, only : norbs, Uc, Uv, Jz, Js, Jp

     use m_spmat, only : cumat
  
     implicit none
  
! local varibales
! orbital index
     integer :: alpha, betta
     integer :: delta, gamma

! band index and spin index
     integer :: aband, bband 
     integer :: dband, gband 
     integer :: aspin, bspin 
     integer :: dspin, gspin 

! dummy variables
     real(dp) :: dtmp
  
! initialize cumat to zero
     cumat = czero
  
! loop for creation operators
     alphaloop: do alpha=1,norbs-1
     bettaloop: do betta=alpha+1,norbs
  
! loop for annihilation operators
        gammaloop: do gamma=1,norbs-1
        deltaloop: do delta=gamma+1,norbs
            aband = (alpha+1)/2; aspin = mod(alpha,2)
            bband = (betta+1)/2; bspin = mod(betta,2)
            gband = (gamma+1)/2; gspin = mod(gamma,2)
            dband = (delta+1)/2; dspin = mod(delta,2)
  
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
                 
            cumat(alpha, betta, delta, gamma) = dtmp
  
        enddo deltaloop ! over delta={gamma+1,norbs} loop
        enddo gammaloop ! over gamma={1,norbs-1} loop
     enddo bettaloop ! over betta={alpha+1,norbs} loop
     enddo alphaloop ! over alpha={1,norbs-1} loop
  
     return
  end subroutine atomic_mkcumat_kanamori

!!>>> atomic_mkcumat_slater: make Coulomb interation U, according to 
!!>>> Slater-Cordon parameterized Hamiltonian
  subroutine atomic_mkcumat_slater()
     use constants, only : dp, zero, half
     use control, only : nband, norbs, F0, F2, F4, F6
     use m_spmat, only : cumat
  
     implicit none
  
! local variables
! Slater-Cordon parameters
     real(dp), allocatable :: slater_cordon(:)

! gaunt coefficients
     real(dp), allocatable :: gaunt(:,:,:)
  
! orbital momentum quantum number
     integer :: l

! loop index
     integer :: i
     integer :: alpha, betta
     integer :: delta, gamma
     integer :: aband, aspin
     integer :: bband, bspin
     integer :: dband, dspin
     integer :: gband, gspin

! dummy variables
     real(dp) :: res
  
  
! allocate memory for slater_cordon and gaunt and then build them
     if (nband == 5) then
         l = 2
         allocate(slater_cordon(0:2*l))     
         slater_cordon = zero
         slater_cordon(0) = F0
         slater_cordon(2) = F2
         slater_cordon(4) = F4
         allocate(gaunt(-l:l, -l:l, 0:2*l))
         call atomic_gaunt_5band(gaunt) 
     elseif(nband == 7) then
         l = 3
         allocate(slater_cordon(0:2*l))     
         slater_cordon = zero
         slater_cordon(0) = F0
         slater_cordon(2) = F2
         slater_cordon(4) = F4
         slater_cordon(6) = F6
         allocate(gaunt(-l:l, -l:l, 0:2*l))
         call atomic_gaunt_7band(gaunt)
     else
         call s_print_error('atomic_make_cumat_slater', 'not implemented for this nband!')
     endif
  
! make Coulomb interaction U matrix
     do alpha=1,norbs
     do betta=1,norbs
         aband = (alpha-1)/2-l
         bband = (betta-1)/2-l
         aspin = mod(alpha, 2)
         bspin = mod(betta, 2)
  
         do delta=1,norbs
         do gamma=1,norbs
             dband = (delta-1)/2-l
             gband = (gamma-1)/2-l
             dspin = mod(delta, 2)
             gspin = mod(gamma, 2)
  
             if ((alpha .eq. betta) .or. (delta .eq. gamma)) cycle
  
             if ((aband + bband) .ne. (dband + gband)) cycle
             if ((aspin .ne. gspin) .or. (bspin .ne. dspin)) cycle
  
             res = zero
             do i=0, 2*l, 2
                 res = res + gaunt(aband, gband, i) * gaunt(dband, bband, i) * slater_cordon(i)
             enddo
             cumat(alpha, betta, delta, gamma) = res
         enddo ! over gamma={1,norbs} loop
         enddo ! over delta={1,norbs} loop
  
     enddo ! over betta={1,norbs} loop
     enddo ! over alpha={1,norbs} loop
  
     cumat = half * cumat
  
! deallocate memory
     if (allocated(slater_cordon)) deallocate(slater_cordon) 
     if (allocated(gaunt))         deallocate(gaunt)
  
     return
  end subroutine atomic_mkcumat_slater
!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_gaunt_5band
!!!           atomic_gaunt_7band
!!! source  : atomic_gaunt.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : make gaunt coefficients
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_gaunt_5band: build gaunt coefficients for 5 band case
  subroutine atomic_gaunt_5band(gaunt)
     use constants, only : dp, zero, one
     
! external variables
     real(dp), intent(out) :: gaunt(-2:2, -2:2, 0:4)
  
     gaunt = zero
  
     gaunt(-2, -2, 0) = one
     gaunt(-1, -1, 0) = one
     gaunt(0,   0, 0) = one
     gaunt(1,   1, 0) = one
     gaunt(2,   2, 0) = one
  
     gaunt(-2, -2, 2) = -sqrt(4.0_dp/49.0_dp) 
     gaunt(-2, -1, 2) =  sqrt(6.0_dp/49.0_dp);   gaunt(-1, -2, 2) = gaunt(-2, -1, 2) * (-1)**(-2+1) 
     gaunt(-2,  0, 2) = -sqrt(4.0_dp/49.0_dp);   gaunt(0,  -2, 2) = gaunt(-2,  0, 2) * (-1)**(-2-0)
     gaunt(-1, -1, 2) =  sqrt(1.0_dp/49.0_dp)
     gaunt(-1,  0, 2) =  sqrt(1.0_dp/49.0_dp);   gaunt(0,  -1, 2) = gaunt(-1,  0, 2) * (-1)**(-1-0)
     gaunt(-1,  1, 2) = -sqrt(6.0_dp/49.0_dp);   gaunt(1,  -1, 2) = gaunt(-1,  1, 2) * (-1)**(-1-1)
     gaunt(0,   0, 2) =  sqrt(4.0_dp/49.0_dp)
     gaunt(1,  -1, 2) = -sqrt(6.0_dp/49.0_dp);   gaunt(-1,  1, 2) = gaunt(1,  -1, 2) * (-1)**(1+1)
     gaunt(1,   0, 2) =  sqrt(1.0_dp/49.0_dp);   gaunt(0,   1, 2) = gaunt(1,   0, 2) * (-1)**(1-0)
     gaunt(1,   1, 2) =  sqrt(1.0_dp/49.0_dp)
     gaunt(2,   0, 2) = -sqrt(4.0_dp/49.0_dp);   gaunt(0,   2, 2) = gaunt(2,   0, 2) * (-1)**(2-0)
     gaunt(2,   1, 2) =  sqrt(6.0_dp/49.0_dp);   gaunt(1,   2, 2) = gaunt(2,   1, 2) * (-1)**(2-1)
     gaunt(2,   2, 2) = -sqrt(4.0_dp/49.0_dp)
     
     gaunt(-2, -2, 4) =  sqrt( 1.0_dp/441.0_dp)
     gaunt(-2, -1, 4) = -sqrt( 5.0_dp/441.0_dp); gaunt(-1, -2, 4) = gaunt(-2, -1, 4) * (-1)**(-2+1)
     gaunt(-2,  0, 4) =  sqrt(15.0_dp/441.0_dp); gaunt(0,  -2, 4) = gaunt(-2,  0, 4) * (-1)**(-2-0)
     gaunt(-2,  1, 4) = -sqrt(35.0_dp/441.0_dp); gaunt(1,  -2, 4) = gaunt(-2,  1, 4) * (-1)**(-2-1)
     gaunt(-2,  2, 4) =  sqrt(70.0_dp/441.0_dp); gaunt(2,  -2, 4) = gaunt(-2,  2, 4) * (-1)**(-2-2)
     gaunt(-1, -1, 4) = -sqrt(16.0_dp/441.0_dp)
     gaunt(-1,  0, 4) =  sqrt(30.0_dp/441.0_dp); gaunt(0,  -1, 4) = gaunt(-1,  0, 4) * (-1)**(-1-0)
     gaunt(-1,  1, 4) = -sqrt(40.0_dp/441.0_dp); gaunt(1,  -1, 4) = gaunt(-1,  1, 4) * (-1)**(-1-1)
     gaunt( 0,  0, 4) =  sqrt(36.0_dp/441.0_dp)
     gaunt( 1,  0, 4) =  sqrt(30.0_dp/441.0_dp); gaunt(0,   1, 4) = gaunt(1,   0, 4) * (-1)**(1-0)
     gaunt( 1,  1, 4) = -sqrt(16.0_dp/441.0_dp)
     gaunt( 2, -1, 4) = -sqrt(35.0_dp/441.0_dp); gaunt(-1,  2, 4) = gaunt(2,  -1, 4) * (-1)**(2+1)
     gaunt( 2,  0, 4) =  sqrt(15.0_dp/441.0_dp); gaunt( 0,  2, 4) = gaunt(2,   0, 4) * (-1)**(2-0)
     gaunt( 2,  1, 4) = -sqrt( 5.0_dp/441.0_dp); gaunt( 1,  2, 4) = gaunt(2,   1, 4) * (-1)**(2-1)
     gaunt( 2,  2, 4) =  sqrt( 1.0_dp/441.0_dp)
     
     return
  end subroutine atomic_gaunt_5band

!!>>> atomic_gaunt_7band: build gaunt coefficients for 7 band case
  subroutine atomic_gaunt_7band(gaunt)
     use constants, only: dp, zero, one
   
     implicit none
  
! external variables
     real(dp), intent(out) :: gaunt(-3:3, -3:3, 0:6)

     gaunt = zero
  
     gaunt(-3, -3, 0) = one
     gaunt(-2, -2, 0) = one
     gaunt(-1, -1, 0) = one
     gaunt( 0,  0, 0) = one
     gaunt( 1,  1, 0) = one
     gaunt( 2,  2, 0) = one
     gaunt( 3,  3, 0) = one
  
     gaunt(-3, -3, 2) = -sqrt(25.0_dp/225.0_dp)
     gaunt(-3, -2, 2) =  sqrt(25.0_dp/225.0_dp);  gaunt(-2, -3, 2) = gaunt(-3, -2, 2) * (-1.0)**(-3+2)
     gaunt(-3, -1, 2) = -sqrt(10.0_dp/225.0_dp);  gaunt(-1, -3, 2) = gaunt(-3, -1, 2) * (-1.0)**(-3+1)
     gaunt(-2, -1, 2) =  sqrt(15.0_dp/225.0_dp);  gaunt(-1, -2, 2) = gaunt(-2, -1, 2) * (-1.0)**(-2+1)
     gaunt(-2,  0, 2) = -sqrt(20.0_dp/225.0_dp);  gaunt( 0, -2, 2) = gaunt(-2,  0, 2) * (-1.0)**(-2-0)
     gaunt(-1, -1, 2) =  sqrt( 9.0_dp/225.0_dp)
     gaunt(-1,  0, 2) =  sqrt( 2.0_dp/225.0_dp);  gaunt( 0, -1, 2) = gaunt(-1,  0, 2) * (-1.0)**(-1-0)
     gaunt(-1,  1, 2) = -sqrt(24.0_dp/225.0_dp);  gaunt( 1, -1, 2) = gaunt(-1,  1, 2) * (-1.0)**(-1-1)
     gaunt( 0,  0, 2) =  sqrt(16.0_dp/225.0_dp)
     gaunt( 1,  0, 2) =  sqrt( 2.0_dp/225.0_dp);  gaunt( 0,  1, 2) = gaunt( 1,  0, 2) * (-1.0)**( 1-0)
     gaunt( 1,  1, 2) =  sqrt( 9.0_dp/225.0_dp)
     gaunt( 2,  0, 2) = -sqrt(20.0_dp/225.0_dp);  gaunt( 0,  2, 2) = gaunt( 2,  0, 2) * (-1.0)**( 2-0)
     gaunt( 2,  1, 2) =  sqrt(15.0_dp/225.0_dp);  gaunt( 1,  2, 2) = gaunt( 2,  1, 2) * (-1.0)**( 2-1)
     gaunt( 3,  1, 2) = -sqrt(10.0_dp/225.0_dp);  gaunt( 1,  3, 2) = gaunt( 3,  1, 2) * (-1.0)**( 3-1)
     gaunt( 3,  2, 2) =  sqrt(25.0_dp/225.0_dp);  gaunt( 2,  3, 2) = gaunt( 3,  2, 2) * (-1.0)**( 3-2)
     gaunt( 3,  3, 2) = -sqrt(25.0_dp/225.0_dp)
  
     gaunt(-3, -3, 4) =  sqrt( 9.0_dp/1089.0_dp)
     gaunt(-3, -2, 4) = -sqrt(30.0_dp/1089.0_dp); gaunt(-2, -3, 4) = gaunt(-3, -2, 4) * (-1.0)**(-3+2)
     gaunt(-3, -1, 4) =  sqrt(54.0_dp/1089.0_dp); gaunt(-1, -3, 4) = gaunt(-3, -1, 4) * (-1.0)**(-3+1)
     gaunt(-3,  0, 4) = -sqrt(63.0_dp/1089.0_dp); gaunt( 0, -3, 4) = gaunt(-3,  0, 4) * (-1.0)**(-3-0)
     gaunt(-3,  1, 4) =  sqrt(42.0_dp/1089.0_dp); gaunt( 1, -3, 4) = gaunt(-3,  1, 4) * (-1.0)**(-3-1)
     gaunt(-2, -2, 4) = -sqrt(49.0_dp/1089.0_dp)
     gaunt(-2, -1, 4) =  sqrt(32.0_dp/1089.0_dp); gaunt(-1, -2, 4) = gaunt(-2, -1, 4) * (-1.0)**(-2+1)
     gaunt(-2,  0, 4) = -sqrt( 3.0_dp/1089.0_dp); gaunt( 0, -2, 4) = gaunt(-2,  0, 4) * (-1.0)**(-2-0)
     gaunt(-2,  1, 4) = -sqrt(14.0_dp/1089.0_dp); gaunt( 1, -2, 4) = gaunt(-2,  1, 4) * (-1.0)**(-2-1)
     gaunt(-2,  2, 4) =  sqrt(70.0_dp/1089.0_dp); gaunt( 2, -2, 4) = gaunt(-2,  2, 4) * (-1.0)**(-2-2)
     gaunt(-1, -1, 4) =  sqrt( 1.0_dp/1089.0_dp)
     gaunt(-1,  0, 4) =  sqrt(15.0_dp/1089.0_dp); gaunt( 0, -1, 4) = gaunt(-1,  0, 4) * (-1.0)**(-1-0)
     gaunt(-1,  1, 4) = -sqrt(40.0_dp/1089.0_dp); gaunt( 1, -1, 4) = gaunt(-1,  1, 4) * (-1.0)**(-1-1)
     gaunt( 0,  0, 4) =  sqrt(36.0_dp/1089.0_dp)
     gaunt( 1,  0, 4) =  sqrt(15.0_dp/1089.0_dp); gaunt( 0,  1, 4) = gaunt( 1,  0, 4) * (-1.0)**( 1-0)
     gaunt( 1,  1, 4) =  sqrt( 1.0_dp/1089.0_dp)
     gaunt( 2, -1, 4) = -sqrt(14.0_dp/1089.0_dp); gaunt(-1,  2, 4) = gaunt( 2, -1, 4) * (-1.0)**( 2+1)
     gaunt( 2,  0, 4) = -sqrt( 3.0_dp/1089.0_dp); gaunt( 0,  2, 4) = gaunt( 2,  0, 4) * (-1.0)**( 2-0)
     gaunt( 2,  1, 4) =  sqrt(32.0_dp/1089.0_dp); gaunt( 1,  2, 4) = gaunt( 2,  1, 4) * (-1.0)**( 2-1)
     gaunt( 2,  2, 4) = -sqrt(49.0_dp/1089.0_dp)
     gaunt( 3, -1, 4) =  sqrt(42.0_dp/1089.0_dp); gaunt(-1,  3, 4) = gaunt( 3, -1, 4) * (-1.0)**( 3+1)
     gaunt( 3,  0, 4) = -sqrt(63.0_dp/1089.0_dp); gaunt( 0,  3, 4) = gaunt( 3,  0, 4) * (-1.0)**( 3-0)
     gaunt( 3,  1, 4) =  sqrt(54.0_dp/1089.0_dp); gaunt( 1,  3, 4) = gaunt( 3,  1, 4) * (-1.0)**( 3-1)
     gaunt( 3,  2, 4) = -sqrt(30.0_dp/1089.0_dp); gaunt( 2,  3, 4) = gaunt( 3,  2, 4) * (-1.0)**( 3-2)
     gaunt( 3,  3, 4) =  sqrt( 9.0_dp/1089.0_dp)
  
     gaunt(-3, -3, 6) = -sqrt(   25.0_dp/184041_dp)
     gaunt(-3, -2, 6) =  sqrt(  175.0_dp/184041_dp); gaunt(-2, -3, 6) = gaunt(-3, -2, 6) * (-1.0)**(-3+2)
     gaunt(-3, -1, 6) = -sqrt(  700.0_dp/184041_dp); gaunt(-1, -3, 6) = gaunt(-3, -1, 6) * (-1.0)**(-3+1)
     gaunt(-3,  0, 6) =  sqrt( 2100.0_dp/184041_dp); gaunt( 0, -3, 6) = gaunt(-3,  0, 6) * (-1.0)**(-3-0)
     gaunt(-3,  1, 6) = -sqrt( 5250.0_dp/184041_dp); gaunt( 1, -3, 6) = gaunt(-3,  1, 6) * (-1.0)**(-3-1)
     gaunt(-3,  2, 6) =  sqrt(11550.0_dp/184041_dp); gaunt( 2, -3, 6) = gaunt(-3,  2, 6) * (-1.0)**(-3-2)
     gaunt(-3,  3, 6) = -sqrt(23100.0_dp/184041_dp); gaunt( 3, -3, 6) = gaunt(-3,  3, 6) * (-1.0)**(-3-3)
     gaunt(-2, -2, 6) =  sqrt(  900.0_dp/184041_dp)
     gaunt(-2, -1, 6) = -sqrt( 2625.0_dp/184041_dp); gaunt(-1, -2, 6) = gaunt(-2, -1, 6) * (-1.0)**(-2+1)
     gaunt(-2,  0, 6) =  sqrt( 5600.0_dp/184041_dp); gaunt( 0, -2, 6) = gaunt(-2,  0, 6) * (-1.0)**(-2-0)
     gaunt(-2,  1, 6) = -sqrt( 9450.0_dp/184041_dp); gaunt( 1, -2, 6) = gaunt(-2,  1, 6) * (-1.0)**(-2-1)
     gaunt(-2,  2, 6) =  sqrt(12600.0_dp/184041_dp); gaunt( 2, -2, 6) = gaunt(-2,  2, 6) * (-1.0)**(-2-2)
     gaunt(-1, -1, 6) = -sqrt( 5625.0_dp/184041_dp)
     gaunt(-1,  0, 6) =  sqrt( 8750.0_dp/184041_dp); gaunt( 0, -1, 6) = gaunt(-1,  0, 6) * (-1.0)**(-1-0)
     gaunt(-1,  1, 6) = -sqrt(10500.0_dp/184041_dp); gaunt( 1, -1, 6) = gaunt(-1,  1, 6) * (-1.0)**(-1-1)
     gaunt( 0,  0, 6) =  sqrt(10000.0_dp/184041_dp)
     gaunt( 1,  0, 6) =  sqrt( 8750.0_dp/184041_dp); gaunt( 0,  1, 6) = gaunt( 1,  0, 6) * (-1.0)**( 1-0)
     gaunt( 1,  1, 6) = -sqrt( 5625.0_dp/184041_dp)
     gaunt( 2, -1, 6) = -sqrt( 9450.0_dp/184041_dp); gaunt(-1,  2, 6) = gaunt( 2, -1, 6) * (-1.0)**( 2+1)
     gaunt( 2,  0, 6) =  sqrt( 5600.0_dp/184041_dp); gaunt( 0,  2, 6) = gaunt( 2,  0, 6) * (-1.0)**( 2-0)
     gaunt( 2,  1, 6) = -sqrt( 2625.0_dp/184041_dp); gaunt( 1,  2, 6) = gaunt( 2,  1, 6) * (-1.0)**( 2-1)
     gaunt( 2,  2, 6) =  sqrt(  900.0_dp/184041_dp)
     gaunt( 3, -2, 6) =  sqrt(11550.0_dp/184041_dp); gaunt(-2,  3, 6) = gaunt( 3, -2, 6) * (-1.0)**( 3+2)
     gaunt( 3, -1, 6) = -sqrt( 5250.0_dp/184041_dp); gaunt(-1,  3, 6) = gaunt( 3, -1, 6) * (-1.0)**( 3+1)
     gaunt( 3,  0, 6) =  sqrt( 2100.0_dp/184041_dp); gaunt( 0,  3, 6) = gaunt( 3,  0, 6) * (-1.0)**( 3-0)
     gaunt( 3,  1, 6) = -sqrt(  700.0_dp/184041_dp); gaunt( 1,  3, 6) = gaunt( 3,  1, 6) * (-1.0)**( 3-1)
     gaunt( 3,  2, 6) =  sqrt(  175.0_dp/184041_dp); gaunt( 2,  3, 6) = gaunt( 3,  2, 6) * (-1.0)**( 3-2)
     gaunt( 3,  3, 6) = -sqrt(   25.0_dp/184041_dp)
  
     return
  end subroutine atomic_gaunt_7band
