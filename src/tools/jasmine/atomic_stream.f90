!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_config
!!!           atomic_check_config
!!!           atomic_read_cmat
!!!           atomic_read_emat
!!!           atomic_read_tmat
!!!           atomic_make_spmat
!!!           atomic_make_fock
!!!           atomic_make_natural
!!!           atomic_2natural_case1
!!!           atomic_2natural_case2
!!!           atomic_2natural_case3
!!!           atomic_2natural_case4
!!! source  : atomic_stream.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!!           10/20/2014 by li huang
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> atomic_config: read config parameters from file atom.config.in
  subroutine atomic_config()
     use constants, only : dp, mytmp
     use parser, only : p_create, p_destroy, p_parse, p_get

     use control ! ALL

     implicit none

! local variables
! file status
     logical :: exists

! setup default values
     itask  = 1           ! type of task
     ictqmc = 1           ! type of CTQMC algorithm
     icu    = 1           ! type of Coulomb interaction
     icf    = 0           ! type of crystal field
     isoc   = 0           ! type of spin-orbital coupling (SOC)

     nband = 1            ! number of bands
     nspin = 2            ! number of spins
     norbs = nband*nspin  ! number of orbits
     ncfgs = 2**norbs     ! number of many-body configurations

     Uc = 2.00_dp         ! intraorbital Coulomb interaction
     Uv = 2.00_dp         ! interorbital Coulomb interaction
     Jz = 0.00_dp         ! Hund's exchange interaction
     Js = 0.00_dp         ! spin-flip interaction
     Jp = 0.00_dp         ! pair-hopping interaction

     Ud = 2.00_dp         ! Ud
     Jh = 0.00_dp         ! Jh
     F0 = 0.00_dp         ! F0
     F2 = 0.00_dp         ! F2
     F4 = 0.00_dp         ! F4
     F6 = 0.00_dp         ! F6

     mune   = 0.00_dp     ! chemical potential
     lambda = 0.00_dp     ! spin-orbit coupling parameter

! read in input file if possible
! reset file status
     exists = .false.

! inquire the input file status: atomic.config.in
     inquire( file = "atom.config.in", exist = exists )

! read parameters from atom.config.in
     if ( exists .eqv. .true. ) then
! create the file parser
         call p_create()

! parse the config file
         call p_parse('atom.config.in')

! extract parameters
         call p_get('itask' ,  itask)
         call p_get('ictqmc', ictqmc)
         call p_get('icu'   ,    icu)
         call p_get('icf'   ,    icf)
         call p_get('isoc'  ,   isoc)

         call p_get('nband' ,  nband)
         call p_get('nspin' ,  nspin) ! not useful
         call p_get('norbs' ,  norbs) ! not useful
         call p_get('ncfgs' ,  ncfgs) ! not useful

         call p_get('Uc'    ,     Uc)
         call p_get('Uv'    ,     Uv)
         call p_get('Jz'    ,     Jz)
         call p_get('Js'    ,     Js)
         call p_get('Jp'    ,     Jp)

         call p_get('Ud'    ,     Ud)
         call p_get('Jh'    ,     Jh)
         call p_get('F0'    ,     F0) ! not useful
         call p_get('F2'    ,     F2) ! not useful
         call p_get('F4'    ,     F4) ! not useful
         call p_get('F6'    ,     F6) ! not useful

         call p_get('mune'  ,   mune)
         call p_get('lambda', lambda)

! destroy the parser
         call p_destroy()

! calculate the norbs and ncfgs
         norbs = nband * nspin
         ncfgs = 2 ** norbs

! calculate F0, F2, F4, F6 here
         F0 = Ud
         select case (nband)
             case (5)
                 F2 = Jh * 14.0_dp / 1.625_dp
                 F4 = 0.625_dp * F2

             case (7)
                 F2 = Jh * 6435.0_dp / (286.0_dp + (195.0_dp * 451.0_dp / 675.0_dp) &
                         + (250.0_dp * 1001.0_dp / 2025.0_dp))
                 F4 = 451.0_dp / 675.0_dp * F2
                 F6 = 1001.0_dp / 2025.0_dp * F2

         end select
     else
         call s_print_error('atomic_config', 'file atom.config.in does not exist!')
     endif ! back if ( exists .eqv. .true. ) block

     return
  end subroutine atomic_config

!!>>> atomic_check_config: check the validity of input config parameters
  subroutine atomic_check_config()
     use constants, only : mystd, zero

     use control ! ALL

! local variables
! status flag for whether all of the parameters are OK
     logical :: lpass

! initialize lpass
     lpass = .true.

! check itask
     if ( itask /= 1 .and. itask /= 2 ) then
         write(mystd,'(2X,a)') 'ERROR: itask must be 1 or 2!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( itask /= 1 .and. itask /= 2 ) block

! check ictqmc
     if ( ictqmc /= 1 .and. ictqmc /= 2 .and. ictqmc /= 3 .and. ictqmc /= 4 .and. ictqmc /= 5 ) then
         write(mystd,'(2X,a)') 'ERROR: ictqmc must be one of 1, 2, 3, 4, 5!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ictqmc /= 1 .and. ictqmc /= 2 .and. ictqmc /= 3 .and. ictqmc /= 4 .and. ictqmc /= 5 ) block

     if ( ictqmc == 3 .and. isoc == 1 ) then
         write(mystd,'(2X,a)') 'ERROR: GQNs (N,Sz) algorithm is NOT supported for SOC case!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ictqmc == 3 .and. isoc == 1 ) block

     if ( ictqmc == 4 .and. isoc == 1 ) then
         write(mystd,'(2X,a)') 'ERROR: GQNs (N,Sz,Ps) algorithm is NOT supported for SOC case!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ictqmc == 4 .and. isoc == 1 ) block

     if ( ictqmc == 4 .and. icu == 2 ) then
         write(mystd,'(2X,a)') 'ERROR: GQNs (N,Sz,Ps) algorithm is NOT supported for Slater-Cordon type interaction U!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ictqmc == 4 .and. icu == 2 ) block

     if ( ictqmc == 5 .and. isoc == 0 ) then
         write(mystd,'(2X,a)') 'ERROR: GQNs (N,Jz) algorithm is ONLY supported for SOC case!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ictqmc == 5 .and. isoc == 0 ) block

     if ( ictqmc == 5 .and. isoc == 1 .and. icf /= 0 ) then
         write(mystd,'(2X,a)') 'ERROR: GQNs (N,Jz) algorithm is NOT supported for SOC plus CF case!'
         write(mystd, *)
         lpass = .false.
     endif ! back if ( ictqmc == 5 .and. isoc == 1 .and. icf /= 0 ) block

! check icu
     if ( icu /= 1 .and. icu /= 2 ) then
         write(mystd,'(2X,a)') 'ERROR: icu must be 1 or 2!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( icu /= 1 .and. icu /= 2 ) block

     if ( icu == 2 .and. nband /= 5 .and. nband /= 7 ) then
         write(mystd,'(2X,a)') 'ERROR: only support Slater-Cordon type Coulomb interaction for nband=5 or 7!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( icu == 2 .and. nband /= 5 .and. nband /= 7 ) block

! check icf
     if ( icf /= 0 .and. icf /= 1 .and. icf /= 2 ) then
         write(mystd,'(2X,a)') 'ERROR: icf must be one of 0, 1, 2!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( icf /= 0 .and. icf /= 1 .and. icf /= 2 ) block

! check isoc
     if ( isoc /= 0 .and. isoc /= 1 ) then
         write(mystd,'(2X,a)') 'ERROR: isoc must be 0 or 1!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( isoc /= 0 .and. isoc /= 1 ) block

     if ( isoc == 1 .and. nband /= 3 .and. nband /= 5 .and. nband /= 7 ) then
         write(mystd,'(2X,a)') 'ERROR: only support SOC for nband=3, 5, or 7!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( isoc == 1 .and. nband /= 3 .and. nband /= 5 .and. nband /= 7 ) block

! check nband
     if ( nband <= 0 ) then
         write(mystd,'(2X,a)') 'ERROR: number of bands must be larger than zero!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( nband <= 0 ) block

! check Uc, Uv, Jz, Js, Jp, Ud, JH
     if ( Uc < zero .or. Uv < zero ) then
         write(mystd,'(2X,a)') 'ERROR: Uc and Uv must be larger than zero!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( Uc < zero .or. Uv < zero ) block

     if ( Jz < zero .or. Js < zero .or. Jp < zero ) then
         write(mystd,'(2X,a)') 'ERROR: Jz, Js, and Jp must be larger than zero!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( Jz < zero .or. Js < zero .or. Jp < zero ) block

     if ( Ud < zero .or. Jh < zero ) then
         write(mystd,'(2X,a)') 'ERROR: Ud and Jh must be larger than zero!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( Ud < zero .or. Jh < zero ) block

     if ( lpass .eqv. .false. ) then
         call s_print_error('atomic_check_config','invalid parameters found in atom.config.in file!')
     endif ! back if ( lpass .eqv. .false. ) block

     return
  end subroutine atomic_check_config

!!>>> atomic_read_cmat: read crystal field from file atomic.cf.in
  subroutine atomic_read_cmat()
     use, intrinsic :: iso_fortran_env, only : iostat_end
     use constants, only : dp, zero, mytmp

     use m_spmat, only : cmat

     implicit none

! local variables
! file status
     logical  :: exists

! iostat
     integer  :: ierr

! dummy variables
     integer  :: i
     integer  :: j
     real(dp) :: raux

! we shall read crystal field cmat from file atom.cf.in
! inquire file at first
     inquire(file = 'atom.cf.in', exist = exists)

     if ( exists .eqv. .true. ) then
! open file atom.cf.in
         open(mytmp, file='atom.cf.in', form='formatted', status='unknown')

! read the data until EOF
         do
             read(mytmp, *, iostat = ierr) i, j, raux
             if ( ierr == iostat_end ) EXIT
! crystal field is actually real
             cmat(i,j) = dcmplx(raux, zero)
         enddo ! over do while loop

! close data file
         close(mytmp)
     else
         call s_print_error('atomic_read_cmat', 'file atomic.cf.in does not exist!')
     endif ! back if ( exists .eqv. .true. ) block

     return
  end subroutine atomic_read_cmat

!!>>> atomic_read_emat: read onsite impurity level from file atomic.eimp.in
  subroutine atomic_read_emat()
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use m_spmat, only : emat

     implicit none

! local variables
! file status
     logical  :: exists

! loop index
     integer  :: i

! dummy variables
     integer  :: i1
     integer  :: i2
     real(dp) :: raux

! we shall read emat from file atomic.eimp.in
! inquire file at first
     inquire(file = 'atom.eimp.in', exist = exists)

     if ( exists .eqv. .true. ) then
! open file atom.eimp.in
         open(mytmp, file='atom.eimp.in', form='formatted', status='unknown')

! read the data file
         do i=1,norbs
             read(mytmp, *) i1, i2, raux
! emat is actually real in natural basis
             emat(i,i) = dcmplx(raux, zero)
         enddo ! over i={1,norbs} loop

! close data file
         close(mytmp)
     else
         call s_print_error('atomic_read_emat', 'file atomic.eimp.in does not exist!')
     endif ! back if ( exists .eqv. .true. ) block

     return
  end subroutine atomic_read_emat

!!>>> atomic_read_tmat: read the transformation matrix tmat from file atomic.tmat.in
  subroutine atomic_read_tmat()
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use m_spmat, only : tmat

     implicit none

! local variables
! file status
     logical :: exists

! loop index
     integer :: i
     integer :: j

! dummy variables
     integer :: i1
     integer :: i2
     real(dp) :: raux

! we shall read transformation matrix tmat from file atomic.umat.in
! inquire file at first
     inquire(file = 'atom.tmat.in', exist = exists)

     if ( exists .eqv. .true. ) then
! open file atom.tmat.in
         open(mytmp, file='atom.tmat.in', form='formatted', status='unknown')

! read the data file
         do i=1,norbs
             do j=1,norbs
                 read(mytmp,*) i1, i2, raux
                 tmat(j,i) = dcmplx(raux, zero)
             enddo ! over j={1,norbs} loop
         enddo ! over i={1,norbs} loop

! close data file
         close(mytmp)
     else
         call s_print_error('atomic_read_tmat', 'file atomic.tmat.in does not exist')
     endif ! back if ( exists .eqv. .true. ) block

     return
  end subroutine atomic_read_tmat

!!>>> atomic_make_spmat: make single particle related matrices, including
!!>>> crystal field (CF), spin-orbit coupling (SOC), and Coulomb interaction
!!>>> U tensor
  subroutine atomic_make_spmat()
     use constants, only : two, czero

     use control, only : itask
     use control, only : icu, icf, isoc
     use control, only : nband
     use control, only : lambda
     use m_spmat, only : cmat, smat

     implicit none

! make crystal field and spin-orbital coupling
! method 1: make them inside
     if ( itask == 1 ) then
! 1A: make crysal field
! we read the non-zero elements of crystal field from file atom.cf.in.
! the crystal field is defined on real orbital basis. at present, we only
! support real crystal field, so, the elements in this file provided by
! user must be real
         if ( icf > 0 ) then
             call atomic_read_cmat()
         else
             cmat = czero
         endif ! back if ( icf > 0 ) block

! 1B: make spin-orbit coupling
! make an atomic on-site SOC, $\lambda * L * S$
! it is defined on the complex orbital basis
         if ( isoc > 0 ) then
             select case (nband)

                 case (3)
                     call atomic_make_smat3(smat)
! for 3 bands system, there is a minus sign
                     smat = -smat * lambda / two

                 case (5)
                     call atomic_make_smat5(smat)
                     smat = smat * lambda / two

                 case (7)
                     call atomic_make_smat7(smat)
                     smat = smat * lambda / two

                 case default
                     call s_print_error('atomic_make_spmat', 'not implemented!')

             end select
         else
             smat = czero
         endif ! back if ( isoc > 0 ) block
! method 2: make them outside
     else
! read the emat (CF+SOC) matrices on natural basis, this matrix should be
! a diagonal matrix, and the elements must be real
         call atomic_read_emat()

! read the transformation matrices used to transfer emat from original
! basis to natural basis. without SOC, the original basis is the real
! orbital basis. with SOC, the original basis is the complex orbital basis
! at present, we just only support real numbers of this tmat
         call atomic_read_tmat()
     endif ! back if ( itask == 1 ) block

! make Coulomb interaction U
     if ( icu == 1 ) then
! Kanamori parameters type
! it is defined on real orbital basis
         call atomic_make_umatK()
     else
! Slater-Cordon parameters type
! it is defined on complex orbital basis
         call atomic_make_umatS()
     endif ! back if ( icu == 1 ) block

     return
  end subroutine atomic_make_spmat

!!>>> atomic_make_fock: make Fock basis for full Hilbert space
  subroutine atomic_make_fock()
     use control, only : norbs, ncfgs
     use m_full, only : dim_sub_n, dec_basis, bin_basis, index_basis

     implicit none

! local variables
! loop index
     integer :: i, j, k

! basis counter
     integer :: basis_count

! number of electrons for Fock state
     integer :: nelec

! initialize them
     dim_sub_n = 0
     dec_basis = 0
     bin_basis = 0
     index_basis = 0

! it is a number of combination C_{norbs}^{i}
     do i=0,norbs
         call s_combination(i, norbs, dim_sub_n(i))
     enddo

! construct decimal form and index of Fock basis
     basis_count = 0
     do i=0, norbs
         do j=0, 2**norbs-1
             nelec = 0
             do k=1,norbs
                 if( btest(j, k-1) ) nelec = nelec + 1
             enddo
             if ( nelec .eq. i ) then
                 basis_count = basis_count + 1
                 dec_basis(basis_count) = j
                 index_basis(j) = basis_count
             endif
         enddo
     enddo

! construct binary form of Fock basis
     do i=1,ncfgs
         do j=1,norbs
             if( btest(dec_basis(i), j-1) ) bin_basis(j, i) = 1
         enddo
     enddo

! dump Fock basis to file "atom.basis.dat" for reference
     call atomic_dump_fock()

     return
  end subroutine atomic_make_fock

!!>>> make natural basis, the natural basis is the basis on which the
!!>>> impurity energy matrix is diagonal
  subroutine atomic_make_natural()
     use constants, only : dp, czero, mystd
     use control, only : norbs, itask, icf, isoc, icu
     use m_spmat, only : umat, tmat

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
     call atomic_dump_emat()

! we need transform Coulomb interaction U
! for non-soc case, the tran_umat is defined as from real orbital basis to natural basis
     if (isoc==0) then
! for Slater-Cordon parameters Coulomb interaction U
! we first need to transfrom umat from complex orbital basis to real orbital basis
         if ( icu == 2 ) then
             call atomic_make_tmat_c2r( umat_c2r )
             call atomic_tran_umat( umat_c2r, umat, tmp_mat )
             umat = tmp_mat
         endif
! for soc case, the tran_umat is defined as from complex orbital basis to natural basis
     else
! for Kanamori parameters Coulomb interaction U
! we first need to transfrom umat from real orbital basis to complex orbital basis
         if ( icu == 1 ) then
             call atomic_make_tmat_r2c( umat_r2c )
             call atomic_tran_umat( umat_r2c, umat, tmp_mat )
             umat = tmp_mat
         endif
     endif

! finally, transform umat to natural basis
     call atomic_tran_umat(tmat, umat, tmp_mat)
     umat = tmp_mat

     return
  end subroutine atomic_make_natural

!!>>> atomic_2natural_case1: make natural basis for no crystal or
!!>>> diagonal crystal without spin-orbital coupling
  subroutine atomic_2natural_case1()
     use constants, only : mystd, zero, cone
     use control, only : norbs, mune
     use m_spmat, only : cmat, emat, tmat

     implicit none

! local variables
     integer :: i

! set emat
     emat = cmat

! add chemical potential to eimpmat
     do i=1, norbs
         emat(i,i) = emat(i,i) + mune
     enddo
! for this case, the natural basis is the real orbital basis
! so, the tran_umat is a unity matrix
     tmat = zero
     do i=1, norbs
         tmat(i,i) = cone
     enddo

     write(mystd,'(2X,a)') 'jasmine >>> natural basis is: real orbital basis'
     write(mystd, *)

     call atomic_dump_natural('# natural basis is real orbital, umat: real to natural')

     return
  end subroutine atomic_2natural_case1

!!>>> atomic_2natural_case2: make natural basis for non-diagonal
!!>>> crystal field without spin-orbital coupling
  subroutine atomic_2natural_case2()
     use constants, only : mystd, dp, czero
     use control, only : norbs, nband, mune
     use m_spmat, only : cmat, emat, tmat

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
     integer :: i, j

! set eimpmat
     emat = cmat

! get eimp for nospin freedom
     do i=1, norbs/2
         do j=1, norbs/2
             eimp_nospin(j,i) = emat(2*j-1,2*i-1)
         enddo
     enddo

! diagonalize eimp_nospin to get natural basis
     eimp_nospin_real = real(eimp_nospin)
     call s_eig_sy(nband, nband, eimp_nospin_real, eigval, eigvec)

     eimp_nospin = czero
     do i=1, nband
         eimp_nospin(i,i) = eigval(i)
     enddo
     umat_nospin = dcmplx(eigvec)

     do i=1, nband
         do j=1, nband
             emat(2*j-1,2*i-1) = eimp_nospin(j,i)
             emat(2*j,2*i)     = eimp_nospin(j,i)
             tmat(2*j-1,2*i-1) = umat_nospin(j,i)
             tmat(2*j,2*i)     = umat_nospin(j,i)
         enddo
     enddo

! add chemical potential to eimpmat
     do i=1, norbs
         emat(i,i) = emat(i,i) + mune
     enddo

     write(mystd, '(2X,a)') 'jasmine >>> natural basis is: linear combination of real orbitals '
     write(mystd, *)

     call atomic_dump_natural('# natural basis is linear combination of real orbitals, umat: real to natural')

     return
  end subroutine atomic_2natural_case2

!!>>> atomic_2natural_case3: make natural basis for the case without
!!>>> crystal field and with spin-orbital coupling
!!>>> for this special case, the natural basis is |J^2,Jz>
  subroutine atomic_2natural_case3()
     use constants, only : dp, mystd
     use control, only : norbs, mune
     use m_spmat, only : emat, smat, tmat

     implicit none

! local variables
! umat from complex orbital basis to |j^2, jz> basis
     complex(dp) :: umat_c2j(norbs, norbs)

! loop inex
     integer :: i

! set eimpmat
     emat = smat

     call atomic_make_tmat_c2j(umat_c2j)

! for soc case, the tran_umat is from complex orbital basis to natural basis
     tmat = umat_c2j

! transform sp_eimp_mat to natural basis
     call atomic_tran_repr_cmpl(norbs, emat, tmat)

! add chemical potential to eimpmat
     do i=1, norbs
         emat(i,i) = emat(i,i) + mune
     enddo

     write(mystd, '(2X,a)') 'jasmine >>> natural basis is: |j2,jz> '
     write(mystd, *)

     call atomic_dump_natural('# natural basis is |j2,jz>, umat: complex to natural')

     return
  end subroutine atomic_2natural_case3

!!>>> atomic_2natural_case4: make natural basis for the case with
!!>>> crystal field and with spin-orbital coupling
  subroutine atomic_2natural_case4()
     use constants, only : dp, mystd, eps6
     use control, only : norbs, mune

     use m_spmat, only : cmat, smat, emat, tmat

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
!     logical :: lreal

! get umat_r2c
     call atomic_make_tmat_r2c(umat_r2c)

! transfrom cfmat to complex orbital basis
     call atomic_tran_repr_cmpl(norbs, cmat, umat_r2c)

! check whether cfmat is real, if not, we cann't make natural basis
     if ( any( abs( aimag(cmat) ) > eps6 ) ) then
         call s_print_error('atomic_2natural_case4', 'crystal field on &
             complex orbital basis is not real, cannot make natural basis !')
     endif

! set eimpmat
     emat = smat + cmat

     tmp_mat = real(emat)

     call s_eig_sy(norbs, norbs, tmp_mat, eigval, eigvec)

     umat_c2n = eigvec
     tmat = umat_c2n

! transform eimpmat to natural basis
     call atomic_tran_repr_cmpl(norbs, emat, umat_c2n)

! add chemical poential to eimpmat
     do i=1, norbs
         emat(i,i) = emat(i,i) + mune
     enddo

     write(mystd, '(2X,a)') 'jasmine >>> natural basis is: linear combination of complex orbitals '
     write(mystd, *)

     call atomic_dump_natural('# natural basis is linear combination of complex orbitals, umat: complex to natural')

     return
  end subroutine atomic_2natural_case4
