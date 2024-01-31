!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : atomic_setup_param
!!!           atomic_check_param
!!!           atomic_input_cmat
!!!           atomic_input_emat
!!!           atomic_input_tmat
!!!           atomic_build_fock
!!!           atomic_build_spmat
!!!           atomic_build_natural
!!!           atomic_alloc_array
!!!           atomic_final_array
!!! source  : atomic_stream.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           01/31/2024 by li huang (last modified)
!!! purpose : read essential input data from the external files, build
!!!           the Fock basis, construct single particle matrices and
!!!           natural eigenbasis, etc.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> config atomic eigenvalue problem solver                          <<<
!!========================================================================

!!
!! @sub atomic_setup_param
!!
!! read control parameters from file atom.config.in
!!
  subroutine atomic_setup_param()
     use constants, only : dp

     use parser, only : p_create
     use parser, only : p_parse
     use parser, only : p_get
     use parser, only : p_destroy

     use control ! ALL

     implicit none

!! local variables
     ! file status, if the file atom.config.in exists
     logical :: exists

!! [body

     ! setup general control flags
     !--------------------------------------------------------------------
     ibasis = 1           ! source of the natural eigenbasis
     ictqmc = 1           ! how to diagonalize atomic Hamiltonian
     icu    = 1           ! type of Coulomb interaction
     icf    = 0           ! type of crystal field splitting (CFS)
     isoc   = 0           ! type of spin-orbit coupling (SOC)
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     ! setup common variables for atomic Hamiltonian
     !--------------------------------------------------------------------
     nband  = 1           ! number of bands
     nspin  = 2           ! number of spins
     norbs  = nband*nspin ! number of orbitals
     ncfgs  = 2**norbs    ! number of many-body configurations
     nmini  = 0           ! minimal total occupancy N to be kept
     nmaxi  = norbs       ! maximal total occupancy N to be kept
     !--------------------------------------------------------------------
     Uc     = 2.00_dp     ! intraorbital Coulomb interaction
     Uv     = 2.00_dp     ! interorbital Coulomb interaction
     Jz     = 0.00_dp     ! Hund's exchange interaction
     Js     = 0.00_dp     ! spin-flip interaction
     Jp     = 0.00_dp     ! pair-hopping interaction
     !--------------------------------------------------------------------
     Ud     = 2.00_dp     ! Coulomb interaction parameter
     Jh     = 0.00_dp     ! Hund's exchange parameter
     !--------------------------------------------------------------------
     mune   = 0.00_dp     ! chemical potential
     lambda = 0.00_dp     ! spin-orbit coupling strength
     !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     ! read the input file if available
     ! reset file status
     exists = .false.

     ! inquire the status of input file: atomic.config.in
     inquire( file = "atom.config.in", exist = exists )

     ! read control parameters from atom.config.in
     if ( exists .eqv. .true. ) then
         ! create the file parser
         call p_create()

         ! parse the config file
         call p_parse('atom.config.in')

         ! extract parameters
         call p_get('ibasis', ibasis)
         call p_get('ictqmc', ictqmc)
         call p_get('icu'   ,    icu)
         call p_get('icf'   ,    icf)
         call p_get('isoc'  ,   isoc)

         call p_get('nband' ,  nband)
         call p_get('nspin' ,  nspin) ! not useful
         call p_get('norbs' ,  norbs) ! not useful
         call p_get('ncfgs' ,  ncfgs) ! not useful

         ! calculate norbs and ncfgs
         norbs = nband * nspin
         ncfgs = 2 ** norbs

         ! setup nmini and nmaxi
         nmini = 0
         nmaxi = norbs

         ! continue to extract parameters
         call p_get('nmini' ,  nmini)
         call p_get('nmaxi' ,  nmaxi)

         call p_get('Uc'    ,     Uc)
         call p_get('Uv'    ,     Uv)
         call p_get('Jz'    ,     Jz)
         call p_get('Js'    ,     Js)
         call p_get('Jp'    ,     Jp)

         call p_get('Ud'    ,     Ud)
         call p_get('Jh'    ,     Jh)

         call p_get('mune'  ,   mune)
         call p_get('lambda', lambda)

         ! destroy the parser
         call p_destroy()
     else
         call s_print_exception('atomic_setup_param', &
             & 'file atom.config.in does not exist!')
     endif ! back if ( exists .eqv. .true. ) block

!! body]

     return
  end subroutine atomic_setup_param

!!
!! @sub atomic_check_config
!!
!! check validity and consistency of the control parameters
!!
  subroutine atomic_check_param()
     use constants, only : zero
     use constants, only : mystd

     use control ! ALL

!! local variables
     ! whether all the control parameters are valid
     logical :: lpass

!! [body

     ! initialize lpass
     lpass = .true.

     ! check ibasis
     if ( ibasis < 1 .or. ibasis > 2 ) then
         write(mystd,'(2X,a)') 'ERROR: ibasis must be 1 or 2!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ibasis < 1 .or. ibasis > 2 ) block

     ! check ictqmc
     if ( ictqmc < 1 .or. ictqmc > 5 ) then
         write(mystd,'(2X,a)') 'ERROR: ictqmc must be one of 1, &
             & 2, 3, 4, 5!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ictqmc < 0 .or. ictqmc > 5 ) block
     !
     if ( ictqmc == 3 .and. isoc == 1 ) then
         write(mystd,'(2X,a)') 'ERROR: subspace diagonalization &
             & algorithm using GQNs (N,Sz) is NOT supported for &
             & SOC case!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ictqmc == 3 .and. isoc == 1 ) block
     !
     if ( ictqmc == 4 .and. isoc == 1 ) then
         write(mystd,'(2X,a)') 'ERROR: subspace diagonalization &
             & algorithm using GQNs (N,Sz,Ps) is NOT supported  &
             & for SOC case!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ictqmc == 4 .and. isoc == 1 ) block
     !
     if ( ictqmc == 4 .and. icu == 2 ) then
         write(mystd,'(2X,a)') 'ERROR: subspace diagonalization &
             & algorithm using GQNs (N,Sz,Ps) is NOT supported  &
             & for Slater-Cordon type interaction U!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ictqmc == 4 .and. icu == 2 ) block
     !
     if ( ictqmc == 5 .and. isoc == 0 ) then
         write(mystd,'(2X,a)') 'ERROR: subspace diagonalization &
             & algorithm using GQNs (N,Jz) is ONLY supported    &
             & for SOC case!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ictqmc == 5 .and. isoc == 0 ) block
     !
     if ( ictqmc == 5 .and. isoc == 1 .and. icf /= 0 ) then
         write(mystd,'(2X,a)') 'ERROR: subspace diagonalization &
             & algorithm using GQNs (N,Jz) is NOT supported for &
             & SOC + CFS case!'
         write(mystd, *)
         lpass = .false.
     endif ! back if ( ictqmc == 5 .and. isoc == 1 .and. icf /= 0 ) block

     ! check icu
     if ( icu < 1 .or. icu > 3 ) then
         write(mystd,'(2X,a)') 'ERROR: icu must be 1 or 2 or 3!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( icu < 1 .or. icu > 3 ) block
     !
     if ( icu == 2 .and. nband /= 5 .and. nband /= 7 ) then
         write(mystd,'(2X,a)') 'ERROR: Slater-Cordon type Coulomb &
             & interaction is only suitable for 5- or 7-band system!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( icu == 2 .and. nband /= 5 .and. nband /= 7 ) block
     !
     if ( icu == 3 .and. nband /= 5 ) then
         write(mystd,'(2X,a)') 'ERROR: anisotropic Hunds rule exchange &
             & in Kanamori type Coulomb interaction is only suitable   &
             & for 5-band system!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( icu == 3 .and. nband /= 5 ) block

     ! check icf
     if ( icf < 0 .or. icf > 2 ) then
         write(mystd,'(2X,a)') 'ERROR: icf must be one of 0, 1, 2!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( icf < 0 .or. icf > 2 ) block

     ! check isoc
     if ( isoc < 0 .or. isoc > 1 ) then
         write(mystd,'(2X,a)') 'ERROR: isoc must be 0 or 1!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( isoc < 0 .or. isoc > 1 ) block
     !
     if ( isoc == 1 .and. nband /= 3 .and. nband /= 5 .and. nband /= 7 ) then
         write(mystd,'(2X,a)') 'ERROR: the SOC setup is only possible &
             & for 3-, 5-, or 7-band system!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( isoc == 1 .and. nband /= 3 .and. nband /= 5 .and. nband /= 7 ) block

     ! check nband
     if ( nband <= 0 .or. nband >= 8 ) then
         write(mystd,'(2X,a)') 'ERROR: number of bands should be a &
             & positive integer (1 <= nband <= 7)!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( nband <= 0 .or. nband >= 8 ) block
     !
     if ( nband >= 5 .and. ictqmc == 1 ) then
         write(mystd,'(2X,a)') 'ERROR: when number of bands is larger &
             & than 4, the direct diagonalization algorithm is NOT    &
             & supported any more!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( nband >= 5 .and. ictqmc <= 1 ) block

     ! check nspin
     if ( nspin /= 2 ) then
         write(mystd,'(2X,a)') 'ERROR: number of spin projections must be 2!'
         write(mystd,*)
         lpass = .false.
     endif

     ! check norbs
     if ( norbs /= nspin * nband ) then
         write(mystd,'(2X,a)') 'ERROR: number of bands is not compatible &
             & with number of orbitals!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( norbs /= nspin * nband ) block

     ! check ncfgs
     if ( ncfgs /= 2 ** norbs ) then
         write(mystd,'(2X,a)') 'ERROR: number of many body states is not &
             & compatible with number of orbitals!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( ncfgs /= 2 ** norbs ) block

     ! check nmini and nmaxi
     if ( nmini < 0 ) then
         nmini = 0
         write(mystd,'(2X,a)') 'WARNING: nmini < 0, enforce to be zero!'
     endif ! back if ( nmini < 0 ) block
     !
     if ( nmaxi > norbs ) then
         nmaxi = norbs
         write(mystd,'(2X,a)') 'WARNING: nmaxi > norbs, enforce to be norbs!'
         write(mystd,*)
     endif ! back if ( nmaix > norbs ) block

     ! check Uc, Uv, Jz, Js, Jp
     if ( Uc < zero .or. Uv < zero ) then
         write(mystd,'(2X,a)') 'ERROR: Uc and Uv must be larger than &
             & or equal to zero!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( Uc < zero .or. Uv < zero ) block
     !
     if ( Jz < zero .or. Js < zero .or. Jp < zero ) then
         write(mystd,'(2X,a)') 'ERROR: Jz, Js, and Jp must be larger &
             & than or equal to zero!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( Jz < zero .or. Js < zero .or. Jp < zero ) block

     ! check Ud and Jh
     if ( Ud < zero .or. Jh < zero ) then
         write(mystd,'(2X,a)') 'ERROR: Ud and Jh must be larger than &
             & or equal to zero!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( Ud < zero .or. Jh < zero ) block

     ! check lambda
     if ( lambda < zero ) then
         write(mystd,'(2X,a)') 'ERROR: lambda must be larger than or &
             & equal to zero!'
         write(mystd,*)
         lpass = .false.
     endif ! back if ( lambda < zero ) block

     ! final assertion
     if ( lpass .eqv. .false. ) then
         call s_print_error('atomic_check_param','invalid parameters &
             & found in atom.config.in file!')
     endif ! back if ( lpass .eqv. .false. ) block

!! body]

     return
  end subroutine atomic_check_param

!!========================================================================
!!>>> setup atomic Hamiltonian                                         <<<
!!========================================================================

!!
!! @sub atomic_input_cmat
!!
!! read crystal field splitting from file atomic.cmat.in
!!
  subroutine atomic_input_cmat()
     use, intrinsic :: iso_fortran_env, only : iostat_end

     use constants, only : dp
     use constants, only : zero
     use constants, only : mytmp

     use control, only : norbs

     use m_spmat, only : cmat

     implicit none

!! local variables
     ! file status
     logical  :: exists

     ! iostat
     integer  :: ierr

     ! dummy variables
     integer  :: i1
     integer  :: i2
     real(dp) :: raux

!! [body

     ! we shall read crystal field splitting into matrix cmat from
     ! file atom.cmat.in. note that crystal field splitting could
     ! be non-diagonal
     !
     ! inquire file's status at first
     inquire( file = 'atom.cmat.in', exist = exists )
     if ( exists .eqv. .false. ) then
         call s_print_error('atomic_input_cmat', &
             & 'file atomic.cmat.in does not exist!')
     endif ! back if ( exists .eqv. .false. ) block

     ! open file atom.cmat.in
     open(mytmp, file='atom.cmat.in', form='formatted', status='unknown')

     ! read the data until EOF
     do
         read(mytmp,*,iostat = ierr) i1, i2, raux
         if ( ierr == iostat_end ) EXIT
         !
         ! crystal field splitting is actually real
         call s_assert( i1 <= norbs .and. i2 <= norbs )
         cmat(i1,i2) = dcmplx(raux, zero)
     enddo ! over do while loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_input_cmat

!!
!! @sub atomic_input_emat
!!
!! read onsite impurity level from file atomic.emat.in
!!
  subroutine atomic_input_emat()
     use constants, only : dp
     use constants, only : zero
     use constants, only : mytmp

     use control, only : norbs

     use m_spmat, only : emat

     implicit none

!! local variables
     ! file status
     logical  :: exists

     ! loop index
     integer  :: i

     ! dummy variables
     integer  :: i1
     integer  :: i2
     real(dp) :: raux

!! [body

     ! we shall read onsite impurity level into matrix emat from
     ! file atomic.emat.in. note that emat is actually real and
     ! diagonal in natural eigenbasis
     !
     ! inquire file's status at first
     inquire( file = 'atom.emat.in', exist = exists )
     if ( exists .eqv. .false. ) then
         call s_print_error('atomic_input_emat', &
             & 'file atomic.emat.in does not exist!')
     endif ! back if ( exists .eqv. .false. ) block

     ! open file atom.emat.in
     open(mytmp, file='atom.emat.in', form='formatted', status='unknown')

     ! read the data file
     do i=1,norbs
         read(mytmp,*) i1, i2, raux
         emat(i,i) = dcmplx(raux, zero)
     enddo ! over i={1,norbs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_input_emat

!!
!! @sub atomic_input_tmat
!!
!! read the transformation matrix tmat from file atomic.tmat.in
!!
  subroutine atomic_input_tmat()
     use constants, only : dp
     use constants, only : zero
     use constants, only : mytmp

     use control, only : norbs
     use m_spmat, only : tmat

     implicit none

!! local variables
     ! file status
     logical :: exists

     ! loop index
     integer :: i
     integer :: j

     ! dummy variables
     integer :: i1
     integer :: i2
     real(dp) :: raux

!! [body

     ! we shall read transformation matrix tmat from file atomic.tmat.in.
     ! it is used to transform the single particle matrices from original
     ! basis to natural eigenbasis 
     !
     ! inquire file's status at first
     inquire( file = 'atom.tmat.in', exist = exists )
     if ( exists .eqv. .false. ) then
         call s_print_error('atomic_input_tmat', &
             & 'file atomic.tmat.in does not exist')
     endif ! back if ( exists .eqv. .false. ) block

     ! open file atom.tmat.in
     open(mytmp, file='atom.tmat.in', form='formatted', status='unknown')

     ! read the data file
     do i=1,norbs
         do j=1,norbs
             read(mytmp,*) i1, i2, raux
             ! tmat is actually real
             tmat(i,j) = dcmplx(raux, zero)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_input_tmat

!!========================================================================
!!>>> build basis for atomic Hamiltonian                               <<<
!!========================================================================

!!
!! @sub atomic_build_fock
!!
!! make Fock basis in the full Hilbert space
!!
  subroutine atomic_build_fock()
     use constants, only : mystd

     use control, only : norbs, ncfgs

     use m_fock, only : dim_sub_n
     use m_fock, only : bin_basis
     use m_fock, only : dec_basis
     use m_fock, only : ind_basis

     implicit none

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

     ! counter for Fock states
     integer :: state_count

     ! number of electrons for the current Fock state
     integer :: nelec

!! [body

     ! initialize them
     dim_sub_n = 0
     bin_basis = 0
     dec_basis = 0
     ind_basis = 0

     ! evaluate dim_sub_n
     ! it is a number of combination C_{norbs}^{i}
     do i=0,norbs
         call s_combination(i, norbs, dim_sub_n(i))
         !
         write(mystd,'(4X,a)', advance = 'no') 'number of Fock states: '
         write(mystd,'(2X,i4)', advance = 'no') dim_sub_n(i)
         write(mystd,'(1X,a,i2,a)') '( N = ', i, ' )'
     enddo ! over i={0,norbs} loop

     ! construct decimal form and index of Fock basis
     state_count = 0
     !
     do i=0,norbs
         do j=0,2**norbs-1 ! actually go through every Fock state
             nelec = 0
             do k=1,norbs
                 if ( btest(j, k-1) ) nelec = nelec + 1
             enddo ! over k={1,norbs} loop
             if ( nelec == i ) then
                 state_count = state_count + 1
                 dec_basis(state_count) = j
                 ind_basis(j) = state_count
             endif ! back if ( nelec == i ) block
         enddo ! over j={0,2**norbs-1} loop
     enddo ! over i={0,norbs} loop
     !
     call s_assert( state_count == ncfgs )

     ! construct binary form of Fock basis
     do i=1,ncfgs
         do j=1,norbs
             if ( btest(dec_basis(i), j-1) ) bin_basis(j,i) = 1
         enddo ! over j={1,norbs} loop
         !
         write(mystd,'(4X,a,i6)', advance = 'no') 'Fock state: ', i
         write(mystd,'(2X,a,i6)', advance = 'no') 'decimal: ', dec_basis(i)
         write(mystd,'(2X,a,*(i1))') 'binary: ', bin_basis(:,i)
     enddo ! over i={1,ncfgs} loop

     ! dump Fock basis to file atom.fock.dat for reference
     call atomic_dump_fock()

!! body]

     return
  end subroutine atomic_build_fock

!!
!! make single particle matrices, including the crystal field splitting
!! (CFS), spin-orbit coupling (SOC), and Coulomb interaction U etc.
!!
!! when constructing these matrices, we should define a single particle
!! basis at first. there are four basis sets that we adopt in the jasmine
!! code. now let us take a 5-orbitals system as an example to illustrate
!! the four basis sets.
!!
!! (1) real orbital basis (the real spherical harmonics)
!!     |dxy,up>,   |dyz,up>,   |dz2,up>,   |dxz,up>,   |dx2-y2,up>
!!     |dxy,down>, |dyz,down>, |dz2,down>, |dxz,down>, |dx2-y2,down>
!!
!! (2) complex orbital basis (the complex spherical functions)
!!     it is eigenstate of operators l^2 and l_z, |l,m,spin>
!!     for d electron system, l = 2, m = \pm 2, \pm 1, 0
!!     |2,-2,up>,   |2,-1,up>,   |2,0,up>,   |2,1,up>,   |2,2,up>
!!     |2,-2,down>, |2,-1,down>, |2,0,down>, |2,1,down>, |2,2,down>
!!
!! (3) j^2 - j_z diagonal basis
!!     it is eigenstate of operators j^2 and j_z, |j,m_j>
!!     for d electron system, j = 3/2 or 5/2, m_j = -j, -j+1, ..., j-1, j
!!     |3/2,-3/2>, |3/2,-1/2>, |3/2,1/2> |3/2,3/2>
!!     |5/2,-5/2>, |5/2,-3/2>, |5/2,-1/2> |5/2,1/2>, |5/2,3/2>, |5/2,5/2>
!!
!! (4) the natural eigenbasis, on which the onsite energy of impurity is
!!     diagonal. we have to diagonalize H_{CFS} + H_{SOC} to obtain
!!     the natural eigenbasis
!!
!! note that the CFS is always defined in real orbital basis, SOC is
!! always defined in complex orbital basis. the Coulomb interaction U
!! should be defined in real orbital basis or complex orbital basis,
!! which depends on the form of Coulomb interaction. thus, we often
!! need to transform them between two different basis sets
!!

!!
!! @sub atomic_build_spmat
!!
!! try to make various single particle matrices, including the crystal
!! field splitting (CFS), the spin-orbit coupling (SOC), and the Coulomb
!! interaction tensor (U)
!!
  subroutine atomic_build_spmat()
     use constants, only : two
     use constants, only : czero
     use constants, only : mystd

     use control, only : ibasis
     use control, only : icu, icf, isoc
     use control, only : nband
     use control, only : lambda

     use m_spmat, only : cmat
     use m_spmat, only : smat

     implicit none

!! [body

     ! make crystal field splitting and spin-orbit coupling
     !
     ! method 1: make them inside
     if ( ibasis == 1 ) then

         ! 1A: make crysal field splitting
         write(mystd,'(4X,a)') 'make crystal field splitting term'
         !
         ! we read the non-zero elements of crystal field splitting
         ! from file atom.cmat.in. the crystal field splitting must
         ! be defined on real orbital basis. so far, we only support
         ! real crystal field splitting. thus, the elements in this
         ! file provided by users must be real
         if ( icf > 0 ) then
             call atomic_input_cmat()
         else
             cmat = czero
         endif ! back if ( icf > 0 ) block

         ! 1B: make spin-orbit coupling
         write(mystd,'(4X,a)') 'make spin-orbit coupling term'
         !
         ! make an atomic spin-orbit coupling, $\lambda * L * S$
         ! it is defined on the complex orbital basis
         if ( isoc > 0 ) then
             select case (nband)

                 case (3) ! 3-band system
                     call atomic_make_smat3(smat)
                     ! for 3 band system, there is a minus sign
                     smat = -smat * lambda / two

                 case (5) ! 5-band system
                     call atomic_make_smat5(smat)
                     smat = smat * lambda / two

                 case (7) ! 7-band system
                     call atomic_make_smat7(smat)
                     smat = smat * lambda / two

                 case default
                     call s_print_error('atomic_build_spmat', &
                         & 'not implemented!')

             end select
         else
             smat = czero
         endif ! back if ( isoc > 0 ) block

     ! method 2: make them outside
     else

         ! read the matrix emat (CFS + SOC) on natural eigenbasis,
         ! this matrix must be a diagonal matrix, and the elements
         ! must be real
         write(mystd,'(4X,a)') 'make crystal field splitting + &
             & spin-orbiit coupling terms'
         !
         call atomic_input_emat()

     endif ! back if ( ibasis == 1 ) block

     ! make Coulomb interaction U
     write(mystd,'(4X,a)') 'make Coulomb interaction term'
     !
     ! Kanamori parameters type
     ! it is defined on real orbital basis
     if ( icu == 1 .or. icu == 3 ) then
     !
         call atomic_make_umatK()
     !
     ! Slater-Cordon parameters type
     ! it is defined on complex orbital basis
     else
     !
         call atomic_make_umatS()
     !
     endif ! back if ( icu == 1 .or. icu == 3 ) block

!! body]

     return
  end subroutine atomic_build_spmat

!!
!! @sub atomic_build_natural
!!
!! make natural eigenbasis, on which the impurity energy matrix
!! should is diagonal
!!
  subroutine atomic_build_natural()
     use constants, only : dp
     use constants, only : czero
     use constants, only : mystd

     use control, only : ibasis
     use control, only : icu, icf, isoc
     use control, only : norbs

     use m_spmat, only : umat
     use m_spmat, only : tmat

     implicit none

!! local variables
     ! transformation matrix from real orbital basis to
     ! complex orbital basis
     complex(dp) :: tmat_r2c(norbs,norbs)

     ! transformation matrix from complex orbital basis
     ! to real orbital basis
     complex(dp) :: tmat_c2r(norbs,norbs)

     ! dummy Coulomb interaction matrix
     complex(dp) :: umat_tmp(norbs,norbs,norbs,norbs)

!! [body

     ! initialize them
     umat_tmp = czero
     tmat_r2c = czero
     tmat_c2r = czero

     ! make transformation matrix (tmat). it is used to transfer matrix
     ! emat from original basis to natural eigenbasis. the onsite energy
     ! of impurity (emat) should be updated as well
     !
     ! A: make tmat internally for different cases
     if ( ibasis == 1 ) then

         write(mystd,'(4X,a)') 'make transformation matrix internally'

         ! no spin-orbit coupling
         ! no crystal field splitting or it is diagonal
         !
         ! the original basis is the real orbital basis
         ! the natural eigenbasis is the real orbital basis
         if      ( isoc == 0 .and. icf <  2 ) then
             call atomic_natural_basis1()

         ! no spin-orbit coupling
         ! non-diagonal crystal field splitting
         !
         ! the original basis is the real orbital basis
         ! the natural eigenbasis is linear combination of real orbitals
         else if ( isoc == 0 .and. icf == 2 ) then
             call atomic_natural_basis2()

         ! with spin-orbit coupling
         ! no crystal field splitting
         !
         ! the original basis is the complex orbital basis
         ! the natural eigenbasis is the j^2 - j_z diagonal basis
         else if ( isoc == 1 .and. icf == 0 ) then
             call atomic_natural_basis3()

         ! with spin-orbit coupling
         ! with crystal field splitting
         !
         ! the original basis is the complex orbital basis
         ! the natural eigenbasis is linear combination of complex orbitals
         else if ( isoc == 1 .and. icf >  0 ) then
             call atomic_natural_basis4()

         endif ! back if      ( isoc == 0 .and. icf <  2 ) block

     ! B: read the transformation matrix (tmat)
     else

         write(mystd,'(4X,a)') 'read transformation matrix'

         ! note that emat is already built in atomic_build_spmat(),
         ! here we just read tmat
         call atomic_input_tmat()

     endif ! back if ( ibasis == 1 ) block

     ! dump emat for reference
     call atomic_dump_emat()

     ! dump tmat for reference
     call atomic_dump_tmat()

     ! we need transform Coulomb interaction U
     !
     ! for non-SOC case, the transformation matrix is defined as from
     ! real orbital basis to natural eigenbasis. so we have to make sure
     ! the Coulomb interaction U is at real orbital basis
     if ( isoc == 0 ) then

         write(mystd,'(4X,a)') 'transform Coulomb interaction to &
             &real orbital basis'

         ! for Slater-Cordon parameterized Coulomb interaction U,
         ! since it is defined at complex orbital basis, we first
         ! need to transfrom umat from complex orbital basis to real
         ! orbital basis
         !
         ! for Kanamori parameterized Coulomb interaction U, since
         ! it is already defined at real orbital basis, we do not
         ! need to transform it now
         if ( icu == 2 ) then
             call atomic_make_tmat_c2r(tmat_c2r)
             call atomic_tran_umat(tmat_c2r, umat, umat_tmp)
             umat = umat_tmp
         endif ! back if ( icu == 2 ) block

     ! for SOC case, the transformation matrix is defined as from complex
     ! orbital basis to natural eigenbasis. so we have to make sure the
     ! Coulomb interaction U is at complex orbital basis
     else

         write(mystd,'(4X,a)') 'transform Coulomb interaction to &
             &complex orbital basis'

         ! for Slater-Cordon parameterized Coulomb interaction U,
         ! since it is already defined at complex orbital basis, we
         ! do not need to transform it now
         !
         ! for Kanamori parameterized Coulomb interaction U, since
         ! it is defined at real orbital basis, so we have to transfrom
         ! umat from real orbital basis to complex orbital basis
         if ( icu == 1 .or. icu == 3 ) then
             call atomic_make_tmat_r2c(tmat_r2c)
             call atomic_tran_umat(tmat_r2c, umat, umat_tmp)
             umat = umat_tmp
         endif ! back if ( icu == 1 .or. icu == 3 ) block

     endif ! back if ( isoc == 0 ) block

     ! finally, transform umat from original basis to natural eigenbasis.
     ! the transformation matrix is just tmat
     write(mystd,'(4X,a)') 'transform Coulomb interaction to &
             &natural eigenbasis'
     !
     call atomic_tran_umat(tmat, umat, umat_tmp)
     umat = umat_tmp

     ! write the U matrix as reference
     call atomic_dump_umat()

!! body]

     !STOP

     return
  end subroutine atomic_build_natural

!!========================================================================
!!>>> manage memory for atomic eigenvalue problem solver               <<<
!!========================================================================

!!
!! @sub atomic_alloc_array
!!
!! allocate memory for global variables and then initialize them
!!
  subroutine atomic_alloc_array()
     use m_fock, only : cat_alloc_fock_basis
     use m_spmat, only : cat_alloc_spmat

     implicit none

!! [body

     ! allocate memory for Fock basis
     call cat_alloc_fock_basis()

     ! allocate memory for single particle matrices
     call cat_alloc_spmat()

!! body]

     return
  end subroutine atomic_alloc_array

!!
!! @sub atomic_final_array
!!
!! garbage collection for this code, please refer to atomic_alloc_array
!!
  subroutine atomic_final_array()
     use m_fock, only : cat_free_fock_basis
     use m_spmat, only : cat_free_spmat

     implicit none

!! [body

     ! deallocate memory for single particle matrices
     call cat_free_spmat()

     ! deallocate memory for Fock basis
     call cat_free_fock_basis()

!! body]

     return
  end subroutine atomic_final_array
