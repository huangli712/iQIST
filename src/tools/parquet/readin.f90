     subroutine analyze_config()
     use constants
     use control

     implicit none

! local variables
! used to check whether the input file (solver.ctqmc.in) exists
     logical :: exists

!=========================================================================
! setup dynamical mean field theory self-consistent engine related common variables
!=========================================================================
     isscf  = 2            ! non-self-consistent (1) or self-consistent mode (2)
     issun  = 2            ! without symmetry    (1) or with symmetry   mode (2)
     isspn  = 1            ! spin projection, PM (1) or AFM             mode (2)
     isbin  = 2            ! without binning     (1) or with binning    mode (2)
     isort  = 1            ! normal measurement  (1) or legendre polynomial  (2) or chebyshev polynomial (3)
     isvrt  = 1            ! without vertex      (1) or with vertex function (2)
!-------------------------------------------------------------------------
     nband  = 1            ! number of correlated bands
     nspin  = 2            ! number of spin projection
     norbs  = nspin*nband  ! number of correlated orbitals (= nband * nspin)
     ncfgs  = 2**norbs     ! number of atomic states
     nzero  = 128          ! maximum number of non-zero elements in sparse matrix style
     niter  = 20           ! maximum number of DMFT + CTQMC self-consistent iterations
!-------------------------------------------------------------------------
     U      = 4.00_dp      ! U : average Coulomb interaction
     Uc     = 4.00_dp      ! Uc: intraorbital Coulomb interaction
     Uv     = 4.00_dp      ! Uv: interorbital Coulomb interaction, Uv = Uc - 2 * Jz for t2g system
     Jz     = 0.00_dp      ! Jz: Hund's exchange interaction in z axis (Jz = Js = Jp = J)
     Js     = 0.00_dp      ! Js: spin-flip term
     Jp     = 0.00_dp      ! Jp: pair-hopping term
!-------------------------------------------------------------------------
     mune   = 2.00_dp      ! chemical potential or fermi level
     beta   = 8.00_dp      ! inversion of temperature
     part   = 0.50_dp      ! coupling parameter t for Hubbard model
     alpha  = 0.70_dp      ! mixing parameter for self-consistent engine
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!=========================================================================
! setup continuous time quantum Monte Carlo quantum impurity solver related common variables
!=========================================================================
     lemax  = 32           ! maximum order for legendre polynomial
     legrd  = 20001        ! number of mesh points for legendre polynomial
     chmax  = 32           ! maximum order for chebyshev polynomial
     chgrd  = 20001        ! number of mesh points for chebyshev polynomial
!-------------------------------------------------------------------------
     mkink  = 1024         ! maximum perturbation expansions order
     mfreq  = 8193         ! maximum number of matsubara frequency
!-------------------------------------------------------------------------
     nffrq  = 32           ! number of matsubara frequency for the two-particle green's function
     nbfrq  = 8            ! number of bosonic frequency for the two-particle green's function
     nfreq  = 128          ! maximum number of matsubara frequency sampling by quantum impurity solver
     ntime  = 1024         ! number of time slice
     npart  = 16           ! number of parts that the imaginary time axis is split
     nflip  = 20000        ! flip period for spin up and spin down states
     ntherm = 200000       ! maximum number of thermalization steps
     nsweep = 20000000     ! maximum number of quantum Monte Carlo sampling steps
     nwrite = 2000000      ! output period
     nclean = 100000       ! clean update period
     nmonte = 10           ! how often to sampling the gmat and nmat
     ncarlo = 10           ! how often to sampling the gtau and prob
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! inquire file status: solver.ctqmc.in
         inquire (file = 'solver.ctqmc.in', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
             open(mytmp, file='solver.ctqmc.in', form='formatted', status='unknown')

             read(mytmp,*)
             read(mytmp,*)
             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) isscf                                         !
             read(mytmp,*) issun                                         !
             read(mytmp,*) isspn                                         !
             read(mytmp,*) isbin                                         !
             read(mytmp,*) isort                                         !
             read(mytmp,*) isvrt                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) nband                                         !
             read(mytmp,*) nspin                                         !
             read(mytmp,*) norbs                                         !
             read(mytmp,*) ncfgs                                         !
             read(mytmp,*) nzero                                         !
             read(mytmp,*) niter                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) U                                             !
             read(mytmp,*) Uc                                            !
             read(mytmp,*) Uv                                            !
             read(mytmp,*) Jz                                            !
             read(mytmp,*) Js                                            !
             read(mytmp,*) Jp                                            !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) mune                                          !
             read(mytmp,*) beta                                          !
             read(mytmp,*) part                                          !
             read(mytmp,*) alpha                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) lemax                                         !
             read(mytmp,*) legrd                                         !
             read(mytmp,*) chmax                                         !
             read(mytmp,*) chgrd                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) mkink                                         !
             read(mytmp,*) mfreq                                         !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             read(mytmp,*)
!------------------------------------------------------------------------+
             read(mytmp,*) nffrq                                         !
             read(mytmp,*) nbfrq                                         !
             read(mytmp,*) nfreq                                         !
             read(mytmp,*) ntime                                         !
             read(mytmp,*) npart                                         !
             read(mytmp,*) nflip                                         !
             read(mytmp,*) ntherm                                        !
             read(mytmp,*) nsweep                                        !
             read(mytmp,*) nwrite                                        !
             read(mytmp,*) nclean                                        !
             read(mytmp,*) nmonte                                        !
             read(mytmp,*) ncarlo                                        !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

             close(mytmp)
         endif ! back if ( exists .eqv. .true. ) block

     return
     end subroutine analyze_config

     subroutine analyze_parameterset()
     use constants
     use control
     implicit none
! local variables
! used to check whether the input file (parameterset) exists
     logical :: exists

! temporary real dummy variable
     real(dp) :: rtmp1

! inquire file status: parameterset
         inquire (file = 'parameterset', exist = exists)

! read in parameters, default setting should be overrided
         if ( exists .eqv. .true. ) then
             open(mytmp, file='parameterset', form='formatted', status='unknown')
!------------------------------------------------------------------------+
             read(mytmp,*) rtmp1                                         !
             read(mytmp,*) rtmp1                                         !
             read(mytmp,*) rtmp1                                         !
             read(mytmp,*) rtmp1                                         !
             read(mytmp,*) rtmp1                                         !
             read(mytmp,*) rtmp1                                         !
             read(mytmp,*) lambda                                        !
             read(mytmp,*) t1                                            !
             read(mytmp,*) rtmp1                                         !
             read(mytmp,*) rtmp1                                         !
             read(mytmp,*) mune                                          !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+

         write(*,*) "Converged chemical potential is ",mune

         endif ! back if ( exists .eqv. .true. ) block

!------------------------------------------------------------------------+
!     Redefine parameters for initial test
         ktime=2
!        Define nt_2p here
         nt_2p = nffrq*4*ktime*ktime
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^+


     return
     end subroutine analyze_parameterset


!>>> allocate memory for global variables and then initialize them
     subroutine analyze_setup_array()
     use context
     implicit none

! allocate memory for context module
     call analyze_allocate_memory_gmat()

     call analyze_allocate_memory_two_mat_local()

     return
     end subroutine analyze_setup_array
     
!>>> allocate memory for global variables and then initialize them: two_mat_non_local
     subroutine analyze_setup_array_two_mat_non_local()
     use context
     implicit none

! allocate memory for context module

     call analyze_allocate_memory_two_mat_non_local()

     return
     end subroutine analyze_setup_array_two_mat_non_local


!>>> Read in single- and two- particle functions
     subroutine analyze_readin_array()

     implicit none
!    Read in one-ptle function in matsubara frequency space
     call analyze_read_grnf()
     call analyze_read_sigma()

!    Read in two-ptle green's function in matsubara frequency space
     call analyze_read_pp_singlet()
     call analyze_read_ph_charge()
     call analyze_read_ph_magnetic()

     return
     end subroutine analyze_readin_array


!>>> Read in impurity green's function in matsubara frequency space
     subroutine analyze_read_grnf()
     use constants
     use control
     use context

     implicit none

! external arguments
     real(dp) :: Regrnf1, Regrnf2, Imgrnf1, Imgrnf2

! temporary buffer
     character*10 buf1

! local variables
! loop index
     integer :: ischeck=1
     integer :: itmp1
     integer :: i
     integer :: j

! open data file: solver.grn.dat
     open(mytmp, file='solver.grn.dat', form='formatted', status='unknown')

!! read the first line
!!     read(mytmp,*) buf1

! read it
     do i=1,nband
         do j=1,mfreq!mfreq=8193
            
         if(j.le.nffrq/2+(nbfrq-1).and.i.eq.3) then !nffrq=22, nbfrq=6, nffrq/2+(nbfrq-1)=15
             read(mytmp,'(i5,5f16.8)') itmp1, rmesh(j), &
                                  Regrnf1, Imgrnf1, &
                                  Regrnf2, Imgrnf2

             grnf(j,i-nband+1)=dcmplx(Regrnf1,Imgrnf1)
             grnf(j,i-nband+2)=dcmplx(Regrnf2,Imgrnf2)
         else
             read(mytmp,*) buf1
         end if ! if {j <= nffrq}

         enddo ! over j={1,mfreq} loop
         read(mytmp,*) ! read empty lines
         read(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

! assign the negative frequencies
     do j=-nffrq/2+1-(nbfrq-1), 0
        rmesh(j)=-rmesh(-j+1)
	grnf(j,1) = dconjg(grnf(-j+1,1))
	grnf(j,2) = dconjg(grnf(-j+1,2))
     end do 

!    print out green's function for checking: just postive frequencies
     if(ischeck.eq.1) then
     open(mytmp, file='analyze.grn.dat', form='formatted', status='unknown')
     do i=1,1
         do j=-nffrq/2+1-(nbfrq-1),nffrq/2+nbfrq-1
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(grnf(j,i)), &
                                 aimag(grnf(j,i)), &
                      real(grnf(j,i+1)), &
                     aimag(grnf(j,i+1))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop
! close data file
     close(mytmp)
     end if ! is(ischeck.eq.1)

     return
     end subroutine analyze_read_grnf

!>>> Read in impurity self energy in matsubara frequency space
     subroutine analyze_read_sigma()
     use constants
     use control
     use context

     implicit none

! external arguments
     real(dp) :: Resigma1, Resigma2, Imsigma1, Imsigma2

! temporary buffer
     character*10 buf1

! local variables
! loop index
     integer :: ischeck=1
     integer :: itmp1
     integer :: i
     integer :: j

! open data file: solver.grn.dat
     open(mytmp, file='solver.sgm.dat', form='formatted', status='unknown')

!! read the first line
!!     read(mytmp,*) buf1

! read it
     do i=1,nband
         do j=1,mfreq!mfreq=8193
            
         if(j.le.nffrq/2+(nbfrq-1).and.i.eq.3) then !nffrq=22, nbfrq=6, nffrq/2+(nbfrq-1)=16
             read(mytmp,'(i5,5f16.8)') itmp1, rmesh(j), &
                                  Resigma1, Imsigma1, &
                                  Resigma2, Imsigma2

             sigf(j,i-nband+1)=dcmplx(Resigma1,Imsigma1)
             sigf(j,i-nband+2)=dcmplx(Resigma2,Imsigma2)
         else
             read(mytmp,*) buf1
         end if ! if {j <= nffrq}

         enddo ! over j={1,mfreq} loop
         read(mytmp,*) ! read empty lines
         read(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

! assign the negative frequencies
     do j=-nffrq/2+1-(nbfrq-1), 0
!        rmesh(j)=-rmesh(-j+1) ! rmesh is already defined in analyze_read_grnf
	sigf(j,1) = dconjg(sigf(-j+1,1))
	sigf(j,2) = dconjg(sigf(-j+1,2))
     end do      
     
!    print out self energy for checking: only the positive frequencies
     if(ischeck.eq.1) then
     open(mytmp, file='analyze.sgm.dat', form='formatted', status='unknown')
     do i=1,1
         do j=-nffrq/2+1-(nbfrq-1),nffrq/2+nbfrq-1
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(sigf(j,i)), &
                                 aimag(sigf(j,i)), &
                      real(sigf(j,i+1)), &
                     aimag(sigf(j,i+1))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop
! close data file
     close(mytmp)
     end if ! is(ischeck.eq.1)

     return
     end subroutine analyze_read_sigma

!>>> Read in two-ptle green's function in matsubara frequency space
!>>> Read in particle-particle singlet pairing function
     subroutine analyze_read_pp_singlet()
     use constants
     use control
     use context
     implicit none
! local variables

     real(dp) :: Re_pp_s, Im_pp_s
! loop index
     integer :: ischeck=1
     integer :: itmp1,itmp2
     integer :: i
     integer :: j
     integer :: k
! temporary buffer
     character*10 buf1

     open(mytmp, file='solver.twop.pp.pairing.dat', form='formatted', status='unknown')
        do k=0,nbfrq-1
	  read(mytmp,*) buf1
        do j=-nffrq/2+1,nffrq/2
        do i=-nffrq/2+1,nffrq/2
        read(mytmp,'(2i5,2f18.8)') itmp1,itmp2, Re_pp_s, Im_pp_s

        g2_pp_s(j,i,k)=dcmplx(Re_pp_s,Im_pp_s)
        if(k.ne.0) g2_pp_s(-j+1,-i+1,-k) = dconjg(g2_pp_s(j,i,k))
        end do ! over j={1,nffrq} loop
        end do ! over i={1,nffrq} loop
	 read(mytmp,*) ! read empty lines
!!         read(mytmp,*)
        end do ! over k={0,nbfrq-1} loop

        
!    print out twop_pp_singlet for checking:
     if(ischeck.eq.1) then
     open(mytmp, file='analyze.twop.pp.pairing.dat', form='formatted', status='unknown')
     do k=0,nbfrq-1
	write(mytmp,'(i5)') k
     do j=-nffrq/2+1,nffrq/2
         do i=-nffrq/2+1,nffrq/2
             write(mytmp,'(2i5,2f18.8)') 2*j-1, 2*i-1, &
                                  real(g2_pp_s(j,i,k)), &
                                 aimag(g2_pp_s(j,i,k))
         enddo ! over i={1,nffrq} loop
     enddo ! over j={1,nffrq} loop
            write(mytmp,*) ' '
     enddo ! over k={0,nbfrq-1} loop
! close data file
     close(mytmp)
     end if ! is(ischeck.eq.1)

     return
     end subroutine analyze_read_pp_singlet

!>>> Read in particle-hole charge function
     subroutine analyze_read_ph_charge()
     use constants
     use control
     use context
     implicit none
! local variables

     real(dp) :: Re_ph, Im_ph
! loop index
     integer :: ischeck=1
     integer :: itmp1,itmp2
     integer :: i
     integer :: j
     integer :: k
! temporary buffer
     character*10 buf1


     open(mytmp, file='solver.twop.ph.charge.dat', form='formatted', status='unknown')
!       only one bosonic frequency for p-p channel
        do k=0,nbfrq-1
         read(mytmp,*) buf1
        do j=-nffrq/2+1,nffrq/2
        do i=-nffrq/2+1,nffrq/2
         read(mytmp,'(2i5,2f18.8)') itmp1,itmp2, Re_ph, Im_ph

        Fc_ph_c(j,i,k)= dcmplx(Re_ph,Im_ph)
        if(k.ne.0) Fc_ph_c(-j+1,-i+1,-k)=dconjg(Fc_ph_c(j,i,k))
        end do ! over j={1,nffrq} loop
        end do ! over i={1,nffrq} loop
         read(mytmp,*) ! read empty lines
         read(mytmp,*)
        end do ! over k={0,nbfrq-1} loop


!    print out twop_pp_singlet for checking:
     if(ischeck.eq.1) then
     open(mytmp, file='analyze.twop.ph.charge.dat', form='formatted', status='unknown')
     do k=-nbfrq+1,nbfrq-1
	write(mytmp,'(i5)') k
       do j=-nffrq/2+1,nffrq/2
         do i=-nffrq/2+1,nffrq/2
             write(mytmp,'(2i5,2f18.8)') 2*j-1, 2*i-1, &
                                  real(Fc_ph_c(j,i,k)), &
                                 aimag(Fc_ph_c(j,i,k))
         enddo ! over i={1,nffrq} loop
       enddo ! over j={1,nffrq} loop
       write(mytmp,*) ' '
     enddo ! over k={0,nbfrq-1} loop
! close data file
     close(mytmp)
     end if ! is(ischeck.eq.1)

     return
     end subroutine analyze_read_ph_charge

!>>> Read in particle-hole magnetic function
     subroutine analyze_read_ph_magnetic()
     use constants
     use control
     use context
     implicit none
! local variables

     real(dp) :: Re_ph, Im_ph
! loop index
     integer :: ischeck=1
     integer :: itmp1,itmp2
     integer :: i
     integer :: j
     integer :: k
! temporary buffer
     character*10 buf1


     open(mytmp, file='solver.twop.ph.magnet.dat', form='formatted', status='unknown')
!       only one bosonic frequency for p-p channel
        do k=0,nbfrq-1
         read(mytmp,*) buf1
        do j=-nffrq/2+1,nffrq/2
        do i=-nffrq/2+1,nffrq/2
         read(mytmp,'(2i5,2f20.10)') itmp1,itmp2, Re_ph, Im_ph

        Fc_ph_m(j,i,k)= dcmplx(Re_ph,Im_ph)
        if(k.ne.0) Fc_ph_m(-j+1,-i+1,-k)=dconjg(Fc_ph_m(j,i,k))
        end do ! over j={1,nffrq} loop
        end do ! over i={1,nffrq} loop
         read(mytmp,*) ! read empty lines
         read(mytmp,*)
        end do ! over k={0,nbfrq-1} loop

!    print out twop_pp_singlet for checking:
     if(ischeck.eq.1) then
     open(mytmp, file='analyze.twop.ph.magnet.dat', form='formatted', status='unknown')
     do k=-nbfrq+1,nbfrq-1
	write(mytmp,'(i5)') k
       do j=-nffrq/2+1,nffrq/2
         do i=-nffrq/2+1,nffrq/2
             write(mytmp,'(2i5,2f20.10)') 2*j-1, 2*i-1, &
                                  real(Fc_ph_m(j,i,k)), &
                                 aimag(Fc_ph_m(j,i,k))
         enddo ! over i={1,nffrq} loop
       enddo ! over j={1,nffrq} loop
       write(mytmp,*) ' '
     enddo ! over k={0,nbfrq-1} loop
! close data file
     close(mytmp)
     end if ! is(ischeck.eq.1)

     return
     end subroutine analyze_read_ph_magnetic

!>>> garbage collection for this program, please refer to analyze_setup_array
  subroutine analyze_final_array_two_mat_local()
     use context
     implicit none
! deallocate memory for context module
     call analyze_deallocate_memory_two_mat_local()
     return
  end subroutine analyze_final_array_two_mat_local
  
  
!>>> garbage collection for this program, please refer to analyze_setup_array
  subroutine analyze_final_array()
     use context

     implicit none

! deallocate memory for context module

     call analyze_deallocate_memory_gmat()
     call analyze_deallocate_memory_two_mat_non_local()

     return
  end subroutine analyze_final_array
