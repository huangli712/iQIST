!=========+=========+=========+=========+=========+=========+=========+>>>
! A test program for dynamical mean field theory (DMFT) self-consistent  !
! engine plus hybridization expansion version continuous time quantum    !
! Monte Carlo (CTQMC) quantum impurity solver                            !
! author  : li huang                                                     !
! version : v2014.01.13T                                                 !
! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK               !
! comment : this impurity solver is based on segment picture formalism   !
!           any question, please contact with huangli712@gmail.com       !
!=========+=========+=========+=========+=========+=========+=========+>>>

# if !defined (API)

  program ctqmc_main
     use constants
     use control

     use mmpi

     implicit none

! local variables
! loop index
     integer :: iter

! convergence flag
     logical :: convergence

! initialize mpi envirnoment
# if defined (MPI)

! initialize the mpi execution environment
     call mp_init()

! determines the rank of the calling process in the communicator
     call mp_comm_rank(myid)

! determines the size of the group associated with a communicator
     call mp_comm_size(nprocs)

# endif  /* MPI */

! print the running header for continuous time quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_header()
     endif

! setup the important parameters for continuous time quantum Monte Carlo
! quantum impurity solver and dynamical mean field theory self-consistent
! engine
     call ctqmc_config()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_summary()
     endif

! allocate memory and initialize
     call ctqmc_setup_array()

! prepare initial hybridization function, init self-consistent iteration
     call ctqmc_selfer_init()

!-------------------------------------------------------------------------
! note: running mode                                                     !
!-------------------------------------------------------------------------
!    if isscf == 1 .and. isbin == 1                                      !
!        call ctqmc_impurity_solver only, normal mode                    !
!                                                                        !
!    if isscf == 1 .and. isbin == 2                                      !
!        call ctqmc_impurity_solver only, binner mode                    !
!                                                                        !
!    if isscf == 2 .and. isbin == 1                                      !
!        call ctqmc_impurity_solver, normal mode                         !
!        plus                                                            !
!        call ctqmc_dmft_selfer                                          !
!        until convergence                                               !
!                                                                        !
!    if isscf == 2 .and. isbin == 2                                      !
!        call ctqmc_impurity_solver, normal mode                         !
!        plus                                                            !
!        call ctqmc_dmft_selfer                                          !
!        until convergence                                               !
!        plus                                                            !
!        call ctqmc_impurity_solver, binner mode                         !
!-------------------------------------------------------------------------

!=========================================================================
!>>> DMFT ITERATION BEGIN                                              <<<
!=========================================================================

! case A: one-shot non-self-consistent mode
!-------------------------------------------------------------------------
! it is suitable for local density approximation plus dynamical mean field
! theory calculation
     if ( isscf == 1 .and. isbin == 1 ) then

! set the iter number
         iter = niter

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a,i3,a)') 'AZALEA >>> DMFT iter:', iter, ' <<< SELFING'
         endif

! call the continuous time quantum Monte Carlo quantum impurity solver, to
! build the impurity green's function and self-energy function
         call ctqmc_impurity_solver(iter)

     endif ! back if ( isscf == 1 .and. isbin == 1 ) block

! case B: self-consistent mode
!-------------------------------------------------------------------------
! it is suitable for lattice model hamiltonian plus dynamical mean field
! theory calculation
     DMFT_CTQMC_ITERATION: do iter=1,niter

! check the running mode
         if ( isscf == 1 ) then
             EXIT DMFT_CTQMC_ITERATION ! jump out the iteration
         endif

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a,i3,a)') 'AZALEA >>> DMFT iter:', iter, ' <<< SELFING'
         endif

! call the continuous time quantum Monte Carlo quantum impurity solver, to
! build the impurity green's function and self-energy function
         call ctqmc_impurity_solver(iter)

! call the self-consistent engine for dynamical mean field theory, to build
! the bath weiss's function and hybridization function
         call ctqmc_dmft_selfer()

! check convergence for dynamical mean field theory iteration
         convergence = .false.
         call ctqmc_dmft_conver(iter, convergence)

! now convergence is achieved
         if ( convergence .eqv. .true. ) then
             EXIT DMFT_CTQMC_ITERATION ! jump out the iteration
         endif

     enddo DMFT_CTQMC_ITERATION ! over iter={1,niter} loop

! case C: binner mode
!-------------------------------------------------------------------------
! perform quantum Monte Carlo data binning
     if ( isbin == 2 ) then

! set the iter number
         iter = 999

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a,i3,a)') 'AZALEA >>> DMFT iter:', iter, ' <<< BINNING'
         endif

! accumulate the quantum Monte Carlo data
         call ctqmc_impurity_solver(iter)

     endif ! back if ( isbin == 2 ) block

!=========================================================================
!>>> DMFT ITERATION END                                                <<<
!=========================================================================

! deallocate memory and finalize
     call ctqmc_final_array()

! print the footer for continuous time quantum Monte Carlo quantum impurity
! solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_footer()
     endif

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

  end program ctqmc_main

# endif  /* API */

!>>> initialize the ctqmc quantum impurity solver
  subroutine cat_init_ctqmc(I_mpi, I_solver)
     use api
     use control

     implicit none

! external arguments
! type structure of mpi
     type (T_mpi), intent(in) :: I_mpi

! type structure of generic solver
     type (T_segment_azalea), intent(in) :: I_solver

! setup I_mpi
     nprocs = I_mpi%nprocs
     myid   = I_mpi%myid
     master = I_mpi%master
     cid    = I_mpi%cid
     cx     = I_mpi%cx
     cy     = I_mpi%cy

! setup I_solver: integer parameters
     isscf  = I_solver%isscf
     issun  = I_solver%issun
     isspn  = I_solver%isspn
     isbin  = I_solver%isbin
     nband  = I_solver%nband
     nspin  = I_solver%nspin
     norbs  = I_solver%norbs
     ncfgs  = I_solver%ncfgs
     niter  = I_solver%niter
     mkink  = I_solver%mkink
     mfreq  = I_solver%mfreq
     nfreq  = I_solver%nfreq
     ntime  = I_solver%ntime
     nflip  = I_solver%nflip
     ntherm = I_solver%ntherm
     nsweep = I_solver%nsweep
     nwrite = I_solver%nwrite
     nclean = I_solver%nclean
     nmonte = I_solver%nmonte
     ncarlo = I_solver%ncarlo

! setup I_solver: real parameters
     U      = I_solver%U
     Uc     = I_solver%Uc
     Uv     = I_solver%Uv
     Jz     = I_solver%Jz
     Js     = I_solver%Js
     Jp     = I_solver%Jp
     mune   = I_solver%mune
     beta   = I_solver%beta
     part   = I_solver%part
     alpha  = I_solver%alpha

! print the running header for continuous time quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_header()
     endif

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_summary()
     endif

! allocate memory and initialize
     call ctqmc_setup_array()

! prepare initial hybridization function, init self-consistent iteration
     call ctqmc_selfer_init()

     return
  end subroutine cat_init_ctqmc

!>>> execute the ctqmc quantum impurity solver
  subroutine cat_exec_ctqmc(iter)
     use control

     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! call the continuous time quantum Monte Carlo quantum impurity solver, to
! build the impurity green's function and self-energy function
     call ctqmc_impurity_solver(iter)

     return
  end subroutine cat_exec_ctqmc

!>>> stop the ctqmc quantum impurity solver
  subroutine cat_stop_ctqmc()
     use control

     implicit none

! deallocate memory and finalize
     call ctqmc_final_array()

! print the footer for continuous time quantum Monte Carlo quantum impurity
! solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_footer()
     endif

     return
  end subroutine cat_stop_ctqmc

!>>> setup the hybridization function
  subroutine cat_set_hybf(size_t, hybf_t)
     use control
     use context

     implicit none

! external arguments
! size of hybf
     integer, intent(in) :: size_t

! hybridization function
     complex(dp), intent(in) :: hybf_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(hybf) ) then
         call ctqmc_print_error('cat_set_hybf', 'wrong dimension size of hybf_t')
     endif

! copy data
     hybf = reshape(hybf_t,(/mfreq,norbs,norbs/))

     return
  end subroutine cat_set_hybf

!>>> setup the symmetry vector
  subroutine cat_set_symm(size_t, symm_t)
     use control
     use context

     implicit none

! external arguments
! size of symm
     integer, intent(in) :: size_t

! symmetry vector
     integer, intent(in) :: symm_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(symm) ) then
         call ctqmc_print_error('cat_set_symm', 'wrong dimension size of symm_t')
     endif

! copy data
     symm = symm_t

     return
  end subroutine cat_set_symm

!>>> setup the impurity level
  subroutine cat_set_eimp(size_t, eimp_t)
     use control
     use context

     implicit none

! external arguments
! size of eimp
     integer, intent(in) :: size_t

! impurity level
     real(dp), intent(in) :: eimp_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(eimp) ) then
         call ctqmc_print_error('cat_set_eimp', 'wrong dimension size of eimp_t')
     endif

! copy data
     eimp = eimp_t

     return
  end subroutine cat_set_eimp

!>>> extract the impurity green's function
  subroutine cat_get_grnf(size_t, grnf_t)
     use control
     use context

     implicit none

! external arguments
! size of grnf
     integer, intent(in) :: size_t

! impurity green's function
     complex(dp), intent(out) :: grnf_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(grnf) ) then
         call ctqmc_print_error('cat_get_grnf', 'wrong dimension size of grnf_t')
     endif

! copy data
     grnf_t = reshape(grnf, (/mfreq*norbs*norbs/))

     return
  end subroutine cat_get_grnf

!>>> extract the self-energy function
  subroutine cat_get_sigf(size_t, sigf_t)
     use control
     use context

     implicit none

! external arguments
! size of sigf
     integer, intent(in) :: size_t

! self-energy function
     complex(dp), intent(out) :: sigf_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(sig2) ) then
         call ctqmc_print_error('cat_get_sigf', 'wrong dimension size of sigf_t')
     endif

! copy data
     sigf_t = reshape(sig2, (/mfreq*norbs*norbs/))

     return
  end subroutine cat_get_sigf
