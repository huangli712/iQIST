!!!=========+=========+=========+=========+=========+=========+=========+!
!!! NARCISSUS @ iQIST                                                    !
!!!                                                                      !
!!! A test program for dynamical mean field theory (DMFT) self-consistent!
!!! engine plus hybridization expansion version continuous time quantum  !
!!! Monte Carlo (CTQMC) quantum impurity solver                          !
!!! author  : Li Huang (at IOP/CAS & SPCLab/CAEP & UNIFR)                !
!!! version : v2015.01.06T                                               !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : this impurity solver is based on segment picture formalism !
!!!           dynamical screening effect is included                     !
!!!           any question, please contact with huangli712@gmail.com     !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!
!!
!! WARNING
!! =======
!!
!! If you want to obtain an executable program, please go to src/build/,
!! open make.sys and comment out the API flag. On the contrary, if you
!! want to compile narcissus as a library, please activate the API flag.
!!
!! Introduction
!! ============
!!
!! The narcissus code is a hybridization expansion version continuous time
!! quantum Monte Carlo quantum impurity solver. It adopts the segment
!! picture, and implements many useful features, such as the orthogonal
!! polynomial representation and the measurement of two-particle Green's
!! function, etc. The narcissus code can be also used to solve the Hubbard
!! -Holstein model and study the dynamical screening effect. It can be
!! considered as an impurity solver for extended dynamical mean field
!! theory as well. So it is a bit less efficient than the azalea code. And
!! it can be used as a standard to benchmark the other ctqmc impurity
!! solvers. The narcissus code also includes a mini dynamical mean field
!! theory engine which implements the self-consistent equation for Bethe
!! lattice in paramagnetic state. So you can use it to perform dynamical
!! mean field theory calculations quickly. Enjoy it.
!!
!! Usage
!! =====
!!
!! # ./ctqmc or bin/narcissus.x
!!
!! Input
!! =====
!!
!! solver.ctqmc.in (optional)
!! solver.eimp.in (optional)
!! solver.hyb.in (optional)
!! solver.ktau.in (optional)
!!
!! Output
!! ======
!!
!! terminal output
!! solver.green.bin.*
!! solver.green.dat
!! solver.grn.dat
!! solver.hybri.dat
!! solver.hyb.dat
!! solver.wss.dat
!! solver.sgm.dat
!! solver.hub.dat
!! solver.hist.dat
!! solver.prob.dat
!! solver.nmat.dat
!! solver.schi.dat
!! solver.ochi.dat
!! solver.twop.dat
!! solver.vrtx.dat
!! solver.pair.dat
!! solver.status.dat
!! solver.kernel.dat
!! etc.
!!
!! Running mode
!! ============
!!
!! case 1: isscf == 1 .and. isbin == 1
!! -----------------------------------
!!
!! call ctqmc_impurity_solver only, normal mode
!!
!! case 2: isscf == 1 .and. isbin == 2
!! -----------------------------------
!!
!! call ctqmc_impurity_solver only, binner mode
!!
!! case 3: isscf == 2 .and. isbin == 1
!! -----------------------------------
!!
!! call ctqmc_impurity_solver, normal mode
!! plus
!! call ctqmc_dmft_selfer
!! until convergence
!!
!! case 4: isscf == 2 .and. isbin == 2
!! -----------------------------------
!!
!! call ctqmc_impurity_solver, normal mode
!! plus
!! call ctqmc_dmft_selfer
!! until convergence
!! plus
!! call ctqmc_impurity_solver, binner mode
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!

# if !defined (API)

  program ctqmc_main
     use constants, only : mystd
     use mmpi, only : mp_init, mp_finalize
     use mmpi, only : mp_comm_rank, mp_comm_size
     use mmpi, only : mp_barrier

     use control, only : isscf, isbin
     use control, only : niter
     use control, only : nprocs, myid, master

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
     endif ! back if ( myid == master ) block

! setup the important parameters for continuous time quantum Monte Carlo
! quantum impurity solver and dynamical mean field theory self-consistent
! engine
     call ctqmc_config()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_summary()
     endif ! back if ( myid == master ) block

! allocate memory and initialize
     call ctqmc_setup_array()

! prepare initial hybridization function, init self-consistent iteration
     call ctqmc_selfer_init()

!!========================================================================
!!>>> DMFT ITERATION BEGIN                                             <<<
!!========================================================================

! case A: one-shot non-self-consistent mode
!-------------------------------------------------------------------------
! it is suitable for local density approximation plus dynamical mean field
! theory calculation
     if ( isscf == 1 .and. isbin == 1 ) then

! set the iter number
         iter = niter

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a,i3,a)') 'NARCISSUS >>> DMFT iter:', iter, ' <<< SELFING'
         endif ! back if ( myid == master ) block

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
         endif ! back if ( isscf == 1 ) block

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a,i3,a)') 'NARCISSUS >>> DMFT iter:', iter, ' <<< SELFING'
         endif ! back if ( myid == master ) block

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
         endif ! back if ( convergence .eqv. .true. ) block

     enddo DMFT_CTQMC_ITERATION ! over iter={1,niter} loop

! case C: binner mode
!-------------------------------------------------------------------------
! perform quantum Monte Carlo data binning
     if ( isbin == 2 ) then

! set the iter number
         iter = 999

! write the iter to screen
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,a,i3,a)') 'NARCISSUS >>> DMFT iter:', iter, ' <<< BINNING'
         endif ! back if ( myid == master ) block

! accumulate the quantum Monte Carlo data
         call ctqmc_impurity_solver(iter)

     endif ! back if ( isbin == 2 ) block

!!========================================================================
!!>>> DMFT ITERATION END                                               <<<
!!========================================================================

! deallocate memory and finalize
     call ctqmc_final_array()

! print the footer for continuous time quantum Monte Carlo quantum impurity
! solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_footer()
     endif ! back if ( myid == master ) block

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

  end program ctqmc_main

# endif  /* API */




!!>>> cat_solver_id: return the solver identity
  subroutine cat_solver_id(I_solver_id)
     use api, only : solver_id_narcissus

     implicit none

! external arguments
! solver identity
     integer, intent(out) :: I_solver_id

     I_solver_id = solver_id_narcissus

     return
  end subroutine cat_solver_id

!!>>> cat_solver_status: return the solver status
  subroutine cat_solver_status(I_solver_status)
     use api, only : solver_is_ready_narcissus

     implicit none

! external arguments
! solver status
     integer, intent(out) :: I_solver_status

     I_solver_status = solver_is_ready_narcissus
     if ( I_solver_status == 0 ) then
         call s_print_error('cat_solver_status','sorry, the current solver is not ready!')
     endif ! back if ( I_solver_status == 0 ) block

     return
  end subroutine cat_solver_status

# if !defined (F2PY)

!!>>> cat_init_ctqmc: initialize the ctqmc quantum impurity solver
!!>>> fortran version
  subroutine cat_init_ctqmc(I_mpi, I_solver)
     use api, only : T_mpi, T_segment_narcissus

     use control ! ALL

     implicit none

! external arguments
! type structure of mpi
     type (T_mpi), intent(in) :: I_mpi

! type structure of generic solver
     type (T_segment_narcissus), intent(in) :: I_solver

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
     isort  = I_solver%isort
     isvrt  = I_solver%isvrt
     isscr  = I_solver%isscr
     nband  = I_solver%nband
     nspin  = I_solver%nspin
     norbs  = I_solver%norbs
     ncfgs  = I_solver%ncfgs
     niter  = I_solver%niter
     lemax  = I_solver%lemax
     legrd  = I_solver%legrd
     chmax  = I_solver%chmax
     chgrd  = I_solver%chgrd
     mkink  = I_solver%mkink
     mfreq  = I_solver%mfreq
     nffrq  = I_solver%nffrq
     nbfrq  = I_solver%nbfrq
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
     lc     = I_solver%lc
     wc     = I_solver%wc
     mune   = I_solver%mune
     beta   = I_solver%beta
     part   = I_solver%part
     alpha  = I_solver%alpha

! print the running header for continuous time quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_header()
     endif ! back if ( myid == master ) block

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_summary()
     endif ! back if ( myid == master ) block

! allocate memory and initialize
     call ctqmc_setup_array()

! prepare initial hybridization function, init self-consistent iteration
     call ctqmc_selfer_init()

     return
  end subroutine cat_init_ctqmc

# else   /* F2PY */

!!>>> cat_init_ctqmc: initialize the ctqmc quantum impurity solver
!!>>> python version
  subroutine cat_init_ctqmc(my_id, num_procs)
     use control, only : nprocs, myid, master

     implicit none

! external arguments
! id for current process
     integer, intent(in) :: my_id

! number of processors
     integer, intent(in) :: num_procs

! initialize mpi envirnoment
     myid = my_id
     nprocs = num_procs

! print the running header for continuous time quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_header()
     endif ! back if ( myid == master ) block

! setup the important parameters for continuous time quantum Monte Carlo
! quantum impurity solver and dynamical mean field theory self-consistent
! engine
     call ctqmc_config()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_summary()
     endif ! back if ( myid == master ) block

! allocate memory and initialize
     call ctqmc_setup_array()

! prepare initial hybridization function, init self-consistent iteration
     call ctqmc_selfer_init()

     return
  end subroutine cat_init_ctqmc

# endif  /* F2PY */

!!>>> cat_exec_ctqmc: execute the ctqmc quantum impurity solver
  subroutine cat_exec_ctqmc(iter)
     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! call the continuous time quantum Monte Carlo quantum impurity solver, to
! build the impurity green's function and self-energy function
     call ctqmc_impurity_solver(iter)

     return
  end subroutine cat_exec_ctqmc

!!>>> cat_stop_ctqmc: stop the ctqmc quantum impurity solver
  subroutine cat_stop_ctqmc()
     use control, only : myid, master

     implicit none

! deallocate memory and finalize
     call ctqmc_final_array()

! print the footer for continuous time quantum Monte Carlo quantum impurity
! solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_footer()
     endif ! back if ( myid == master ) block

     return
  end subroutine cat_stop_ctqmc

!!>>> cat_set_hybf: setup the hybridization function
  subroutine cat_set_hybf(size_t, hybf_t)
     use constants, only : dp

     use control, only : norbs
     use control, only : mfreq
     use context, only : hybf

     implicit none

! external arguments
! size of hybf
     integer, intent(in)     :: size_t

! hybridization function
     complex(dp), intent(in) :: hybf_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(hybf) ) then
         call s_print_error('cat_set_hybf','wrong dimension size of hybf_t')
     endif ! back if ( size_t /= size(hybf) ) block

! copy data
     hybf = reshape(hybf_t,(/mfreq,norbs,norbs/))

     return
  end subroutine cat_set_hybf

!!>>> cat_set_symm: setup the symmetry vector
  subroutine cat_set_symm(size_t, symm_t)
     use context, only : symm

     implicit none

! external arguments
! size of symm
     integer, intent(in) :: size_t

! symmetry vector
     integer, intent(in) :: symm_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(symm) ) then
         call s_print_error('cat_set_symm','wrong dimension size of symm_t')
     endif ! back if ( size_t /= size(symm) ) block

! copy data
     symm = symm_t

     return
  end subroutine cat_set_symm

!!>>> cat_set_eimp: setup the impurity level
  subroutine cat_set_eimp(size_t, eimp_t)
     use constants, only : dp

     use context, only : eimp

     implicit none

! external arguments
! size of eimp
     integer, intent(in)  :: size_t

! impurity level
     real(dp), intent(in) :: eimp_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(eimp) ) then
         call s_print_error('cat_set_eimp','wrong dimension size of eimp_t')
     endif ! back if ( size_t /= size(eimp) ) block

! copy data
     eimp = eimp_t

     return
  end subroutine cat_set_eimp

!!>>> cat_set_ktau: setup the screening function and its first derivates
  subroutine cat_set_ktau(size_t, ktau_t, ptau_t)
     use constants, only : dp

     use context, only : ktau, ptau

     implicit none

! external arguments
! size of ktau
     integer, intent(in)  :: size_t

! screening function K(\tau)
     real(dp), intent(in) :: ktau_t(size_t)

! first derivate of screening function K'(\tau)
     real(dp), intent(in) :: ptau_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(ktau) ) then
         call s_print_error('cat_set_ktau','wrong dimension size of ktau_t')
     endif ! back if ( size_t /= size(ktau) ) block

     if ( size_t /= size(ptau) ) then
         call s_print_error('cat_set_ktau','wrong dimension size of ptau_t')
     endif ! back if ( size_t /= size(ptau) ) block

! copy data
     ktau = ktau_t
     ptau = ptau_t

     return
  end subroutine cat_set_ktau

!!>>> cat_get_grnf: extract the impurity green's function
  subroutine cat_get_grnf(size_t, grnf_t)
     use constants, only : dp

     use control, only : norbs
     use control, only : mfreq
     use context, only : grnf

     implicit none

! external arguments
! size of grnf
     integer, intent(in)      :: size_t

! impurity green's function
     complex(dp), intent(out) :: grnf_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(grnf) ) then
         call s_print_error('cat_get_grnf','wrong dimension size of grnf_t')
     endif ! back if ( size_t /= size(grnf) ) block

! copy data
     grnf_t = reshape(grnf, (/mfreq*norbs*norbs/))

     return
  end subroutine cat_get_grnf

!!>>> cat_get_sigf: extract the self-energy function
  subroutine cat_get_sigf(size_t, sigf_t)
     use constants, only : dp

     use control, only : norbs
     use control, only : mfreq
     use context, only : sig2

     implicit none

! external arguments
! size of sigf
     integer, intent(in)      :: size_t

! self-energy function
     complex(dp), intent(out) :: sigf_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(sig2) ) then
         call s_print_error('cat_get_sigf','wrong dimension size of sigf_t')
     endif ! back if ( size_t /= size(sig2) ) block

! copy data
     sigf_t = reshape(sig2, (/mfreq*norbs*norbs/))

     return
  end subroutine cat_get_sigf

!!>>> cat_get_nmat: extract the occupation number
  subroutine cat_get_nmat(size_t, nmat_t)
     use constants, only : dp

     use control, only : norbs
     use context, only : nmat

     implicit none

! external arguments
! size of nmat
     integer, intent(in)   :: size_t

! occupation number
     real(dp), intent(out) :: nmat_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(nmat) ) then
         call s_print_error('cat_get_nmat','wrong dimension size of nmat_t')
     endif ! back if ( size_t /= size(nmat) ) block

! copy data
     nmat_t = reshape(nmat, (/norbs/))

     return
  end subroutine cat_get_nmat

!!>>> cat_get_nnmat: extract the double occupation number
  subroutine cat_get_nnmat(size_t, nnmat_t)
     use constants, only : dp

     use control, only : norbs
     use context, only : nnmat

     implicit none

! external arguments
! size of nnmat
     integer, intent(in)   :: size_t

! double occupation number
     real(dp), intent(out) :: nnmat_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(nnmat) ) then
         call s_print_error('cat_get_nnmat','wrong dimension size of nnmat_t')
     endif ! back if ( size_t /= size(nnmat) ) block

! copy data
     nnmat_t = reshape(nnmat, (/norbs*norbs/))

     return
  end subroutine cat_get_nnmat
