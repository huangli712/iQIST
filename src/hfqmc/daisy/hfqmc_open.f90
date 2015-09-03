

!!>>> cat_solver_id: return the solver identity
  subroutine cat_solver_id(I_solver_id)
     use dapi, only : solver_id_daisy

     implicit none

! external arguments
! solver identity
     integer, intent(out) :: I_solver_id

     I_solver_id = solver_id_daisy

     return
  end subroutine cat_solver_id

!!>>> cat_solver_status: return the solver status
  subroutine cat_solver_status(I_solver_status)
     use dapi, only : solver_is_ready_daisy

     implicit none

! external arguments
! solver status
     integer, intent(out) :: I_solver_status

     I_solver_status = solver_is_ready_daisy
     if ( I_solver_status == 0 ) then
         call s_print_error('cat_solver_status','sorry, the current solver is not ready!')
     endif ! back if ( I_solver_status == 0 ) block

     return
  end subroutine cat_solver_status

# if !defined (MPY)

!!>>> cat_init_hfqmc: initialize the hfqmc quantum impurity solver
!!>>> fortran version
  subroutine cat_init_hfqmc(I_mpi, I_solver)
     use dapi, only : T_mpi, T_daisy

     use control ! ALL

     implicit none

! external arguments
! type structure of mpi
     type (T_mpi), intent(in) :: I_mpi

! type structure of generic solver
     type (T_daisy), intent(in) :: I_solver

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
     niter  = I_solver%niter
     mstep  = I_solver%mstep
     mfreq  = I_solver%mfreq
     nsing  = I_solver%nsing
     ntime  = I_solver%ntime
     ntherm = I_solver%ntherm
     nsweep = I_solver%nsweep
     nclean = I_solver%nclean
     ncarlo = I_solver%ncarlo

! setup I_solver: real parameters
     Uc     = I_solver%Uc
     Jz     = I_solver%Jz
     mune   = I_solver%mune
     beta   = I_solver%beta
     part   = I_solver%part
     alpha  = I_solver%alpha

! print the running header for Hirsch-Fye quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call hfqmc_print_header()
     endif ! back if ( myid == master ) block

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call hfqmc_print_summary()
     endif ! back if ( myid == master ) block

! allocate memory and initialize
     call hfqmc_setup_array()

! prepare initial bath weiss's function, init self-consistent iteration
     call hfqmc_selfer_init()

     return
  end subroutine cat_init_hfqmc

# else   /* MPY */

!!>>> cat_init_hfqmc: initialize the hfqmc quantum impurity solver
!!>>> python version
  subroutine cat_init_hfqmc(my_id, num_procs)
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

! print the running header for Hirsch-Fye quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call hfqmc_print_header()
     endif ! back if ( myid == master ) block

! setup the important parameters for Hirsch-Fye quantum Monte Carlo
! quantum impurity solver and dynamical mean field theory self-consistent
! engine
     call hfqmc_config()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call hfqmc_print_summary()
     endif ! back if ( myid == master ) block

! allocate memory and initialize
     call hfqmc_setup_array()

! prepare initial bath weiss's function, init self-consistent iteration
     call hfqmc_selfer_init()

     return
  end subroutine cat_init_hfqmc

# endif  /* MPY */

!!>>> cat_exec_hfqmc: execute the hfqmc quantum impurity solver
  subroutine cat_exec_hfqmc(iter)
     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! call the Hirsch-Fye quantum Monte Carlo quantum impurity solver, to
! build the impurity green's function and self-energy function
     call hfqmc_impurity_solver(iter)

     return
  end subroutine cat_exec_hfqmc

!!>>> cat_stop_hfqmc: stop the hfqmc quantum impurity solver
  subroutine cat_stop_hfqmc()
     use control, only : myid, master

     implicit none

! deallocate memory and finalize
     call hfqmc_final_array()

! print the footer for Hirsch-Fye quantum Monte Carlo quantum impurity
! solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call hfqmc_print_footer()
     endif ! back if ( myid == master ) block

     return
  end subroutine cat_stop_hfqmc

!!>>> cat_set_wssf: setup the bath weiss's function
  subroutine cat_set_wssf(size_t, wssf_t)
     use constants, only : dp

     use control, only : norbs
     use control, only : mfreq
     use context, only : wssf

     implicit none

! external arguments
! size of wssf
     integer, intent(in)     :: size_t

! bath weiss's function
     complex(dp), intent(in) :: wssf_t(size_t)

! check whether size_t is correct
     if ( size_t /= size(wssf) ) then
         call s_print_error('cat_set_wssf','wrong dimension size of wssf_t')
     endif ! back if ( size_t /= size(wssf) ) block

! copy data
     wssf = reshape(wssf_t,(/mfreq,norbs/))

     return
  end subroutine cat_set_wssf

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
!!>>> note: the daisy code does not support this function now
  subroutine cat_set_ktau(size_t, ktau_t, ptau_t)
     use constants, only : dp

     implicit none

! external arguments
! size of ktau
     integer, intent(in)  :: size_t

! screening function K(\tau)
     real(dp), intent(in) :: ktau_t(size_t)

! first derivate of screening function K'(\tau)
     real(dp), intent(in) :: ptau_t(size_t)

! to avoid the warning from compiler
     call s_assert( size(ktau_t) == size_t )
     call s_assert( size(ptau_t) == size_t )
     call s_print_error('cat_set_ktau','sorry, this feature is not supported')

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
     grnf_t = reshape(grnf, (/mfreq*norbs/))

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
     sigf_t = reshape(sig2, (/mfreq*norbs/))

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
