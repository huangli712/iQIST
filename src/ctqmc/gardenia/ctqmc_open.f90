!!!-----------------------------------------------------------------------
!!! project : gardenia
!!! program : cat_solver_id
!!!           cat_solver_status <<<---
!!!           cat_init_ctqmc
!!!           cat_exec_ctqmc
!!!           cat_stop_ctqmc    <<<---
!!!           cat_set_hybf
!!!           cat_set_symm
!!!           cat_set_eimp
!!!           cat_set_ktau
!!!           cat_set_uumat     <<<---
!!!           cat_get_grnf
!!!           cat_get_sigf
!!!           cat_get_nmat
!!!           cat_get_nnmat     <<<---
!!! source  : ctqmc_open.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 08/12/2015 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : to provide necessary application programming interface for
!!!           the hybridization expansion version continuous time quantum
!!!           Monte Carlo (CTQMC) quantum impurity solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> status query subroutines                                         <<<
!!========================================================================

!!>>> cat_solver_id: return the solver identity
  subroutine cat_solver_id(I_solver_id)
     use capi, only : solver_id_gardenia

     implicit none

! external arguments
! solver identity
     integer, intent(out) :: I_solver_id

! declare f2py directives
!F2PY intent(out) I_solver_id

     I_solver_id = solver_id_gardenia

     return
  end subroutine cat_solver_id

!!>>> cat_solver_status: return the solver status
  subroutine cat_solver_status(I_solver_status)
     use capi, only : solver_is_ready_gardenia

     implicit none

! external arguments
! solver status
     integer, intent(out) :: I_solver_status

! declare f2py directives
!F2PY intent(out) I_solver_status

     I_solver_status = solver_is_ready_gardenia
     if ( I_solver_status == 0 ) then
         call s_print_error('cat_solver_status','sorry, the current solver is not ready!')
     endif ! back if ( I_solver_status == 0 ) block

     return
  end subroutine cat_solver_status

!!========================================================================
!!>>> flow control subroutines                                         <<<
!!========================================================================

# if !defined (PYAPI)

!!>>> cat_init_ctqmc: initialize the ctqmc quantum impurity solver
!!>>> fortran version
  subroutine cat_init_ctqmc(I_mpi, I_solver)
     use capi, only : T_mpi, T_segment_gardenia

     use control ! ALL

     implicit none

! external arguments
! type structure of mpi
     type (T_mpi), intent(in) :: I_mpi

! type structure of generic solver
     type (T_segment_gardenia), intent(in) :: I_solver

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
     issus  = I_solver%issus
     isvrt  = I_solver%isvrt
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
     mune   = I_solver%mune
     beta   = I_solver%beta
     part   = I_solver%part
     alpha  = I_solver%alpha

! print the running header for continuous time quantum Monte Carlo quantum
! impurity solver and dynamical mean field theory self-consistent engine
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_header()
     endif ! back if ( myid == master ) block

! allocate memory and initialize
     call ctqmc_setup_array()

! prepare initial hybridization function, init self-consistent iteration
     call ctqmc_selfer_init()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_summary()
     endif ! back if ( myid == master ) block

     return
  end subroutine cat_init_ctqmc

# else   /* PYAPI */

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

! declare f2py directives
!F2PY intent(in) my_id
!F2PY intent(in) num_procs

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

! allocate memory and initialize
     call ctqmc_setup_array()

! prepare initial hybridization function, init self-consistent iteration
     call ctqmc_selfer_init()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then ! only master node can do it
         call ctqmc_print_summary()
     endif ! back if ( myid == master ) block

     return
  end subroutine cat_init_ctqmc

# endif  /* PYAPI */

!!>>> cat_exec_ctqmc: execute the ctqmc quantum impurity solver
  subroutine cat_exec_ctqmc(iter)
     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! declare f2py directives
!F2PY intent(in) iter

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

!!========================================================================
!!>>> data setter subroutines                                          <<<
!!========================================================================

!!>>> cat_set_hybf: setup the hybridization function
  subroutine cat_set_hybf(size_t, hybf_t)
     use control, only : norbs
     use control, only : mfreq
     use context, only : hybf

     implicit none

! external arguments
! size of hybf
     integer, intent(in)    :: size_t

! hybridization function
     complex(8), intent(in) :: hybf_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(in) hybf_t
!F2PY depend(size_t) hybf_t

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

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(in) symm_t
!F2PY depend(size_t) symm_t

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
     use context, only : eimp

     implicit none

! external arguments
! size of eimp
     integer, intent(in) :: size_t

! impurity level
     real(8), intent(in) :: eimp_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(in) eimp_t
!F2PY depend(size_t) eimp_t

! check whether size_t is correct
     if ( size_t /= size(eimp) ) then
         call s_print_error('cat_set_eimp','wrong dimension size of eimp_t')
     endif ! back if ( size_t /= size(eimp) ) block

! copy data
     eimp = eimp_t

     return
  end subroutine cat_set_eimp

!!>>> cat_set_ktau: setup the screening function and its first derivates
!!>>> note: the gardenia code does not support this function now
  subroutine cat_set_ktau(size_t, ktau_t, ptau_t)
     implicit none

! external arguments
! size of ktau
     integer, intent(in) :: size_t

! screening function K(\tau)
     real(8), intent(in) :: ktau_t(size_t)

! first derivate of screening function K'(\tau)
     real(8), intent(in) :: ptau_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(in) ktau_t
!F2PY intent(in) ptau_t
!F2PY depend(size_t) ktau_t
!F2PY depend(size_t) ptau_t

! to avoid the warning from compiler
     call s_assert( size(ktau_t) == size_t )
     call s_assert( size(ptau_t) == size_t )
     call s_print_error('cat_set_ktau','sorry, this feature is not supported')

     return
  end subroutine cat_set_ktau

!!>>> cat_set_uumat: setup the Coulomb interaction matrix
  subroutine cat_set_uumat(size_t, uumat_t)
     use control, only : norbs
     use context, only : uumat

     implicit none

! external arguments
! size of uumat
     integer, intent(in) :: size_t

! Coulomb interaction matrix
     real(8), intent(in) :: uumat_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(in) uumat_t
!F2PY depend(size_t) uumat_t

! check whether size_t is correct
     if ( size_t /= size(uumat) ) then
         call s_print_error('cat_set_uumat','wrong dimension size of uumat_t')
     endif ! back if ( size_t /= size(uumat) ) block

! copy data
     uumat = reshape(uumat_t,(/norbs,norbs/))

     return
  end subroutine cat_set_uumat

!!========================================================================
!!>>> data getter subroutines                                          <<<
!!========================================================================

!!>>> cat_get_grnf: extract the impurity green's function
  subroutine cat_get_grnf(size_t, grnf_t)
     use control, only : norbs
     use control, only : mfreq
     use context, only : grnf

     implicit none

! external arguments
! size of grnf
     integer, intent(in)     :: size_t

! impurity green's function
     complex(8), intent(out) :: grnf_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(out) grnf_t
!F2PY depend(size_t) grnf_t

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
     use control, only : norbs
     use control, only : mfreq
     use context, only : sig2

     implicit none

! external arguments
! size of sigf
     integer, intent(in)     :: size_t

! self-energy function
     complex(8), intent(out) :: sigf_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(out) sigf_t
!F2PY depend(size_t) sigf_t

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
     use control, only : norbs
     use context, only : nmat

     implicit none

! external arguments
! size of nmat
     integer, intent(in)  :: size_t

! occupation number
     real(8), intent(out) :: nmat_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(out) nmat_t
!F2PY depend(size_t) nmat_t

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
     use control, only : norbs
     use context, only : nnmat

     implicit none

! external arguments
! size of nnmat
     integer, intent(in)  :: size_t

! double occupation number
     real(8), intent(out) :: nnmat_t(size_t)

! declare f2py directives
!F2PY intent(in) size_t
!F2PY intent(out) nnmat_t
!F2PY depend(size_t) nnmat_t

! check whether size_t is correct
     if ( size_t /= size(nnmat) ) then
         call s_print_error('cat_get_nnmat','wrong dimension size of nnmat_t')
     endif ! back if ( size_t /= size(nnmat) ) block

! copy data
     nnmat_t = reshape(nnmat, (/norbs*norbs/))

     return
  end subroutine cat_get_nnmat
