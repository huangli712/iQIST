### Fortran binding

**> JASMINE** component (atomic eigenvalue problem solver)

Now the following subroutines/functions are supported by the **JASMINE** component.

* *cat_init_atomic*()
* *cat_exec_atomic*()
* *cat_stop_atomic*()

See the simplified source codes for more details.

```fortran
!!>>> cat_init_atomic: initialize the atomic eigenvalue problem solver
!!>>> fortran version
  subroutine cat_init_atomic(I_solver)
     type (T_jasmine), intent(in) :: I_solver
     return
  end subroutine cat_init_atomic

!!>>> cat_exec_atomic: execute the atomic eigenvalue problem solver
  subroutine cat_exec_atomic()
     return
  end subroutine cat_exec_atomic

!!>>> cat_stop_atomic: stop the atomic eigenvalue problem solver
  subroutine cat_stop_atomic()
     return
  end subroutine cat_stop_atomic
```
---

**> CT-HYB** quantum impurity solvers

Now the following subroutines/functions are supported by various CT-HYB quantum impurity solvers.

* *cat_solver_id*()
* *cat_solver_status*()
* *cat_init_ctqmc*()
* *cat_exec_ctqmc*()
* *cat_stop_ctqmc*()
* *cat_set_hybf*()
* *cat_set_symm*()
* *cat_set_eimp*()
* *cat_set_ktau*()
* *cat_set_uumat*()
* *cat_get_grnf*()
* *cat_get_sigf*()
* *cat_get_nmat*()
* *cat_get_nnmat*()

See the simplified source codes for more details.

```fortran
!!========================================================================
!!>>> status query subroutines                                         <<<
!!========================================================================

!!>>> cat_solver_id: return the solver identity
  subroutine cat_solver_id(I_solver_id)
     integer, intent(out) :: I_solver_id
     return
  end subroutine cat_solver_id

!!>>> cat_solver_status: return the solver status
  subroutine cat_solver_status(I_solver_status)
     integer, intent(out) :: I_solver_status
     return
  end subroutine cat_solver_status

!!========================================================================
!!>>> flow control subroutines                                         <<<
!!========================================================================

!!>>> cat_init_ctqmc: initialize the ctqmc quantum impurity solver
!!>>> fortran version
  subroutine cat_init_ctqmc(I_mpi, I_solver)
     type (T_mpi), intent(in) :: I_mpi
     type (T_generic_solver), intent(in) :: I_solver
     return
  end subroutine cat_init_ctqmc

!!>>> cat_exec_ctqmc: execute the ctqmc quantum impurity solver
  subroutine cat_exec_ctqmc(iter)
     integer, intent(in) :: iter
     return
  end subroutine cat_exec_ctqmc

!!>>> cat_stop_ctqmc: stop the ctqmc quantum impurity solver
  subroutine cat_stop_ctqmc()
     return
  end subroutine cat_stop_ctqmc

!!========================================================================
!!>>> data setter subroutines                                          <<<
!!========================================================================

!!>>> cat_set_hybf: setup the hybridization function
  subroutine cat_set_hybf(size_t, hybf_t)
     integer, intent(in)    :: size_t
     complex(8), intent(in) :: hybf_t(size_t)
     return
  end subroutine cat_set_hybf

!!>>> cat_set_symm: setup the symmetry vector
  subroutine cat_set_symm(size_t, symm_t)
     integer, intent(in) :: size_t
     integer, intent(in) :: symm_t(size_t)
     return
  end subroutine cat_set_symm

!!>>> cat_set_eimp: setup the impurity level
  subroutine cat_set_eimp(size_t, eimp_t)
     integer, intent(in) :: size_t
     real(8), intent(in) :: eimp_t(size_t)
     return
  end subroutine cat_set_eimp

!!>>> cat_set_ktau: setup the screening function and its first derivates
  subroutine cat_set_ktau(size_t, ktau_t, ptau_t)
     integer, intent(in) :: size_t
     real(8), intent(in) :: ktau_t(size_t)
     real(8), intent(in) :: ptau_t(size_t)
     return
  end subroutine cat_set_ktau

!!>>> cat_set_uumat: setup the Coulomb interaction matrix
  subroutine cat_set_uumat(size_t, uumat_t)
     integer, intent(in) :: size_t
     real(8), intent(in) :: uumat_t(size_t)
     return
  end subroutine cat_set_uumat

!!========================================================================
!!>>> data getter subroutines                                          <<<
!!========================================================================

!!>>> cat_get_grnf: extract the impurity green's function
  subroutine cat_get_grnf(size_t, grnf_t)
     integer, intent(in)     :: size_t
     complex(8), intent(out) :: grnf_t(size_t)
     return
  end subroutine cat_get_grnf

!!>>> cat_get_sigf: extract the self-energy function
  subroutine cat_get_sigf(size_t, sigf_t)
     integer, intent(in)     :: size_t
     complex(8), intent(out) :: sigf_t(size_t)
     return
  end subroutine cat_get_sigf

!!>>> cat_get_nmat: extract the occupation number
  subroutine cat_get_nmat(size_t, nmat_t)
     integer, intent(in)  :: size_t
     real(8), intent(out) :: nmat_t(size_t)
     return
  end subroutine cat_get_nmat

!!>>> cat_get_nnmat: extract the double occupation number
  subroutine cat_get_nnmat(size_t, nnmat_t)
     integer, intent(in)  :: size_t
     real(8), intent(out) :: nnmat_t(size_t)
     return
  end subroutine cat_get_nnmat
```

---

**> HF-QMC** quantum impurity solvers

Now the following subroutines/functions are supported by various HF-QMC quantum impurity solvers.

* *cat_solver_id()*
* *cat_solver_status()*
* *cat_init_hfqmc()*
* *cat_exec_hfqmc()*
* *cat_stop_hfqmc()*
* *cat_set_wssf()*
* *cat_set_symm()*
* *cat_set_eimp()*
* *cat_set_ktau()*
* *cat_get_grnf()*
* *cat_get_sigf()*
* *cat_get_nmat()*
* *cat_get_nnmat()*

See the simplified source codes for more details.

```fortran
!!========================================================================
!!>>> status query subroutines                                         <<<
!!========================================================================

!!>>> cat_solver_id: return the solver identity
  subroutine cat_solver_id(I_solver_id)
     integer, intent(out) :: I_solver_id
     return
  end subroutine cat_solver_id

!!>>> cat_solver_status: return the solver status
  subroutine cat_solver_status(I_solver_status)
     integer, intent(out) :: I_solver_status
     return
  end subroutine cat_solver_status

!!========================================================================
!!>>> flow control subroutines                                         <<<
!!========================================================================

!!>>> cat_init_hfqmc: initialize the hfqmc quantum impurity solver
!!>>> fortran version
  subroutine cat_init_hfqmc(I_mpi, I_solver)
     type (T_mpi), intent(in) :: I_mpi
     type (T_daisy), intent(in) :: I_solver
     return
  end subroutine cat_init_hfqmc

!!>>> cat_exec_hfqmc: execute the hfqmc quantum impurity solver
  subroutine cat_exec_hfqmc(iter)
     integer, intent(in) :: iter
     return
  end subroutine cat_exec_hfqmc

!!>>> cat_stop_hfqmc: stop the hfqmc quantum impurity solver
  subroutine cat_stop_hfqmc()
     return
  end subroutine cat_stop_hfqmc

!!========================================================================
!!>>> data setter subroutines                                          <<<
!!========================================================================

!!>>> cat_set_wssf: setup the bath weiss's function
  subroutine cat_set_wssf(size_t, wssf_t)
     integer, intent(in)    :: size_t
     complex(8), intent(in) :: wssf_t(size_t)
     return
  end subroutine cat_set_wssf

!!>>> cat_set_symm: setup the symmetry vector
  subroutine cat_set_symm(size_t, symm_t)
     integer, intent(in) :: size_t
     integer, intent(in) :: symm_t(size_t)
     return
  end subroutine cat_set_symm

!!>>> cat_set_eimp: setup the impurity level
  subroutine cat_set_eimp(size_t, eimp_t)
     integer, intent(in) :: size_t
     real(8), intent(in) :: eimp_t(size_t)
     return
  end subroutine cat_set_eimp

!!>>> cat_set_ktau: setup the screening function and its first derivates
!!>>> note: the daisy code does not support this function now
  subroutine cat_set_ktau(size_t, ktau_t, ptau_t)
     integer, intent(in) :: size_t
     real(8), intent(in) :: ktau_t(size_t)
     real(8), intent(in) :: ptau_t(size_t)
     return
  end subroutine cat_set_ktau

!!========================================================================
!!>>> data getter subroutines                                          <<<
!!========================================================================

!!>>> cat_get_grnf: extract the impurity green's function
  subroutine cat_get_grnf(size_t, grnf_t)
     integer, intent(in)     :: size_t
     complex(8), intent(out) :: grnf_t(size_t)
     return
  end subroutine cat_get_grnf

!!>>> cat_get_sigf: extract the self-energy function
  subroutine cat_get_sigf(size_t, sigf_t)
     integer, intent(in)     :: size_t
     complex(8), intent(out) :: sigf_t(size_t)
     return
  end subroutine cat_get_sigf

!!>>> cat_get_nmat: extract the occupation number
  subroutine cat_get_nmat(size_t, nmat_t)
     integer, intent(in)  :: size_t
     real(8), intent(out) :: nmat_t(size_t)
     return
  end subroutine cat_get_nmat

!!>>> cat_get_nnmat: extract the double occupation number
  subroutine cat_get_nnmat(size_t, nnmat_t)
     integer, intent(in)  :: size_t
     real(8), intent(out) :: nnmat_t(size_t)
     return
  end subroutine cat_get_nnmat
```