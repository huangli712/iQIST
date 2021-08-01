### Python binding

The Python binding is derived directly from the corresponding Fortran binding using the *f2py* package. So it is not strange that they are almost identical. However, due to the limitations of the *f2py* package, you have to pay special attentions to the following three APIs:

* *cat_init_atomic*()
* *cat_init_ctqmc*()
* *cat_init_hfqmc*()

In the Fortran binding, the input parameters for the above subroutines are predefined Fortran types. But unfortunately, the *f2py* does not support them. So we have to customize a *f2py*-compatible version of them. In the following, we will only show the *f2py*-compatible version of them. As for the other subroutines, please check section [Fortran binding](for_bind.md).

**> JASMINE** component (atomic eigenvalue problem solver)

The definition for the *cat_init_atomic*() subroutine is as follows:

```fortran
!!>>> cat_init_atomic: initialize the atomic eigenvalue problem solver
!!>>> python version
  subroutine cat_init_atomic()
     return
  end subroutine cat_init_atomic
```

> NOTE:

> Now the user has to setup control parameters for the **JASMINE** component via the *atom.config.in* file. With the help of *script/u_atomic.py*, it is a trivial task.

---

**> CT-HYB** quantum impurity solvers

The definition for the *cat_init_ctqmc*() subroutine is as follows:

```fortran
!!>>> cat_init_ctqmc: initialize the ctqmc quantum impurity solver
!!>>> python version
  subroutine cat_init_ctqmc(my_id, num_procs)
     integer, intent(in) :: my_id
     integer, intent(in) :: num_procs
     return
  end subroutine cat_init_ctqmc
```

> NOTE:

> Now the user has to setup control parameters for the **CT-HYB** components via the *solver.ctqmc.in* file. With the help of *script/u_ctqmc.py*, it is a trivial task.

---

**> HF-QMC** quantum impurity solver

The definition for the *cat_init_hfqmc*() subroutine is as follows:

```fortran
!!>>> cat_init_hfqmc: initialize the hfqmc quantum impurity solver
!!>>> python version
  subroutine cat_init_hfqmc(my_id, num_procs)
     integer, intent(in) :: my_id
     integer, intent(in) :: num_procs
     return
  end subroutine cat_init_hfqmc
```

> NOTE:

> Now the user has to setup control parameters for the **HF-QMC** components via the *solver.hfqmc.in* file. With the help of *script/u_hfqmc.py*, it is a trivial task.