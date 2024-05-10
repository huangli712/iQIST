# Parameter: isbnd

**Definition**

> Key control flag, it is used to determine the symmetry of bands/orbitals.

**Type**

> Integer

**Default value**

> 1

**Component**

> ALL

**Behavior**

> There are two possible values for the *isbnd* parameter so far:
>
> * *isbnd* = 1, the bands are not symmetrized.
> * *isbnd* = 2, the bands are symmetrized according to symmetry matrix ```symm```.
>
> The quantum impurity solvers will symmetrize the relevant physical quantities according to the *isbnd* parameter.

**Comment**

> There is still one question. What's the symmetry matrix and where is it?
>
> In fact, we should call it symmetry vector, instead of symmetry matrix. It is an integer vector whose size is exact *norbs*. The quantum impurity solvers will use it to perform symmetrization. Let's suppose this vector is called ```symm```. For example, there is a 5-band model, the symmetry vector is as follows:
>
> ```fortran
> symm(01) = 1
> symm(02) = 2
> symm(03) = 3
> symm(04) = 4
> symm(05) = 4
> symm(06) = 4
> symm(07) = 5
> symm(08) = 5
> symm(09) = 5
> symm(10) = 1
> ```
>
> Guess which orbitals are degenerated?
>
> ehh, the orbitals 1 and 10, since ```symm(01)``` and ```symm(10)``` are all 1. Similarly, the orbitals 4, 5, and 6 are degenerated. The orbitals 7, 8, and 9 are degenerated as well.
>
> Yes, you are right. The quantum impurity solvers will consider that those orbitals which have the same symmetry number (the element of ```symm```) are degenerated.
>
> By default, the elements in ```symm``` vector are 1. In other words, if *isbnd* = 2, all of the orbitals are degenerated. If you want to re-assign the degenerated orbitals, you can edit the *solver.eimp.in* file. See [solver.eimp.in](in_eimp.md) for more details.
