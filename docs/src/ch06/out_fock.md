### atom.fock.dat

**Introduction**

The *atom.fock.dat* file contains the information about the Fock basis (occupation number basis).

**Format**

The format of the *atom.fock.dat* file is as follows:

---

*column 1*: unified index ``i``, integer

*column 2*: decimal number, integer

*column 3*: index of the Fock state, integer

*column 4*: binary representation of the Fock state, integer vector

---

**Code**

The corresponding Fortran code block for the writing of the *atom.fock.dat* file is as follows:

```fortran
! open file atom.fock.dat to write
     open(mytmp, file='atom.fock.dat', form='formatted', status='unknown')

! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | decimal | index | binary'
     write(mytmp,'(75a1)') dash ! dashed line

! write the data
     do i=1,ncfgs
         write(mytmp,'(i6)',advance='no') i
         write(mytmp,'(i6)',advance='no') dec_basis(i)
         write(mytmp,'(i6)',advance='no') ind_basis(dec_basis(i))
         write(mytmp,'(4X,14i1)') bin_basis(:,i)
     enddo ! over i={1,ncfgs} loop

! close data file
     close(mytmp)
```