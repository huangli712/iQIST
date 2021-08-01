### toolbox/makedos

**Introduction**

This code will generate density of states $$\rho(\epsilon)$$ for the following lattice models:

* Gaussian density of states, $$d = \infty$$ cubic lattice
* $$d = 3$$ cubic lattice
* Semi-circular density of states, bethe lattice
* Lorentzian density of states

**Usage**

```
$ ./mdos
```

**Input**

See the terminal prompt.

**Output**

* *dos.gauss.dat*
* *dos.cubic.dat*
* *dos.bethe.dat*
* *dos.loren.dat*

**Comment**

The users can rename the output files to *solver.anydos.in* which can be read by the *ctqmc\_dmft\_anydos()* subroutine in the *ctqmc\_dmft.f90* file.

See also [solver.anydos.in](../ch04/in_anydos.md) for more details.