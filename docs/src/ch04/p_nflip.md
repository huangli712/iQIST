### Parameter: nflip

**Definition**

Flip period for spin up and spin down states.

**Type**

Integer

**Default value**

20000

**Component**

ALL, except for the **DAISY** component.

**Behavior**

Care must be taken to prevent the system from being trapped in a state which breaks a symmetry of local Hamiltonian when it should not be. To avoid un-physical trapping, we introduce "flip" moves, which exchange the operators corresponding, for example, to up and down spins in a given band.

In this code, nowadays the following flip schemes are supported:

* ```cflip``` = 1, flip inter-orbital spins randomly,

* ```cflip``` = 2, flip intra-orbital spins one by one,

* ```cflip``` = 3, flip intra-orbital spins globally.

Here ```cflip``` is an internal variable, instead of an input parameters, you can not setup it in the *solver.ctqmc.in* file. 

So the question is: how to control ```cflip``` via the *nflip* parameter?

Now we use the sign of *nflip* to control flip schemes

* *nflip* = 0, means infinite long period to do flip

* *nflip* > 0, combine ```cflip``` = 2 (80%) and ```cflip``` = 3 (20%)

* *nflip* < 0, combine ```cflip``` = 1 (80%) and ```cflip``` = 3 (20%)

**Comment**

If *nflip /= 0*, the absolute value of *nflip* is the flip period.

!!! note

    When *cflip = 1*, the symmetry of all orbitals must be taken into consideration, otherwise the code may be trapped by a deadlock.