# Parameter: nflip

**Definition**

> Flip period for spin up and spin down states.

**Type**

> Integer

**Default value**

> 20000

**Component**

> ALL

**Behavior**

Care must be taken to prevent the system from being trapped in a state which breaks a symmetry of local Hamiltonian when it should not be. To avoid un-physical trapping, we introduce "flip" moves, which exchange the operators corresponding, for example, to up and down spins in a given band.

In this code, nowadays the following flip schemes are supported:

* ```cflip``` = 1, flip intra-orbital spins one by one,

* ```cflip``` = 2, flip intra-orbital spins globally.

Here ```cflip``` is an internal variable, instead of an input parameters, you can not setup it in the *solver.ctqmc.in* file.

So the question is: how to control ```cflip``` via the *nflip* parameter?

Now we use the sign of *nflip* to control flip schemes

* *nflip* = 0, means infinite long period to do flip

* *nflip* > 0, flip intra-orbital spins one by one (90%) and globally (10%)

* *nflip* < 0, flip intra-orbital spins globally (90%) and one by one (10%)

**Comment**

If *nflip /= 0*, the absolute value of *nflip* is the flip period.
