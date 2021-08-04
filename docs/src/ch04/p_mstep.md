### Parameter: mstep

**Definition**

We can use the delayed Green's function updating algorithm[^1] to accelerate the Hirsch-Fye quantum Monte Carlo impurity solver (the **DAISY** component in the iQIST software package). Here, the *mstep* parameter denotes the maximum number of delayed update steps which is the key control parameter for this algorithm.

**Type**

Integer

**Default value**

16

**Component**

Only for the **DAISY** component.

**Behavior**

There are two possible choices for *mstep* parameter so far:

* *mstep* = 1, using traditional update algorithm, low efficiency.

* *mstep* > 1, using delayed update algorithm to improve the computational efficiency.

**Comment**

*mstep* = 16 is an optimal choice. Do not change it rashly.

**Reference**

[^1]: Phani K. V. V. Nukala, Thomas A. Maier, Michael S. Summers, Gonzalo Alvarez, and Thomas C. Schulthess, *Phys. Rev. B* **80**, 195111 (2009)