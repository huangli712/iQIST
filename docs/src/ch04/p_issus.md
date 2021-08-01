### Parameter: issus

**Definition**

It is a key control flag, which is used to determine whether we should measure the charge or spin susceptibility and which one we should measure.

**Type**

Integer

**Default value**

1

**Component**

Only for the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **CAMELLIA**, and **MANJUSHAKA** components.

**Behavior**

We just use the following algorithm to judge which susceptibility should be calculated:

* *issus* parameter is converted from a decimal representation to a binary representation at first. For example, 10$$_{10}$$ is converted to 1010$$_{2}$$, 15$$_{10}$$ is converted to 1111$$_{2}$$, etc.

* Then we examine the bits from rightmost to leftmost one by one. If it is 1, then we do the calculation specified by this bit. If it is 0, then we ignore the calculation specified by this bit. For example, we just use the second bit (from right side to left side) to represent the calculation of spin-spin correlation function. Therefore, if *issus* parameter is 10$$_{10}$$ (1010$$_{2}$$), we will calculate the spin-spin correlation function since the second bit is 1. Supposing *issus* parameter is 13$$_{10}$$ (1101$$_{2}$$), we will not calculate it as the second bit is 0.

The following are the definitions of bit representation [*p* is the bit order (from right to left side)]:

* *p* = 1, do nothing

* *p* = 2, calculate spin-spin correlation function (*in time space*)

* *p* = 3, calculate orbital-orbital correlation function (*in time space*)

* *p* = 4, calculate spin-spin correlation function (*in frequency space*)

* *p* = 5, calculate orbital-orbital correlation function (*in frequency space*)

* *p* = 6, calculate $$\langle k^2 \rangle - \langle k \rangle^2$$

* *p* = 7, calculate fidelity susceptibility matrix

* *p* = 8, reserved

* *p* = 9, reserved

Now let's consider a concrete example.

Supposed that *issus* =  469$$_{10}$$, then the corresponding binary representation is 111010101$$_{2}$$. And then we can easily obtain that:

* *p* = 9, bit = 1

* *p* = 8, bit = 1

* *p* = 7, bit = 1

* *p* = 6, bit = 0

* *p* = 5, bit = 1

* *p* = 4, bit = 0

* *p* = 3, bit = 1

* *p* = 2, bit = 0

* *p* = 1, bit = 1

According to the previous definitions about the bit representation, we can conclude that the quantum impurity solvers should calculate the orbital-orbital correlation function (in time space and frequency space) and the fidelity susceptibility simultaneously.

For the **CAMELLIA**, **LAVENDER**, and **MANJUSHAKA** components, the bits at *p* = 2, 3, 4, and 5 are not ready.

**Comment**

The spin-spin correlation function will be written into the *solver.schi.dat* and *solver.sfom.dat* files. The orbital-orbital correlation function will be written into the *solver.ochi.dat* and *solver.ofom.dat* files. The kinetic energy fluctuation will be written into the *solver.kmat.dat* file. The fidelity susceptibility will be written into the *solver.lmat.dat* file.

The Green's function $$G(\tau)$$ and $$G(i\omega_n)$$, hybridization function $$\Delta(\tau)$$ and $$\Delta(i\omega_n)$$, and self-energy function $$\Sigma(i\omega_n)$$ will be calculated and output by the quantum impurity solvers in the iQIST software package implicitly. However, the susceptibilities and high-order correlation functions won't. You have to use the *issus* and *isvrt* parameters to active the related calculations.

See [isvrt](p_isvrt.md) for more details.