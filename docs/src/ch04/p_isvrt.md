### Parameter: isvrt

**Definition**

It is a key control flag, which is used to determine whether we should measure the high order correlation functions and which one we should measure.

**Type**

Integer

**Default value**

1

**Component**

Only for the **GARDENIA**, **NARCISSUS**, **LAVENDER**, **CAMELLIA**, and **MANJUSHAKA** components.

**Behavior**

We just use the following algorithm to judge which correlation function should be calculated:

* *isvrt* parameter is converted from a decimal representation to a binary representation at first. For example, 10``_{10}`` is converted to 1010``_{2}``, 15``_{10}`` is converted to 1111``_{2}``, etc.

* Then we examine the bits from rightmost to leftmost one by one. If it is 1, then we do the calculation specified by this bit. If it is 0, then we ignore the calculation specified by this bit. For example, we just use the second bit (from right side to left side) to represent the calculation of two-particle Green's function. So, if *isvrt* parameter is 10``_{10}`` (1010``_{2}``), we will calculate the two-particle Green's function since the second bit is 1. If *isvrt* parameter is 13``_{10}`` (1101``_{2}``), we will not calculate it since the second bit is 0.

The following are the definitions of bit representation [*p* is the bit order (from right to left side)]:

* *p* = 1, do nothing

* *p* = 2, calculate two-particle Green's function and vertex function (standard algorithm)

* *p* = 3, calculate two-particle Green's function and vertex function (improved algorithm)

* *p* = 4, calculate particle-particle pair susceptibility

* *p* = 5, reserved

* *p* = 6, reserved

* *p* = 7, reserved

* *p* = 8, reserved

* *p* = 9, reserved

Now let's illustrate a concrete example.

Supposed that *isvrt* =  469``_{10}``, so the corresponding binary representation is 111010101``_{2}``. And then we can easily obtain that:

* *p* = 9, bit = 1

* *p* = 8, bit = 1

* *p* = 7, bit = 1

* *p* = 6, bit = 0

* *p* = 5, bit = 1

* *p* = 4, bit = 0

* *p* = 3, bit = 1

* *p* = 2, bit = 0

* *p* = 1, bit = 1

According to the previous definitions about the bit representation, we can conclude that the quantum impurity solvers should calculate the two-particle Green's function and vertex function using an improved algorithm.

If the bits in *p* = 2 or *p* = 3 are set, both the two-particle Green's and vertex functions are computed, but using two different algorithms. You should not set them to 1 at the same time. In order words, if you set the bit at *p* = 2 to 1, then the bit at *p* = 3 must be 0, and vice versa.

If the bit in *p* = 2 is set, the traditional (standard) algorithm is used. If the bit in *p* = 3 is set, the improved estimator for two-particle Green's function is used which should have better accuracy.

**Comment**

As for the two-particle Green's function and vertex function, if the standard algorithm is used, the results are written into the *solver.twop.dat* file. If improved algorithm is adapted, the data are written into the *solver.vrtx.dat* file. The pair susceptibility will be written into the *solver.pair.dat* file.

The Green's function ``G(\tau)`` and ``G(i\omega_n)``, hybridization function ``\Delta(\tau)`` and ``\Delta(i\omega_n)``, and self-energy function ``\Sigma(i\omega_n)`` will be calculated and output by the quantum impurity solvers in the iQIST software package implicitly. However, the susceptibilities and high-order correlation functions won't. You have to use the *issus* and *isvrt* parameters to active the related calculations.

The *isvrt* parameter has nothing to do with the *isort* parameter, but when the bit in *p* = 3 is set to be 1, a prefactor for the improved estimator (please see ```pref``` in *ctqmc\_context.f90*) should be calculated and used to improve the computational accuracy. On the other hand, if *isort* >= 4, the same prefactor will be calculated as well. You can not setup *isscr* = 2 at this time (only be notable for the **NARCISSUS** component).

See [isvrt](p_isvrt.md), [isort](p_isort.md), and [isscr](p_isscr.md) for more details.