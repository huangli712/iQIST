# Random number generators

Fast, reliable, and long period pseudo-random number generators are a key factor for any Monte Carlo simulations. Currently, the most popular random number generator is the Mersenne Twister which was developed in 1998 by Matsumoto and Nishimura[^1]. Its name derives from the fact that its period length is chosen to be a Mersenne prime. In the iQIST software package, we implemented the commonly used version of Mersenne Twister, MT19937. It has a very long period of ``2^{19937}-1``. Of course, if we choose different parameter sets, its period can be shorter or longer.

The Mersenne Twister is a bit slow by today's standards. So in 2006, a variant of Mersenne Twister, the SIMD-oriented Fast Mersenne Twister (SFMT)[^2] was introduced. It was designed to be fast when it runs on 128-bit SIMD. It is almost twice as fast as the original Mersenne Twister and has better statistics properties. We also implemented it in the iQIST software package, and use it as the default random number generator.

**Reference**

[^1]: M. Matsumoto, T. Nishimura, *ACM Transactions on Modeling and Computer Simulation* **8** 3 (1998)

[^2]: Mutsuo Saito and Makoto Matsumoto, *Monte Carlo and Quasi-Monte Carlo Methods*, **Springer**, 607 -- 622 (2008)
