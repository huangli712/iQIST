# Parameter: nfreq

**Definition**

> Number of Matsubara frequency points sampling by continuous-time quantum Monte Carlo quantum impurity solver.

**Type**

> Integer

**Default value**

> 128

**Component**

> ALL

**Behavior**

> Determine the size of involved arrays. The larger *nfreq*, the more computational time is needed.

**Comment**

> There are *mfreq* matsubara frequency points in total. The physical observables on the first *nfreq* points are sampled directly by the Monte Carlo algorithm, however, the rest (*mfreq - nfreq + 1* points) values are evaluated by using Hubbard-I approximation or the other algorithm.
