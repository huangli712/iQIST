# Truncation approximation

As discussed in previous section, although we have used GQNs to split the full Hilbert space with very large dimension into blocks with smaller dimensions [for cases such as 7-band systems with GQNs (``N``, ``J_{z}``) and 5-band systems with GQN (``N``)], the dimensions of some blocks are still too large and the numbers of blocks are too much so that it is still very expensive to evaluate the local trace. Haule proposed to discard some high-energy states because they are rarely visited[^1]. For example, for 7-band system with only 1 electron (like Ce metal), only states with occupancy ``N=0``, 1, 2 will be frequently visited, and states with occupancy ``N>2`` can be truncated completely to reduce the large Hilbert space to a very small one. Of course, this truncation approximation may cause some bias because a frequently visited state may be accessed via an infrequently visited state. Therefore, one should be cautious when adopting the truncation approximation, and for example run some convergence tests.

Currently, we adopted two truncation schemes in our codes. The first scheme relies on the occupation number. We just keep those states whose occupation numbers are close to the nominal valence and skip the other states, as shown in the above example. This scheme is quite robust if the charge fluctuations are small enough, such as in the case of a Mott insulating phase. Another scheme is to dynamically truncate the states with very low probability based on statistics which is recorded during the Monte Carlo sampling. This scheme is not very stable, so one needs to use it with caution.

**Reference**

[^1]: Kristjan Haule, Phys. Rev. B 75, 155113 (2007)
