### Parameter: npart

**Definition**

Number of parts that the imaginary time axis is split.

**Type**

Integer

**Default value**

4

**Component**

Only for the **BEGONIA**, **LAVENDER**, **MANJUSHAKA**, and **PANSY** components.

**Behavior**

All operators in the imaginary time axis are grouped into *npart* parts according to their time values, in each Monte Carlo steps, only those changed parts are carefully dealt with, not all the parts. This trick can accelerate the codes significantly.

**Comment**

$$ 2\sqrt{3 \langle k \rangle \text{nband}} \sim 4\sqrt{3 \langle k \rangle \text{nband}}$$ may be the optimal value for *npart* to achieve maximum performance where $$\langle k \rangle$$ is the averaged perturbation expansion order. See [nband](p_nband.md) for more details.