# How To Choose Suitable Quantum Impurity Solvers?

At a first glance, you may feel puzzled why there are so many continuous-time quantum impurity solvers in the iQIST software package. What are the differences between them? Do we really need them? How to choose a suitable quantum impurity solvers for a given strongly correlated electron problem?

Relax! In the following we will uncover the secrets.

---

**How Many Quantum Monte Carlo Impurity Solvers Are There?**

Hmm, let me think. **Two**. Now the answer is two. These are the newest data. We are not sure whether there is new quantum impurity solver in the future.

The two CT-HYB quantum impurity solvers are as follows:

* **NARCISSUS** (in *iqist/src/ctseg*)
* **MANJUSHAKA** (in *iqist/src/cthyb*)

You probably have found that we always used the flowers to name the quantum impurity solvers (and the other projects). Well, so many flowers. I think at least you will like one of them.

![narcissus image](../assets/narcissus.jpg)

**Figure 1** | Narcissus (source: internet).

![manjushaka image](../assets/manjushaka.jpg)

**Figure 2** | Manjushaka (source: internet).

---

**Do We Really Need These Quantum Impurity Solvers?**

This problem is a bit complicated and equivalent to why we had designed and implemented so many quantum impurity solvers in the iQIST software package.

The answer is "Yes".

First, the state-of-the-art continuous-time quantum impurity solvers become more and more complex. In order to achieve high performance, many tricky algorithms, numerical methods are invented. Some of them conflict with each other. It is almost impossible to implement all of them in a single program. Second, the quantum impurity models in strongly correlated systems are extremely complicated and polytropic. It is also impossible to design such a perfect impurity solver which can solve all of the impurity problems. Third, each quantum impurity solver in the iQIST software package is designed for a specific impurity problem. These quantum impurity solvers share similar computational kernel, but they are highly optimized for the given impurity model.

The price for maintaining nine quantum impurity solvers is high, but it is worth it.

---

**What Are The Differences Between These Quantum Impurity Solvers?**

**NARCISSUS** component
>
>* Model: Density-density interaction, retarded interaction, no SOC.
>* Algorithm: Segment representation.
>* Feature: Full-fledged.
>* Scenario: Used in real research, especially in E-DMFT calculations.
>

**MANJUSHAKA** component

>
>* Model: General interaction, SOC.
>* Algorithm: General matrix representation + subspace + divide-and-conquer algorithm + >* Lazy trace evaluation + dynamical truncation + skip listing algorithm.
>* Feature: Full-fledged.
>* Scenario: Used in real research.
>

**See also**:

* [Features](../ch01/feature.md) // Exact features for the quantum impurity solvers.

---

**How To Choose Suitable Quantum Impurity Solvers For A Given Impurity Problem?**

OK, now let's return to the original problem. How to choose a suitable quantum impurity solver? It depends on the problems what you face. In the following, we will provide some guidelines.
