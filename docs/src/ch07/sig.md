### toolbox/makesig

**Introduction**

To calculate real physical quantities, such as the optical conductivity, Seebeck coefficient, thermopower, etc., the self-energy function on the real axis is an essential input. With the Pade approximation[^1], we can convert the self-energy function from the Matsubara frequency to real frequency axis. We implemented the Pade approximation for $$\Sigma(i\omega_n)$$ in the **HIBISCUS** component.

The *makesig* code is often used to transform self-energy functions from Matsubara frequency representation to real frequency representation via the Pade approximation. The results are very sensitive to the data noises in the self-energy function. So we do not recommend to use this code to perform analytical continuation for the self-energy function. However, the **HIBISCUS**/swing code may be a better choice.

[^1]: H. Vidberg and J. Serene, J. Low Temp. Phys. 29, 179 (1977).

**Usage**

```
$ ./msig
```

**Input**

* See the terminal prompt
* *solver.sgm.dat* (necessary)

**Output**

* *sig.sgm.dat*

**Comment**

N/A