### toolbox/makestd

**Introduction**

The *makestd* code is often used to postprocess the self-energy function data to generate suitable input files for the **HIBISCUS**/swing code.

At first, you have to run the CT-HYB quantum impurity solver for many times to generate a series of *solver.sgm.dat* file. And then this code is used to deal with these *solver.sgm.dat* files to generate *std.sgm.dat* which is just what the **HIBISCUS**/swing code needs.

**Usage**

```
$ ./mstd
```

**Input**

* See the terminal prompt
* *solver.sgm.dat*.* (necessary)

**Output**

* *std.sgm.dat*

**Comment**

N/A