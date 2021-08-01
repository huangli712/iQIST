### atom.eigval.dat

**Introduction**

The *atom.eigval.dat* file contains the eigenvalues of the atomic Hamiltonian. All of the eigenvalues are outputted.

**Format**

The *atom.eigval.dat* file has two different file formats. They depends on the *ictqmc* parameter.

* **File Format A**, should be compatible with the **CAMELLIA**, **BEGONIA**, and **LAVENDER** components
> NOTE:
>
> *ictqmc* = 0, 1.

* **File Format B**, should be compatible with the **MANJUSHAKA** and **PANSY** components
> NOTE:
> 
> *ictqmc* = 2, 3, 4, 5.

See [ictqmc](p_ictqmc.md) parameter for more details.

The file formats of the *atom.eigval.dat* file are a bit complex. It is unnecessary to understand its full details. If you are a maintainer or developer of the **JASMINE** component, please read the source codes in *iqist/src/tools/jasmine/atomic_dump.f90*.

**Code**

N/A