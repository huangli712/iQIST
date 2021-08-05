### atom.eigvec.dat

**Introduction**

The *atom.eigvec.dat* file contains the eigenvectors of the atomic Hamiltonian. Only those elements which are larger than the threshold (> 1.0E-6) are outputted.

**Format**

The *atom.eigvec.dat* file has two different file formats. They depends on the *ictqmc* parameter.

* **File Format A**, should be compatible with the **CAMELLIA**, **BEGONIA**, and **LAVENDER** components

!!! note

    *ictqmc* = 0, 1.

* **File Format B**, should be compatible with the **MANJUSHAKA** and **PANSY** components

!!! note

    *ictqmc* = 2, 3, 4, 5.

See [ictqmc](p_ictqmc.md) parameter for more details.

The file formats of the *atom.eigvec.dat* file are a bit complex. It is unnecessary to understand its full details. If you are a maintainer or developer of the **JASMINE** component, please read the source codes in *iqist/src/tools/jasmine/atomic_dump.f90*.

**Code**

N/A