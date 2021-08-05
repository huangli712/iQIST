### script/d_cmp.py

**Introduction**

The purpose of this script is to compare the output files between two cases when testing solvers. With this script, we can judge whether the new results generated with a new version of iQIST are consistent with the previous generated with an old iQIST code. In this script, the histogram of perturbation expansion order, ``G(\tau)``, ``G(i\omega_n)``, and  ``\Sigma(i\omega_n)`` are used as criteria.

**Type**

Python script

**Usage**

This script needs two parameters, i.e., the directories for old and new results.

```
$ ./d_cmp.py  directory_of_case_A directory_of_case_B
```

**Input**

N/A

**Output**

See the terminal output.

**Comment**

!!! note

    This script is used by the iQIST Developer Team *internally*.