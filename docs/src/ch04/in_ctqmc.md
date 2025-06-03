# solver.ctqmc.in

**Introduction**

The *solver.ctqmc.in* file is the only configuration file for the continuous-time quantum Monte Carlo impurity solver. You can setup the parameters for the quantum impurity solvers in it. Besides the Fortran API, this is the only way to setup the parameters for the quantum impurity solvers. But we have to emphasize that this file is optional. In other words, the CT-QMC impurity solvers can run without it.

There are many input parameters for the CT-QMC impurity solvers. But all of parameters have default values. If the *solver.ctqmc.in* file is absent, the quantum impurity solvers will use the default values. If the *solver.ctqmc.in* file is present, the quantum impurity solvers will read it, parse it, and apply the settings in it to initialize the quantum impurity solvers. The default values for the parameters are designed for a single-band Hubbard model. Thus in most cases, you need a *solver.ctqmc.in* file to override the default settings.

---

**Format**

All of the CT-QMC impurity solvers in the iQIST software package share the same *solver.ctqmc.in* file. In other words, you can exchange *solver.ctqmc.in* files between different CT-QMC quantum impurity solvers. The format of the *solver.ctqmc.in* file adopts the simple "key-value" style. The detailed rules are as follows:

(1) Anything after "#" and "!" character will be treated as comments and will be ignored completely.
```
example:
# this is a comment line
! this is a comment line
nband = 4 # this is in line comment
norbs = 8 ! this in line comment
```

(2) It is not case sensitive.
```
example:
Nband = 4
NORBS = 8
NspiN = 2
```

(3) The key and value pair is separated by "=" or ":" character.
```
example:
nband = 4 ! you can use nband : 4
norbs : 8 ! you can use norbs = 8
```

(4) Any space will be ignored. Any blank lines will be skipped as well.
```
example:
n b a n d = 4 ! it is valid
no   rb s = 8 ! it is valid
```

(5) You can only use one line to define one key-value pair.
```
example
nband = 4 norbs = 8  ! it is not valid
nband = 4, norbs = 8 ! it is not valid
nband = 4; norbs = 8 ! it is not valid
nband =              !
4                    ! it is not valid
```

(6) In the value part, now only integer, real(dp), logical, and character data type are supported.
```
example:
nband = 4        ! integer type
mune  = 4.0      ! real(dp) type
isscf = .true.   ! logical type, you can also use .false., T, F
model = anderson ! character type, do not use "" or '' characters to quote it
```

(7) In the value part, a vector is also support. The items in the vector  should be separated by "," character.
```
example:
nband = 1, 2, 3, 4                   ! 4 items
mune = 0.0, -1.0, 2.0                ! 3 items
isscf = .true., .true., F, T, .true. ! 5 items
model = anderson, hubbard            ! 2 items
```

(8) An empty input file is acceptable.

(9) If one key occurs in the input file for more than 1 times, only the last occurrence is recognized.

!!! warning

    We mention that the quantum impurity solvers will not check whether the settings in the *solver.ctqmc.in* file are reasonable and correct. It is the user's responsibility.

---

**Code**

N/A

!!! tip

    In the [Flink](https://github.com/huangli712/Flink) library, we implement a parser to parse the *solver.ctqmc.in* file.

!!! tip

    The [ZenGui](https://github.com/huangli712/Flink) code can be used to generate the *solver.ctqmc.in* file.
