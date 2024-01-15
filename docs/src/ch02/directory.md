# Directory structures

The uncompressed iQIST software package is somewhat complicated, but have a rich content. We think that it is helpful to introduce the detailed directory structure of it. As has been explained before, in the following, we always use the term **iqist** to denote the top folder which contains the uncompressed iQIST software package.

![dir image](../figures/dir.png)

**Figure** | The directory structures of the iQIST software package.

**Where is the compling system?**

```
iqist/build
```

**Where are the source codes?**

```
iqist/src
```

**Where are the quantum impurity solvers?**

* For CT-QMC quantum impurity solvers

```
iqist/src/ct_hyb1
iqist/src/ct_hyb2
```

**Where is the atomic eigenvalue problem solver?**

```
iqist/src/atomic
```

**Where are the auxiliary tools?**

```
iqist/src/tools
```

**Where are the documentations?**

```
iqist/docs
```

!!! warning

    The official reference manual is released on the following website:
    ```
    https://www.gitbook.com/book/huangli712/iqist/details
    ```

### Where are the tutorials?

```
iqist/tutor
```

### Where are the benchmark cases?

```
iqist/test
```
