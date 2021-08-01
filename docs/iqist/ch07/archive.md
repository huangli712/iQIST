### script/d_archive.sh

**Introduction**

We can use this script to generate a compressed archive from the current repo branch of the iQIST software package.

**Type**

Bash shell script

**Usage**

```
$ ./d_archive.sh
```

**Input**

N/A

**Output**

The name of the output archive should like this:

```
iqist_43e2cbb_1441276643.tar.gz
```

Here ```43e2cbb``` is the abbreviated commit hash tag, and ```1441276643``` is the UNIX timestamp when this commit was committed.

**Comment**

To use this script, the *git* tool must be installed. 

> NOTE:

> This script is used by the iQIST Developer Team *internally*.