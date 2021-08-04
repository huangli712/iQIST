### script/d_sar.sh

**Introduction**

The purpose of this script is to scan a file or directory, and then replace some characters with given characters. So we name it as *sar.sh* (**S**can **A**nd **R**eplace). We can use it to preprocess the *atom.config.in*/*solver.ctqmc.in* files in the *iqist/working/ctqmc/standard* directory.

**Type**

Bash shell script

**Usage**

There is not input parameter for this script. You can execute it without any parameters.

```
$./d_sar.sh
```

Before you start to use this shell script, you have to check and edit
carefully the string match pattern in it.

**Input**

N/A

**Output**

See the terminal output.

**Comment**

For Mac OS X system, the grammar for *sed* is (we don't generate backup)

```
sed -i '' ...
```

However, for Linux-based system, the grammar for *sed* is

```
sed -i ...
```

!!! note:

    This script is used by the iQIST Developer Team *internally*.