### toolbox/maketau

**Introduction**

The *maketau* code is often used to convert *solver.green.dat.``*``* or *solver.green.dat* files to *tau.grn.dat* file, which contains the necessary input data for the **HIBISCUS**/entropy or **HIBISCUS**/stoch codes.

**Usage**

```
./mtau
```

**Input**

* See the terminal prompt
* *solver.green.dat* or *solver.green.dat*.* (necessary)

!!! info

    About *nskip* control parameter:
    If *ctqmc* = 2 or 4, then *nskip* must be 1. If *ctqmc* = 1 or 3, then *nskip* could be positive integer. Be careful, *nskip* can not be arbitrary integer. Notice that mod(*ntime* - 1, *nskip*) must be 0, or else the obtained *tau.grn.dat* should be wrong.
    
    About *solver.green.dat.``*``* files:
    In order to obtain *solver.green.dat.``*``* files, you have to run the CT-HYB/HF-QMC codes for several times, and rename the *solver.green.dat* file to *solver.green.dat.``*``* file manually.

**Output**

* *tau.grn.dat*

**Comment**

N/A