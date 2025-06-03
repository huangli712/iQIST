# Configure Your System

Please add the directory which includes the executable programs of iQIST software package to your system path. Usually, you can modify the *.bashrc* file in your home directory to reach this goal:

```shell
export PATH=iqist/build:$PATH
```

Of course, you have to copy the executable programs to the *iqist/build* at first. If you don't know how to modify the system environment variables, please consult your system administrator.

---

**Where Is The Executable Programs?**

>
> * NARCISSUS Component - iqist/src/ctseg/ctqmc
> * MANJUSHAKA Component - iqist/src/cthyb/ctqmc
> * JASMINE Component - iqist/src/atomic/atomic
>
