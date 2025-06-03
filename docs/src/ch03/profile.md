# Profile The Codes

!!! note

    Trust me, if you don't want to get involved into the tedious development of the iQIST software package, you can skip this section.

Are you satisfying with the efficiency of the iQIST's quantum impurity solvers? Though we have spent a lot of time in improving the efficiency of the iQIST, there are still a lot of codes/algorithms to be optimized.

Well, how to optimize the iQIST software package is not a trivial task. But, it is very important to figure out the hot-spots of the code at first. Usually, we use the following approach to analyze the computational efficiency of iQIST.

1. Modify *iqist/build/make.inc*. Active the "-pg" compiler opinion
2. Rerun the codes.
3. After the execution is finished, you will find a *gmon.out* file in the working directory. Here, we assume that you are running the iQIST codes in a Linux system. If you are using the MacOS system, the situation may be a bit different.
4. Then you can use the GNU *gprof* tool to analyze the *gmon.out* file. You have to ensure the execute program is in the same directory. As for the usage of the *gprof* tool, please google or 'man' it by yourself.
5. Of course, you can use any other advanced GUI tools to do this job. But it is far beyond the scopes of this reference manual.
