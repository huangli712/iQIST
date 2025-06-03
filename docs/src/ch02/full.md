# Build iQIST At One Step

!!! tip

    The iQIST software package depends on the [Flink](https://github.com/huangli712/Flink) library, so please make sure that the Flink library is correctly installed and the environment variable **FLINK** is correctly setup. See [Compiling Environment](envir.md) for more details.

When the system is ready and the compiling environment is configured correctly, you can try to compile the iQIST software package right now. Let's start now.

```shell
$ cd iqist/build
$ editor make.inc
$ make all
```

!!! note

    Here **editor** means any ascii text editor you prefer. The *vim* is our favorite.

After a few minutes (it depends on the performance of your system), if everything is OK you will find some executable programs in the *iqist/src/* directory, which means the compiling is completed. Now all of the components of the iQIST software package are successfully built. Congratulations! You can start to enjoy using iQIST and begin with your great research!

However, don't get frustrated once the compiling fails. It is wise for you to type the following command in the terminal:

```
$ make clean
```

Then you go to check and correct the *make.inc* file, fix the other system-related issues, and repeat the above procedures. We are sure that you will solve this problem successfully.

Do you expect a more controllable compiling procedure? Some other great ideas are given in the next section.
