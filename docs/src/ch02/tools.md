### Build auxiliary tools

In fact, the auxiliary tools include the **JASMINE** and **HIBISCUS** components. You can employ the following commands to compile them.

**Method 1**:
```
$ cd iqist/build
$ editor make.sys
$ make base
$ make capi
$ make tool
$ ./x_setup.sh
```

!!! note 

    1. Here **editor** means any ascii text editor you prefer.
    2. If the base library and the application programming interfaces were already compiled successfully, then you can skip the following commands:
    ```sh
    $ make base
    $ make capi
    ```

> 3. If you want to build the **HIBISCUS** component only, please replace
```
$ make tool
```
with
```
$ make hibiscus
```
If you want to build the **JASMINE** component only, please use
```
$ make jasmine
```

After a few minutes, you will find *jasmine.x*, *entropy.x*, *sac.x*, etc., in the *iqist/build* directory. They are what you need, the executable programs for the **JASMINE** and **HIBISCUS** components. Then you can add them to the system path, and do your great research.

!!! note

    The **JASMINE** component has a single executable program, *jasmine.x*. However, the **HIBISCUS** component has a few executable programs, including *entropy.x*, *sac.x*, *mdos.x*, etc.

**Method 2**:

In this approach, we go to the directory of the **JASMINE** component at first, 'make' it as usual. And then we go to the directory of the **HIBISCUS** component, 'make' it. That is all.

```
$ cd iqist/src/tools/jasmine
$ make
$ cd ../hibiscus
$ make
```

Here, we assume that the compiling system was correctly configured, the base library and the application programming interfaces were already compiled successfully.

!!! note 

    Now the executable program will not be copied into the *iqist/build* directory. You have to go to the *iqist/build* directory, and execute
    ```
    $ ./x_setup.sh
    ```