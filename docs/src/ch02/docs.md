# Build Documentation

The documentations for the iQIST software package are written by using the *Markdown* language and the *Julia's Documenter.jl* package. You can type the following command in the terminal to build the documentations:

```shell
$ cd iqist/docs
$ rm -fr build
$ julia make.jl
```

Then the HTML-version documentation is generated in the *iqist/docs/build* directory. Or you can read the official reference manual in the following website:
    
```text
https://huangli712.github.io/projects/iqist/index.html
```

!!! tip

    For more information about the Julia language and the Documenter.jl package, please visit the following websites:
    * [Julia](https://julialang.org/)
    * [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl)
