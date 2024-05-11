# u_movie.py

**Introduction**

The purpose of this script is to generate the animation movie using the data contained in the *solver.diag.dat* file.

**Type**

Python script

**Usage**

First, please edit the configuration parameters (such as inverse temperature ``\beta``, number of orbitals *norbs*, number of frames *nsweep/nwrite*, etc) in the script carefully, and then execute it in the terminal:

```
$ ./u_movie.py movie.mp4
```

Here *movie.mp4* is the output file. We can use the VLC to play it. If you don't supply any valid filename, the default output should be *diag.mp4*.

!!! note

    Usually, the default setting in this script is not good. You have to adjust them again and again to obtain a good mp4 video.

**Input**

* *solver.diag.dat*

All of the CT-HYB quantum impurity solvers in the iQIST software package can output the *solver.diag.dat* file if the *nwrite* parameter is correctly set.

**Output**

* *movie.mp4*

You can specify the output filename by yourself. The suffix must be *mp4*.

**Comment**

The generated animation movie can be played using the *VLC* video tool. Sometimes the Quicktime player in the Mac OS X can not play it correctly.
