# iQIST (Interacting Quantum Impurity Solver Toolkit)

The iQIST software package includes several quantum impurity solvers which implement the hybridization expansion version continuous-time quantum Monte Carlo algorithm and Hirsch-Fye quantum Monte Carlo algorithm, and corresponding preprocessing and postprocessing tools.

### WARNING

The iQIST is still in heavy development. The codes are extremely unstable. Some features are still experimental. Everything could be changed in the future release. We can not guarantee that it is bug free. So be careful when you are using it and verify your data again and again before you submit your calculated results in any peer-reviewed journal.

### Version

v0.4.0

### License

GNU General Public License Version 3

### Installation
* Full Installation
```sh
$ cd iqist/src/build
$ editor make.sys
$ make all
$ cd ../../bin
$ ./setup.sh
```

* Partial Installation
```sh
$ cd iqist/src/build
$ editor make.sys
$ make common
$ make api
$ make component (component could be azalea, gardenia, narcissus, etc.)
$ cd ../../bin
$ ./setup.sh
```

Enjoy it! 

If you want to know more about the compiling system implemented in the iQIST, please read the manual carefully.

### Document

see iQIST/doc/manual/ug.pdf (We are sorry. Currently this manual is far away from completeness).

### Development

The iQIST software package is developed and maintained by the iQIST Developer Team.

Find a bug? Want to contribute? Want new features? Great! Please contact with us as soon as possible.

### Reference

If you are using iQIST to do some studies and would like to publish your great works, it would be really appreciated if you can cite the following paper:

```sh
iQIST: An open source continuous-time quantum Monte Carlo impurity solver toolkit
Li Huang, Yilin Wang, Zi Yang Meng, Liang Du, Philipp Werner and Xi Dai
arXiv:1409.7573 (2014)
```

### Contact

```sh
Li Huang 
Department of Physics, Fribourg University, Switzerland
email: huangli712 at gmail.com
```

or

```sh
Yilin Wang
Institute of Physics, Chinese Academy of Sciences, Beijing, PRC
email: qhwyl2006 at 126.com 
```

