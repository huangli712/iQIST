using Documenter

ch01 = Any[
            "README" => "ch01/README.md",
            "What's iQIST?" => "ch01/what.md",
            "Motivation" => "ch01/motivation.md",
            "Components" => "ch01/components.md",
            "Features" => "ch01/feature.md",
            "Software architecture" => "ch01/architecture.md",
            "Policy" => "ch01/policy.md",
       ]

ch02 = Any[
            "README" => "ch02/README.md",
            "Obtain" => "ch02/obtain.md",
            "Uncompress" => "ch02/uncompress.md",
            "Directory structures" => "ch02/directory.md",
            "Compiling environment" => "ch02/envir.md",
            "Compiling system" => "ch02/system.md",
            "Build iQIST at one step" => "ch02/full.md",
            "Build iQIST at multiple steps" => "ch02/multi.md",
       ]

ch03 = Any[
            "README" => "ch03/README.md",
            "Configure your system" => "ch03/config.md",
            "iQIST recipes" => "ch03/recipes.md",
            "Prepare input files" => "ch03/create.md",
            "Execute the codes" => "ch03/execute.md",
            "Monitor the codes" => "ch03/monitor.md",
            "Profile the codes" => "ch03/profile.md",
       ]

ch04 = Any[
            "README" => "ch04/README.md",
            "How to choose suitable quantum impurity solvers?" => "ch04/choose.md",
            "Standard input files" => "ch04/input.md",
            "Standard output files" => "ch04/output.md",
            "Parameters" => "ch04/parameters.md",
       ]

#=
* [Applications](ch05/README.md)
   * [DFT + DMFT](ch05/dft_dmft.md)
   * [Ladder dual fermions](ch05/ladder.md)
* [Atomic eigenvalue problem solver](ch06/README.md)
   * [Standard input files](ch06/input.md)
       * [atom.config.in](ch06/in_atom.md)
       * [atom.cmat.in](ch06/in_cmat.md)
       * [atom.emat.in](ch06/in_emat.md)
       * [atom.tmat.in](ch06/in_tmat.md)
   * [Standard output files](ch06/output.md)
       * [Terminal output](ch06/out_term.md)
       * [solver.umat.in](ch06/out_umat1.md)
       * [atom.fock.dat](ch06/out_fock.md)
       * [atom.tmat.dat](ch06/out_tmat.md)
       * [atom.emat.dat](ch06/out_emat.md)
       * [atom.umat.dat](ch06/out_umat2.md)
       * [atom.eigval.dat](ch06/out_val.md)
       * [atom.eigvec.dat](ch06/out_vec.md)
       * [atom.sector.dat](ch06/out_sector.md)
       * [atom.cix](ch06/out_cix.md)
   * [Parameters](ch06/parameters.md)
       * [ibasis](ch06/p_ibasis.md)
       * [ictqmc](ch06/p_ictqmc.md)
       * [icu](ch06/p_icu.md)
       * [icf](ch06/p_icf.md)
       * [isoc](ch06/p_isoc.md)
       * [nband](ch06/p_nband.md)
       * [nspin](ch06/p_nspin.md)
       * [norbs](ch06/p_norbs.md)
       * [ncfgs](ch06/p_ncfgs.md)
       * [nmini](ch06/p_nmini.md)
       * [nmaxi](ch06/p_nmaxi.md)
       * [Uc](ch06/p_uc.md)
       * [Uv](ch06/p_uv.md)
       * [Jz](ch06/p_jz.md)
       * [Js](ch06/p_js.md)
       * [Jp](ch06/p_jp.md)
       * [Ud](ch06/p_ud.md)
       * [Jh](ch06/p_jh.md)
       * [mune](ch06/p_mune.md)
       * [lambda](ch06/p_lambda.md)
* [Auxiliary tools](ch07/README.md)
   * [Maximum entropy method](ch07/mem.md)
   * [Stochastic analytical continuation](ch07/sac.md)
   * [Analytical continuation for self-energy](ch07/swing.md)
   * [Toolbox](ch07/toolbox.md)
       * [toolbox/makechi](ch07/chi.md)
       * [toolbox/makedos](ch07/dos.md)
       * [toolbox/makekra](ch07/kra.md)
       * [toolbox/makescr](ch07/scr.md)
       * [toolbox/makesig](ch07/sig.md)
       * [toolbox/makestd](ch07/std.md)
       * [toolbox/maketau](ch07/tau.md)
       * [toolbox/makeups](ch07/ups.md)
   * [Scripts](ch07/script.md)
       * [script/d_archive.sh](ch07/archive.md)
       * [script/d_check.py](ch07/check.md)
       * [script/d_clean.py](ch07/clean.md)
       * [script/d_cmp.py](ch07/cmp.md)
       * [script/d_sar.sh](ch07/sar.md)
       * [script/d_trailing.sh](ch07/trailing.md)
       * [script/u_animator.py](ch07/animator.md)
       * [script/u_atomic.py](ch07/atomic.md)
       * [script/u_ctqmc.py](ch07/ctqmc.md)
       * [script/u_hfqmc.py](ch07/hfqmc.md)
       * [script/u_reader.py](ch07/reader.md)
       * [script/u_writer.py](ch07/writer.md)
* [Application programming interfaces](ch08/README.md)
   * [Fortran APIs](ch08/fortran.md)
       * [Explicit Fortran types](ch08/for_type.md)
       * [Fortran binding](ch08/for_bind.md)
       * [Examples](ch08/for_ex.md)
   * [Python APIs](ch08/python.md)
       * [Python binding](ch08/py_bind.md)
       * [Examples](ch08/py_ex.md)
* [iQIST in action](ch09/README.md)
   * [Basic applications](ch09/basic.md)
       * [Hello iQIST!](ch09/hello.md)
       * [Mott metal-insulator transition](ch09/mott.md)
   * [Advanced applications I: Complex systems](ch09/complex.md)
       * [General Coulomb interaction](ch09/general.md)
       * [Spin-orbital coupling](ch09/soc.md)
       * [Crystal field splitting](ch09/cfs.md)
       * [Retarded interaction and dynamical screening effect](ch09/screening.md)
   * [Advanced applications II: Accurate measurement of physical observable](ch09/accurate.md)
       * [One-shot and self-consistent calculations](ch09/self.md)
       * [Data binning mode](ch09/binning.md)
       * [Imaginary-time Green's function](ch09/gtau.md)
       * [Matsubara Green's function and self-energy function](ch09/matsubara.md)
       * [Spin-spin correlation function and orbital-orbital correlation function](ch09/chi.md)
       * [Two-particle Green's function and vertex function](ch09/vertex.md)
   * [Advanced applications III: Post-processing procedures](ch09/post.md)
       * [Analytical continuation for imaginary-time Green's function](ch09/mem.md)
       * [Analytical continuation for Matsubara self-energy function](ch09/swing.md)
   * [Practical exercises](ch09/practical.md)
       * [Orbital-selective Mott transition in two-band Hubbard model](ch09/osmt.md)
       * [Orbital Kondo and spin Kondo effects in three-band Anderson impurity model](ch09/kondo.md)
   * [Library mode](ch09/library.md)
       * [Call iQIST from Fortran language](ch09/fortran.md)
       * [Call iQIST from Python language](ch09/python.md)
   * [Code validation](ch09/valid.md)
   * [Successful stories](ch09/story.md)
* [Inside iQIST](ch10/README.md)
   * [Basic theory and methods](ch10/basic.md)
       * [Quantum impurity model](ch10/qim.md)
       * [Principles of continuous-time quantum Monte Carlo algorithm](ch10/ct.md)
       * [Hybridization expansion](ch10/hyb.md)
   * [Algorithms](ch10/algo.md)
       * [Transition probability](ch10/tran.md)
       * [Hubbard-Holstein model](ch10/holstein.md)
       * [Dynamical screening effect](ch10/screening.md)
       * [Physical observable](ch10/obs.md)
       * [Orthogonal polynomial representation](ch10/ortho.md)
       * [Kernel polynomial method](ch10/kpm.md)
       * [Improved estimator for the self-energy function](ch10/sig.md)
       * [Fast matrix update](ch10/fast.md)
       * [Good quantum number, subspace, and symmetry](ch10/symmetry.md)
       * [Krylov subspace iteration](ch10/krylov.md)
       * [Newton-Leja polynomial interpolation](ch10/leja.md)
       * [Truncation approximation](ch10/truncation.md)
       * [Lazy trace evaluation](ch10/lazy.md)
       * [Skip listing algorithm](ch10/skip.md)
       * [Divide-and-conquer algorithm](ch10/dac.md)
       * [Sparse matrix tricks](ch10/sparse.md)
       * [Delayed update algorithm](ch10/delay.md)
       * [Atomic eigenvalue solver](ch10/atomic.md)
       * [Single particle basis](ch10/basis.md)
       * [Spin-orbit coupling](ch10/soc.md)
       * [Coulomb interaction matrix](ch10/coulomb.md)
       * [Maximum entropy method](ch10/mem.md)
       * [Stochastic analytical continuation](ch10/sac.md)
   * [Codes](ch10/code.md)
       * [Development platform](ch10/platform.md)
       * [Common service module library](ch10/csml.md)
       * [Common service subroutine library](ch10/cssl.md)
       * [Fast Fourier transformation](ch10/fft.md)
       * [Interpolation](ch10/inter.md)
       * [Random number generators](ch10/rng.md)
       * [Parallelization](ch10/parallel.md)
       * [A guide to the source codes of the CT-HYB components](ch10/struct.md)
       * [How to add new parameter?](ch10/new_param.md)
       * [How to add new observable?](ch10/new_obs.md)
=#

makedocs(
    sitename="iQIST",
    clean = false,
    authors = "Li Huang",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
    ),
    pages = [
        "Home" => "index.md",
        "Team" => "team.md",
        "Copyright" => "copy.md",
        "Dedication" => "thanks.md",
        "Introduction" => ch01,
        "Installation" => ch02,
        "Getting started" => ch03,
        "Quantum Monte Carlo impurity solvers" => ch04,
        "Appendix" => "appendix/README.md",
        "Glossary" => "glossary.md",
    ],
)
