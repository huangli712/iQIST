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
            "Build iQIST at multiple steps" => Any[
                "README" => "ch02/multi.md",
                "Build base library" => "ch02/base.md",
                "Build application programming interfaces" => "ch02/apis.md",
                "Build quantum impurity solvers" => "ch02/solvers.md",
                "Build applications" => "ch02/apps.md",
                "Build atomic eigenvalue problem solver" => "ch02/atomic.md",
                "Build auxiliary tools" => "ch02/tools.md",
                "Build libraries for Fortran "=> "ch02/fortran.md",
                "Build modules for Python" => "ch02/python.md",
            ],
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
            "Standard input files" => Any[
                "README" => "ch04/input.md",
                "solver.ctqmc.in" => "ch04/in_ctqmc.md",
                "solver.hfqmc.in" => "ch04/in_hfqmc.md",
                "solver.umat.in" => "ch04/in_umat.md",
                "solver.eimp.in" => "ch04/in_eimp.md",
                "solver.anydos.in" => "ch04/in_anydos.md",
                "solver.ktau.in" => "ch04/in_ktau.md",
                "atom.cix" => "ch04/in_atom.md",
            ],
            "Standard output files" => Any[
                "README" => "ch04/output.md",
                "Terminal output" => "ch04/out_term.md",
                "solver.green.dat" => "ch04/out_green.md",
                "solver.weiss.dat" => "ch04/out_weiss.md",
                "solver.hybri.dat" => "ch04/out_hybri.md",
                "solver.grn.dat" => "ch04/out_grn.md",
                "solver.wss.dat" => "ch04/out_wss.md",
                "solver.hyb.dat" => "ch04/out_hyb.md",
                "solver.sgm.dat" => "ch04/out_sgm.md",
                "solver.hub.dat" => "ch04/out_hub.md",
                "solver.hist.dat" => "ch04/out_hist.md",
                "solver.prob.dat" => "ch04/out_prob.md",
                "solver.nmat.dat" => "ch04/out_nmat.md",
                "solver.kmat.dat" => "ch04/out_kmat.md",
                "solver.lmat.dat" => "ch04/out_lmat.md",
                "solver.schi.dat" => "ch04/out_schi.md",
                "solver.sfom.dat" => "ch04/out_sfom.md",
                "solver.ochi.dat" => "ch04/out_ochi.md",
                "solver.ofom.dat" => "ch04/out_ofom.md",
                "solver.twop.dat" => "ch04/out_twop.md",
                "solver.vrtx.dat" => "ch04/out_vrtx.md",
                "solver.pair.dat" => "ch04/out_pair.md",
                "solver.diag.dat" => "ch04/out_diag.md",
                "solver.kernel.dat" => "ch04/out_kern.md",
                "solver.status.dat" => "ch04/out_stat.md",
            ],
            "Parameters" => Any[
                "README" => "ch04/parameters.md",
                "isscf" => "ch04/p_isscf.md",
                "issun" => "ch04/p_issun.md",
                "isspn" => "ch04/p_isspn.md",
                "isbin" => "ch04/p_isbin.md",
                "isort" => "ch04/p_isort.md",
                "issus" => "ch04/p_issus.md",
                "isvrt" => "ch04/p_isvrt.md",
                "isscr" => "ch04/p_isscr.md",
                "ifast" => "ch04/p_ifast.md",
                "itrun" => "ch04/p_itrun.md",
                "nband" => "ch04/p_nband.md",
                "nspin" => "ch04/p_nspin.md",
                "norbs" => "ch04/p_norbs.md",
                "ncfgs" => "ch04/p_ncfgs.md",
                "nzero" => "ch04/p_nzero.md",
                "niter" => "ch04/p_niter.md",
                "lemax" => "ch04/p_lemax.md",
                "legrd" => "ch04/p_legrd.md",
                "chmax" => "ch04/p_chmax.md",
                "chgrd" => "ch04/p_chgrd.md",
                "mkink" => "ch04/p_mkink.md",
                "mstep" => "ch04/p_mstep.md",
                "mfreq" => "ch04/p_mfreq.md",
                "nffrq" => "ch04/p_nffrq.md",
                "nbfrq" => "ch04/p_nbfrq.md",
                "nfreq" => "ch04/p_nfreq.md",
                "nsing" => "ch04/p_nsing.md",
                "ntime" => "ch04/p_ntime.md",
                "nvect" => "ch04/p_nvect.md",
                "nleja" => "ch04/p_nleja.md",
                "npart" => "ch04/p_npart.md",
                "nflip" => "ch04/p_nflip.md",
                "ntherm" => "ch04/p_ntherm.md",
                "nsweep" => "ch04/p_nsweep.md",
                "nwrite" => "ch04/p_nwrite.md",
                "nclean" => "ch04/p_nclean.md",
                "nmonte" => "ch04/p_nmonte.md",
                "ncarlo" => "ch04/p_ncarlo.md",
                "U" => "ch04/p_u.md",
                "Uc" => "ch04/p_uc.md",
                "Uv" => "ch04/p_uv.md",
                "Jz" => "ch04/p_jz.md",
                "Js" => "ch04/p_js.md",
                "Jp" => "ch04/p_jp.md",
                "lc" => "ch04/p_lc.md",
                "wc" => "ch04/p_wc.md",
                "mune" => "ch04/p_mune.md",
                "beta" => "ch04/p_beta.md",
                "part" => "ch04/p_part.md",
                "alpha" => "ch04/p_alpha.md",
            ],
       ]

ch05 = Any[ 
            "README" => "ch05/README.md",
            "DFT + DMFT" => "ch05/dft_dmft.md",
            "Ladder dual fermions" => "ch05/ladder.md",
       ]

ch06 = Any[
            "README" => "ch06/README.md",
            "Standard input files" => Any[
                "README" => "ch06/input.md",
                "atom.config.in" => "ch06/in_atom.md",
                "atom.cmat.in" => "ch06/in_cmat.md",
                "atom.emat.in" => "ch06/in_emat.md",
                "atom.tmat.in" => "ch06/in_tmat.md",
            ],
            "Standard output files" => Any[
                "README" => "ch06/output.md",
                "Terminal output" => "ch06/out_term.md",
                "solver.umat.in" => "ch06/out_umat1.md",
                "atom.fock.dat" => "ch06/out_fock.md",
                "atom.tmat.dat" => "ch06/out_tmat.md",
                "atom.emat.dat" => "ch06/out_emat.md",
                "atom.umat.dat" => "ch06/out_umat2.md",
                "atom.eigval.dat" => "ch06/out_val.md",
                "atom.eigvec.dat" => "ch06/out_vec.md",
                "atom.sector.dat" => "ch06/out_sector.md",
                "atom.cix" => "ch06/out_cix.md",
            ],
            "Parameters" => Any[
                "README" => "ch06/parameters.md",
                "ibasis" => "ch06/p_ibasis.md",
                "ictqmc" => "ch06/p_ictqmc.md",
                "icu" => "ch06/p_icu.md",
                "icf" => "ch06/p_icf.md",
                "isoc" => "ch06/p_isoc.md",
                "nband" => "ch06/p_nband.md",
                "nspin" => "ch06/p_nspin.md",
                "norbs" => "ch06/p_norbs.md",
                "ncfgs" => "ch06/p_ncfgs.md",
                "nmini" => "ch06/p_nmini.md",
                "nmaxi" => "ch06/p_nmaxi.md",
                "Uc" => "ch06/p_uc.md",
                "Uv" => "ch06/p_uv.md",
                "Jz" => "ch06/p_jz.md",
                "Js" => "ch06/p_js.md",
                "Jp" => "ch06/p_jp.md",
                "Ud" => "ch06/p_ud.md",
                "Jh" => "ch06/p_jh.md",
                "mune" => "ch06/p_mune.md",
                "lambda" => "ch06/p_lambda.md",
            ],
       ]

ch07 = Any[
            "README" => "ch07/README.md",
            "Maximum entropy method" => "ch07/mem.md",
            "Stochastic analytical continuation" => "ch07/sac.md",
            "Analytical continuation for self-energy" => "ch07/swing.md",
            "Toolbox" => Any[
                "README" => "ch07/toolbox.md",
                "toolbox/makechi" => "ch07/chi.md",
                "toolbox/makedos" => "ch07/dos.md",
                "toolbox/makekra" => "ch07/kra.md",
                "toolbox/makescr" => "ch07/scr.md",
                "toolbox/makesig" => "ch07/sig.md",
                "toolbox/makestd" => "ch07/std.md",
                "toolbox/maketau" => "ch07/tau.md",
                "toolbox/makeups" => "ch07/ups.md",
            ],
            "Scripts" => Any[
                "README" => "ch07/script.md",
                "script/d_archive.sh" => "ch07/archive.md",
                "script/d_check.py" => "ch07/check.md",
                "script/d_clean.py" => "ch07/clean.md",
                "script/d_cmp.py" => "ch07/cmp.md",
                "script/d_sar.sh" => "ch07/sar.md",
                "script/d_trailing.sh" => "ch07/trailing.md",
                "script/u_animator.py" => "ch07/animator.md",
                "script/u_atomic.py" => "ch07/atomic.md",
                "script/u_ctqmc.py" => "ch07/ctqmc.md",
                "script/u_hfqmc.py" => "ch07/hfqmc.md",
                "script/u_reader.py" => "ch07/reader.md",
                "script/u_writer.py" => "ch07/writer.md",
            ],
       ]

ch08 = Any[
            "README" => "ch08/README.md",
            "Basic applications" => Any[
                "README" => "ch08/basic.md",
                "Hello iQIST!" => "ch08/hello.md",
                "Mott metal-insulator transition" => "ch08/mott.md",
            ],
            "Advanced applications I: Complex systems" => Any[
                "README" => "ch08/complex.md",
                General Coulomb interaction](ch09/general.md)
                Spin-orbital coupling](ch09/soc.md)
                Crystal field splitting](ch09/cfs.md)
                Retarded interaction and dynamical screening effect](ch09/screening.md)
            ],
            "Advanced applications II: Accurate measurements" => Any[
                "README" => "ch08/accurate.md",
                * [One-shot and self-consistent calculations](ch09/self.md)
                * [Data binning mode](ch09/binning.md)
                * [Imaginary-time Green's function](ch09/gtau.md)
                * [Matsubara Green's function and self-energy function](ch09/matsubara.md)
                * [Spin-spin correlation function and orbital-orbital correlation function](ch09/chi.md)
                * [Two-particle Green's function and vertex function](ch09/vertex.md)
            ],
            "Advanced applications III: Post-processing procedures" => Any[
                "ch08/post.md",
                * [Analytical continuation for imaginary-time Green's function](ch09/mem.md)
                * [Analytical continuation for Matsubara self-energy function](ch09/swing.md)
            ],
            "Practical exercises" => Any[
                "ch08/practical.md",
                * [Orbital-selective Mott transition in two-band Hubbard model](ch09/osmt.md)
                * [Orbital Kondo and spin Kondo effects in three-band Anderson impurity model](ch09/kondo.md)
            ],
            "Code validation" => "ch08/valid.md",
            "Successful stories" => "ch08/story.md",
       ]

ch09 = Any[
            "README" => "ch09/README.md",
            "Basic theory and methods" => "ch09/basic.md",
            "Algorithms" => "ch09/algo.md",
            "Codes" => "ch09/code.md",
       ]

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
        "Applications" => ch05,
        "Atomic eigenvalue problem solver" => ch06,
        "Auxiliary tools" => ch07,
        "iQIST in action" => ch08,
        "Inside iQIST" => ch09,
        "Appendix" => "appendix/README.md",
        "Glossary" => "glossary.md",
    ],
)
