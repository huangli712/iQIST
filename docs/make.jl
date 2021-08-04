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
                * [solver.ctqmc.in](ch04/in_ctqmc.md)
                * [solver.hfqmc.in](ch04/in_hfqmc.md)
                * [solver.umat.in](ch04/in_umat.md)
                * [solver.eimp.in](ch04/in_eimp.md)
                * [solver.anydos.in](ch04/in_anydos.md)
                * [solver.ktau.in](ch04/in_ktau.md)
                * [atom.cix](ch04/in_atom.md)
            ],
            "Standard output files" => "ch04/output.md",
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
            "Standard input files" => "ch06/input.md",
            "Standard output files" => "ch06/output.md",
            "Parameters" => "ch06/parameters.md",
       ]

ch07 = Any[
            "README" => "ch07/README.md",
            "Maximum entropy method" => "ch07/mem.md",
            "Stochastic analytical continuation" => "ch07/sac.md",
            "Analytical continuation for self-energy" => "ch07/swing.md",
            "Toolbox" => "ch07/toolbox.md",
            "Scripts" => "ch07/script.md",
       ]

ch08 = Any[
            "README" => "ch08/README.md",
            "Basic applications" => "ch08/basic.md",
            "Advanced applications I: Complex systems" => "ch08/complex.md",
            "Advanced applications II: Accurate measurements" => "ch08/accurate.md",
            "Advanced applications III: Post-processing procedures" => "ch08/post.md",
            "Practical exercises" => "ch08/practical.md",
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
