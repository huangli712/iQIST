using Documenter

makedocs(
    sitename="Flink",
    clean = false,
    authors = "Li Huang",
    format = Documenter.HTML(
        prettyurls = false,
        ansicolor = true,
    ),
    pages = [
        "Home" => "index.md",
        "Introduction" => "intro.md",
        "Installation" => "install.md",
        "Basic usage" => "usage.md",
        "Guide" => Any[
            "make.sys" => "guide/make.md",
            "Constants" => "guide/m_constants.md",
            "Linked list" => "guide/m_linkedlist.md",
            "Message passaging interface" => "guide/m_mpi.md",
            "Configuration parser" => "guide/m_parser.md",
            "Sparse matrix" => "guide/m_sparse.md",
            "Pseudorandom number generator" => "guide/m_spring.md",
            "Stack" => "guide/m_stack.md",
            "Analytical tetrahedron algorithm" => "guide/m_tetra.md",
            "Error and exception" => "guide/s_error.md",
            "Fourier transformation" => "guide/s_fourier.md",
            "Special functions" => "guide/s_function.md",
            "Integration" => "guide/s_integrator.md",
            "Spline interpolation" => "guide/s_spline.md",
            "Utility" => "guide/s_util.md",
            "Vector" => "guide/s_vector.md",
        ],
    ],
)
