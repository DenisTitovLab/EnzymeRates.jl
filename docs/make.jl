using Documenter
using DocumenterCitations
using EnzymeRates

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style = :authoryear)

makedocs(;
    sitename = "EnzymeRates.jl",
    authors = "Denis Titov and contributors",
    modules = [EnzymeRates],
    doctest = true,
    checkdocs = :exports,
    plugins = [bib],
    format = Documenter.HTML(;
        canonical = "https://DenisTitovLab.github.io/EnzymeRates.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Deriving rate equations" => [
            "Rate equations from textbooks" => "deriving/textbooks.md",
            "Rapid equilibrium vs steady state" => "deriving/re_vs_ss.md",
            "The Cha / King–Altman algorithm" => "deriving/cha_king_altman.md",
            "Thermodynamic constraints" => "deriving/thermodynamic_constraints.md",
            "Ping-pong mechanisms" => "deriving/ping_pong.md",
            "Mechanisms with regulators" => "deriving/regulators.md",
            "Mechanisms with allosteric regulators" => "deriving/mwc_allostery.md",
        ],
        "Fitting rate equations" => [
            "Fitting tutorial & data format" => "fitting/tutorial.md",
            "Normalized vs absolute rate" => "fitting/normalized_vs_absolute.md",
            "Loss & optimizers" => "fitting/loss_and_optimizers.md",
        ],
        "Identifying the best rate equation" => [
            "Identify tutorial" => "identify/tutorial.md",
            "Model selection" => "identify/model_selection.md",
            "Running in parallel" => "identify/parallel.md",
            "The enumeration engine" => "identify/enumeration_engine.md",
        ],
        "Developer / Architecture" => "developer.md",
        "API Reference" => "api.md",
        "Roadmap" => "roadmap.md",
        "References" => "references.md",
    ],
)

deploydocs(;
    repo = "github.com/DenisTitovLab/EnzymeRates.jl.git",
    devbranch = "main",
)
