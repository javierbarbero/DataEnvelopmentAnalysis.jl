using Documenter, DataEnvelopmentAnalysis

makedocs(sitename = "DataEnvelopment\nAnalysis.jl",
        authors = "Javier Barbero and José Luis Zofío.",
        pages = [
        "Home" => "index.md",
        "Technical Efficiency Models" => Any[
                "Radial Models" => "technical/radial.md",
                "Directional Distance Function Models" => "technical/directional.md",
                "Additive Models" => "technical/additive.md",
                "Generalized Distance Function Models" => "technical/generalizeddf.md"
                ],
        "Economic Efficiency Models" => Any[
                "Cost Models" => "economic/cost.md",
                "Revenue Models" => "economic/revenue.md",
                "Profit Models" => "economic/profit.md",
                "Profitability Models" => "economic/profitability.md"
                ],
        "Productivity Change Models" => Any[
                "Malmquist Index" => "productivity/malmquist.md"
                ],
        "Bibliography" => "bibliography.md"
        ],
        format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)
        
deploydocs(
    repo = "github.com/javierbarbero/DataEnvelopmentAnalysis.jl.git",
)
