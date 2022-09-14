using Documenter, DataEnvelopmentAnalysis

# Solve nonlinear problem to display Ipopt initial message
begin
        X = [1; 2; 3];
        Y = [1; 1; 1];
        deagdf(X, Y, alpha = 0.5, rts = :VRS)
end
    
DocMeta.setdocmeta!(DataEnvelopmentAnalysis, :DocTestSetup, :(using DataEnvelopmentAnalysis); recursive=true)

makedocs(sitename = "DataEnvelopmentAnalysis",
        authors = "Javier Barbero and José Luis Zofío.",
        pages = [
        "Home" => "index.md",
        "Technical Efficiency Models" => Any[
                
                "Radial Models" => "technical/radial.md",
                "Radial Big Data Models" => "technical/radialbigdata.md",
                "Directional Distance Function Models" => "technical/directional.md",
                "Additive Models" => "technical/additive.md",
                "Generalized Distance Function Models" => "technical/generalizeddf.md",
                "Russell Models" => "technical/russell.md",
                "Enhanced Russell Graph Slack Based Measure" => "technical/enhancedrussell.md",
                "Modified Directional Distance Function" => "technical/modifiedddf.md",
                "Hölder Distance Function" => "technical/holder.md",
                "Reverse Directional Distance Function" => "technical/reverseddf.md",
                "Common functions for technical models" => "technical/commontechnical.md"
                ],
        "Economic Efficiency Models" => Any[
                "Cost Models" => "economic/cost.md",
                "Revenue Models" => "economic/revenue.md",
                "Profit Models" => "economic/profit.md",
                "Profitability Models" => "economic/profitability.md",
                "Common functions for economic models" => "economic/commoneconomic.md"
                ],
        "Productivity Change Models" => Any[
                "Malmquist Index" => "productivity/malmquist.md"
                "Common functions for productivity change models" => "productivity/commonproductivity.md"
                ],
        "Statistical Analysis" => Any[
                "Bootstrap Radial DEA Model" => "statistical/radialboot.md",
                "Returns to Scale Test" => "statistical/rtstest.md"
                ],
        "Configuring Optimizer" => "optimizer.md",
        "Bibliography" => "bibliography.md"
        ],
        format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/javierbarbero/DataEnvelopmentAnalysis.jl.git",
)
