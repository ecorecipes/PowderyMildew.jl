using Documenter
using PowderyMildew

makedocs(;
         modules = [PowderyMildew],
         authors = "Simon Frost",
         sitename = "PowderyMildew.jl",
         checkdocs = :exports,
         format = Documenter.HTML(;
                                  prettyurls = get(ENV, "CI", "false") == "true",
                                  canonical = "https://ecorecipes.github.io/PowderyMildew.jl",
                                  size_threshold = 300 * 1024),
         pages = ["Home" => "index.md",
                  "Tutorials" => ["The model of Hamelin et al. (2016)" => "tutorials/01_hamelin.md",
                                  "Timing and severity" => "tutorials/02_leaf_model.md",
                                  "Spread" => "tutorials/03_spread.md"],
                  "API Reference" => ["Hamelin et al. (2016)" => "api/hamelin.md",
                                      "Leaf-level model" => "api/leaf_model.md",
                                      "Spatial model" => "api/spatial.md"]])

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo = "github.com/ecorecipes/PowderyMildew.jl.git", push_preview = true)
end
