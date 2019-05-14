using Documenter, PressureFieldContact

makedocs(
    modules = [PressureFieldContact],
    root = @__DIR__,
    checkdocs = :exports,
    sitename ="PressureFieldContact.jl",
    authors = "Ryan Elandt and contributors.",
    pages = [
        "Home" => "index.md"
        "Quickstart guide" => "quick_start.md"
        "Geometry" => "geometry.md"
        "Friction" => "friction.md"
        "MechanismScenario" => "mechanism_scenario.md"
        "Polygon clipping" => "polygon_clipping.md"
        "Algorithms" => "algorithms.md"
        "Examples" => "examples.md"
      ],
    format = Documenter.HTML(prettyurls = parse(Bool, get(ENV, "CI", "false")))
)

deploydocs(
    repo = "github.com/ryanelandt/PressureFieldContact.jl.git"
)
