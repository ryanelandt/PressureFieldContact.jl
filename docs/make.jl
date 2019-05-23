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
        "Add Contact/Friction" => "friction.md"
        "MechanismScenario" => "mechanism_scenario.md"
        "Algorithms" => "algorithms.md"
        "Polygon Clipping" => "polygon_clipping.md"
        "Examples" => "examples.md"
      ],
    format = Documenter.HTML(prettyurls = parse(Bool, get(ENV, "CI", "false")))
)

deploydocs(
    repo = "github.com/ryanelandt/PressureFieldContact.jl.git"
)
