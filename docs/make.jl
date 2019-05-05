using Documenter, PressureFieldContact

makedocs(
    modules = [PressureFieldContact],
    root = @__DIR__,
    checkdocs = :exports,
    sitename ="PressureFieldContact.jl",
    authors = "Ryan Elandt and contributors.",
    pages = [
        "Home" => "index.md"
        "Algorithms" => "algorithms.md"
      ],
    format = Documenter.HTML(prettyurls = parse(Bool, get(ENV, "CI", "false")))
)

deploydocs(
    repo = "github.com/ryanelandt/PressureFieldContact.jl.git"
)
