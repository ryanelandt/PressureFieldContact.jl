using Documenter, SoftContact

makedocs(
    modules = [SoftContact],
    root = @__DIR__,
    checkdocs = :exports,
    sitename ="SoftContact.jl",
    authors = "Ryan Elandt and contributors.",
    pages = [
        "Home" => "index.md"
        "Algorithms" => "algorithms.md"
      ],
    format = Documenter.HTML(prettyurls = parse(Bool, get(ENV, "CI", "false")))
)

deploydocs(
    repo = "github.com/ryanelandt/SoftContact.jl.git"
)
