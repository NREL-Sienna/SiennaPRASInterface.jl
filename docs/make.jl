using Documenter, PRASInterface
import OrderedCollections: OrderedDict

pages = OrderedDict(
    "Welcome Page" => "index.md",
    "Quick Start Guide" => "quick_start_guide.md",
    "Tutorials" => "tutorials/intro_page.md",
    "Public API Reference" => "api/public.md",
    "Internal API Reference" => "api/internal.md",
)

makedocs(
    modules=[PRASInterface],
    format=Documenter.HTML(prettyurls=haskey(ENV, "GITHUB_ACTIONS")),
    sitename="PRASInterface.jl",
    authors="Surya Dhulipala, Joseph McKinsey, José Daniel Lara",
    pages=Any[p for p in pages],
)

deploydocs(
    repo="github.com/NREL-SIENNA/PRASInterface.git",
    target="build",
    branch="gh-pages",
    devbranch="main",
    devurl="dev",
    push_preview=true,
    versions=["stable" => "v^", "v#.#"],
)
