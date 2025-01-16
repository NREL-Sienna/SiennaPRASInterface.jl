using Documenter, SiennaPRASInterface
import SiennaPRASInterface: PRASCore
import OrderedCollections: OrderedDict

pages = OrderedDict(
    "Welcome Page" => "index.md",
    "Tutorials" =>
        Any["How do I add outage data?" => "tutorials/how_do_i_add_outage_data.md"],
    "Public API Reference" => "api/public.md",
    "Internal API Reference" => "api/internal.md",
)

makedocs(
    modules=[SiennaPRASInterface, PRASCore],
    format=Documenter.HTML(prettyurls=haskey(ENV, "GITHUB_ACTIONS")),
    sitename="SiennaPRASInterface.jl",
    authors="Surya Dhulipala, Joseph McKinsey, José Daniel Lara",
    pages=Any[p for p in pages],
)

deploydocs(
    repo="github.com/NREL-Sienna/SiennaPRASInterface.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="main",
    devurl="dev",
    push_preview=true,
    versions=["stable" => "v^", "v#.#"],
)
