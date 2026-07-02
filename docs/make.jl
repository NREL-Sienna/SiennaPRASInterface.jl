using Documenter
import DataStructures: OrderedDict
using DocumenterInterLinks
using PowerSystems
using SiennaPRASInterface
using SiennaPRASInterface.PRASCore

links = InterLinks(
    "PowerSystems" => "https://sienna-platform.github.io/PowerSystems.jl/stable/",
    "InfrastructureSystems" => "https://sienna-platform.github.io/InfrastructureSystems.jl/stable/",
    "PRASCore" => "https://natlabrockies.github.io/PRAS/stable/",
    "PowerSystemCaseBuilder" => "https://sienna-platform.github.io/PowerSystemCaseBuilder.jl/stable/",
)

include(joinpath(@__DIR__, "make_tutorials.jl"))
make_tutorials()

pages = OrderedDict(
    "Welcome Page" => "index.md",
    "Tutorials" => Any[
        "Resource adequacy workflow" => "tutorials/generated_resource_adequacy_workflow.md",
    ],
    "How to..." => Any[
        "Prepare stochastic outage data for resource adequacy" => "how_to_guides/how_do_i_add_outage_data.md",
        "Configure device mappings for PRAS" => "how_to_guides/how_to_configure_device_mappings.md",
        "Request multiple result specifications" => "how_to_guides/how_to_request_multiple_results.md",
        "Inspect or export a PRAS SystemModel" => "how_to_guides/how_to_inspect_pras_system.md",
    ],
    "Explanation" => Any[
        "Default outage values" => "explanations/default_outage_values.md",
    ],
    "Reference" => Any[
        "Public API Reference" => "api/public.md",
        "Internal API Reference" => "api/internal.md",
    ],
)

makedocs(;
    modules=[SiennaPRASInterface, PRASCore],
    format=Documenter.HTML(;
        prettyurls=haskey(ENV, "GITHUB_ACTIONS"),
        size_threshold=nothing,
    ),
    sitename="SiennaPRASInterface.jl",
    authors="Surya Dhulipala, Joseph McKinsey, José Daniel Lara",
    pages=Any[p for p in pages],
    checkdocs=:none,
    doctest=false,
    plugins=[links],
)

deploydocs(;
    repo="github.com/Sienna-Platform/SiennaPRASInterface.jl.git",
    target="build",
    branch="gh-pages",
    devbranch="main",
    devurl="dev",
    push_preview=true,
    versions=["stable" => "v^", "v#.#"],
)
