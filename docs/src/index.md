# Welcome to SiennaPRASInterface.jl

```@meta
CurrentModule = SiennaPRASInterface
```

## About

`SiennaPRASInterface.jl` is part of the National Laboratory of the Rockies'
[Sienna ecosystem](https://www.nlr.gov/analysis/sienna.html), an open source framework for
scheduling problems and dynamic simulations for power systems.

`SiennaPRASInterface.jl` is a [`Julia`](http://www.julialang.org) package that bridges
[`PowerSystems.jl`](https://github.com/Sienna-Platform/PowerSystems.jl) and the Probabilistic
Resource Adequacy Suite ([`PRAS.jl`](https://natlabrockies.github.io/PRAS/stable/)) from
[NatLabRockies](https://github.com/NatLabRockies/PRAS). It translates a
[`PowerSystems.System`](@extref) into a [`PRASCore.Systems.SystemModel`](@extref) and runs Monte Carlo
resource adequacy studies through [`assess`](@ref).

The package handles:

  - Mapping Sienna components, time series, and outage data into PRAS device and region models
  - Running [`PRASCore.Simulations.SequentialMonteCarlo`](@extref) simulations via a single
    [`assess`](@ref) entry point
  - Computing reliability metrics such as [`PRASCore.Results.LOLE`](@extref) and
    [`PRASCore.Results.EUE`](@extref) from shortfall results

Start with the [Resource adequacy workflow](@ref resource_adequacy_workflow) tutorial for an end-to-end example.

## How to use this documentation

There are four main sections containing different information:

  - **Tutorials** — Detailed walk-throughs to help you *learn* how to run a resource adequacy
    study with `SiennaPRASInterface.jl`
  - **How to...** — Directions to help *guide* your work for a particular task that bridges
    `PowerSystems.jl` and PRAS
  - **Explanation** — Additional details and background information to help you *understand*
    default behaviors and design choices
  - **Reference** — Technical references and API for a quick *look-up* during your work,
    including re-exported PRAS types and functions

`SiennaPRASInterface.jl` strives to follow the [Diátaxis](https://diataxis.fr/) documentation
framework.

## Installation and quick links

  - Install from the Julia package manager:

    ```julia
    using Pkg;
    Pkg.add("SiennaPRASInterface")
    ```

  - [Sienna installation page](https://sienna-platform.github.io/Sienna/SiennaDocs/docs/build/how-to/install/):
    Instructions to install `SiennaPRASInterface.jl` and other Sienna packages

  - [Sienna Documentation Hub](https://sienna-platform.github.io/Sienna/SiennaDocs/docs/build/index.html):
    Links to other Sienna packages' documentation

  - [PRAS documentation](https://natlabrockies.github.io/PRAS/stable/):
    Upstream simulation methods, result specifications, and system model details

  - [PowerSystems outage and contingency data](@extref PowerSystems :label:`outage_and_contingency_data`):
    Background on stochastic outage supplemental attributes used by this interface

!!! note

    `SiennaPRASInterface.jl` depends on [`PowerSystems.jl`](https://sienna-platform.github.io/PowerSystems.jl/stable/)
    for the system data model and on [`PRASCore`](@extref PRASCore :doc:`PRASCore/api`)
    for simulation and results. For most workflows you import `SiennaPRASInterface` and
    `PowerSystems`; PRAS types are re-exported for convenience.

* * *

SiennaPRASInterface has been developed as part of the Transmission Planning Tools Maintenance
project at the U.S. Department of Energy's National Renewable Energy Laboratory
([NREL](https://www.nrel.gov/)) funded by DOE Grid Deployment Office (GDO).
