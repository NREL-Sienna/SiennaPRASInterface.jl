```@meta
CurrentModule = SiennaPRASInterface
DocTestSetup = quote
    using SiennaPRASInterface
end
```

# Public API Reference

`SiennaPRASInterface.jl` exposes translation and assessment methods plus re-exported
[`PRASCore`](@extref PRASCore :doc:`PRASCore/api`) simulation types. For
PRAS theory, simulation options, and result semantics, see the
[PRAS documentation](https://natlabrockies.github.io/PRAS/stable/).

```@contents
Pages = ["public.md"]
Depth = 3
```

```@raw html
&nbsp;
&nbsp;
```

## Sienna PRAS interface

### Assessment entry points

```@docs
assess
```

### System translation

```@docs
generate_pras_system
```

#### Templates and device mappings

```@autodocs
Modules = [SiennaPRASInterface]
Pages = ["formulation_definitions.jl"]
Order = [:type, :function]
Public = true
Private = false
```

### Outage workflows

```@docs
generate_outage_profile!
make_generator_outage_draws!
```

```@raw html
&nbsp;
&nbsp;
```

* * *

## PRAS simulation API (re-exported from PRASCore)

The symbols below are imported from [`PRASCore`](@extref PRASCore :doc:`PRASCore/api`)
and re-exported by `SiennaPRASInterface.jl`. Native PRAS methods such as
[`PRASCore.Simulations.assess`](@extref) operate on a
[`PRASCore.Systems.SystemModel`](@extref); use the
Sienna [`assess`](@ref) overloads to translate from a [`PowerSystems.System`](@extref) automatically.

### System model

See the [PRAS System Model Specification](@extref PRASCore :doc:`PRAS/sysmodelspec`)
and [`PRASCore` API](@extref PRASCore :label:`Systems`) for
`SystemModel` fields and component groups.

### Simulation methods

```@docs
PRASCore.Simulations.SequentialMonteCarlo
```

### Result specifications

```@docs
PRASCore.Results.Shortfall
PRASCore.Results.ShortfallSamples
PRASCore.Results.Surplus
PRASCore.Results.SurplusSamples
PRASCore.Results.Flow
PRASCore.Results.FlowSamples
PRASCore.Results.Utilization
PRASCore.Results.UtilizationSamples
PRASCore.Results.StorageEnergy
PRASCore.Results.StorageEnergySamples
PRASCore.Results.GeneratorStorageEnergy
PRASCore.Results.GeneratorStorageEnergySamples
```

### Availability results

| Symbol | Description |
|:-------|:------------|
| [`PRASCore.Results.GeneratorAvailability`](@extref) | Generator availability time series per sample |
| [`PRASCore.Results.StorageAvailability`](@extref) | Storage availability time series per sample |
| [`PRASCore.Results.GeneratorStorageAvailability`](@extref) | Generator-storage availability per sample |
| [`PRASCore.Results.LineAvailability`](@extref) | Line availability time series per sample |

### Reliability metrics

```@docs
PRASCore.Results.LOLE
PRASCore.Results.EUE
```

Use `val` and `stderror` on reliability metric objects returned by `LOLE` and `EUE`.
