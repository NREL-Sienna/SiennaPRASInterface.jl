# [Configure device mappings for PRAS](@id configure_device_mappings)

```@meta
CurrentModule = SiennaPRASInterface
```

[`RATemplate`](@ref) and [`DeviceRAModel`](@ref) control how [PowerSystems.jl](https://sienna-platform.github.io/PowerSystems.jl/stable/)
components translate into a [`PRASCore.Systems.SystemModel`](@extref).
Use them when the default mappings do not match your study assumptions — for example, when
lumping renewable generation within a region or selecting different time-series field names.

## Default behavior

Calling [`assess`](@ref) without a template uses area-level aggregation and built-in device
mappings for generators, storage, lines, and loads. To inspect the default translation, see
[Inspect or export a PRAS SystemModel](@ref inspect_pras_system).

## Build a custom template

```@example configure_mappings
using SiennaPRASInterface
using PowerSystems
import PowerSystemCaseBuilder # hide
const PSY = PowerSystems # hide
const PSCB = PowerSystemCaseBuilder # hide
sys = PSCB.build_system(PSCB.SPISystems, "RTS_GMLC_Hourly with Static Outage Data") # hide
PSY.set_units_base_system!(sys, PSY.UnitSystem.NATURAL_UNITS) # hide

template = RATemplate(
    PSY.Area,
    deepcopy(SiennaPRASInterface.DEFAULT_DEVICE_MODELS),
)
set_device_model!(
    template,
    DeviceRAModel(
        PSY.RenewableDispatch,
        GeneratorPRAS;
        lump_renewable_generation=true,
    ),
)
template
```

[`set_device_model!`](@ref) appends a mapping. Later mappings take precedence when multiple
models apply to the same device type.

## Formulation types

Each [`DeviceRAModel`](@ref) pairs a Sienna device supertype with a formulation struct:

| Formulation                        | Role                                   |
|:---------------------------------- |:-------------------------------------- |
| [`GeneratorPRAS`](@ref)            | Thermal and renewable generators       |
| [`HybridSystemPRAS`](@ref)         | Hybrid systems with co-located storage |
| [`HydroEnergyReservoirPRAS`](@ref) | Hydro reservoirs with energy limits    |
| [`EnergyReservoirSoC`](@ref)       | Storage state-of-charge devices        |
| [`LinePRAS`](@ref)                 | Transmission branches                  |
| [`AreaInterchangeLimit`](@ref)     | Inter-area transfer limits             |
| [`StaticLoadPRAS`](@ref)           | Static loads                           |

[`GeneratorPRAS`](@ref) accepts keyword arguments for time-series names (for example
`max_active_power="max_active_power"`) and `add_default_transition_probabilities` to inject
[`PowerSystems.GeometricDistributionForcedOutage`](@extref) data when none is attached.

## Run an assessment with the template

```julia
method = SequentialMonteCarlo(; samples=100, seed=1, threaded=false)
shortfalls, = assess(sys, template, method, Shortfall())
LOLE(shortfalls)
```

## See also

  - [Resource adequacy workflow](@ref resource_adequacy_workflow) — tutorial section on [`RATemplate`](@ref)
  - [Prepare stochastic outage data for resource adequacy](@ref prepare_stochastic_outage_data)
  - [PRAS System Model Specification](@extref PRASCore :doc:`PRAS/sysmodelspec`)
