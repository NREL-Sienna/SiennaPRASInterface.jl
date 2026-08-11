# [Inspect or export a PRAS SystemModel](@id inspect_pras_system)

```@meta
CurrentModule = SiennaPRASInterface
```

[`generate_pras_system`](@ref) translates a [`PowerSystems.System`](@extref) into a
[`PRASCore.Systems.SystemModel`](@extref) without running a simulation. Use it to verify mappings,
share a `.pras` file with native PRAS workflows, or debug translation issues before calling
[`assess`](@ref).

## Build a system model

```@example inspect_pras
using SiennaPRASInterface
using PowerSystems
import PowerSystemCaseBuilder # hide
const PSY = PowerSystems # hide
const PSCB = PowerSystemCaseBuilder # hide
function PSY.get_storage_capacity(res::PSY.HydroReservoir) # hide
    return PSY.get_storage_level_limits(res).max # hide
end # hide
sys = PSCB.build_system(PSCB.SPISystems, "RTS_GMLC_Hourly with Static Outage Data") # hide
PSY.set_units_base_system!(sys, PSY.UnitSystem.NATURAL_UNITS) # hide

pras_system = generate_pras_system(sys, PSY.Area)
pras_system
```

The printed summary lists PRAS regions, generators, storage, and interfaces derived from the
Sienna system. For field-level details, see the
[PRAS System Model Specification](@extref PRASCore :doc:`PRAS/sysmodelspec`)
and [`.pras` file format](@extref PRASCore :doc:`SystemModel_HDF5_spec`).

## Export to a `.pras` file

Pass `export_location` with a path ending in `.pras`:

```julia
generate_pras_system(sys, PSY.Area; export_location="rts_gmlc.pras")
```

You can also pass `export_location` to the [`RATemplate`](@ref) overload:

```julia
template = RATemplate(PSY.Area, deepcopy(SiennaPRASInterface.DEFAULT_DEVICE_MODELS))
generate_pras_system(sys, template, "rts_gmlc.pras")
```

Exported files can be read with [PRASFiles](@extref PRASCore :doc:`PRASFiles/api`)
and used in the [PRAS 101 Walkthrough](@extref PRASCore :doc:`examples/pras_walkthrough`).

## Requirements

Translation requires:

  - At least one static time series on the Sienna system
  - At least one component of the chosen aggregation topology (for example [`PowerSystems.Area`](@extref))
  - Outage data on devices, tabulated defaults, or a template with `add_default_transition_probabilities=true`

See [Prepare stochastic outage data for resource adequacy](@ref prepare_stochastic_outage_data)
and [Default outage values](@ref default_outage_values).

## See also

  - [Configure device mappings for PRAS](@ref configure_device_mappings)
  - [Resource adequacy workflow](@ref resource_adequacy_workflow) — tutorial section on inspecting the translated model
