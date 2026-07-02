# [Request multiple result specifications](@id request_multiple_results)

```@meta
CurrentModule = SiennaPRASInterface
```

[`assess`](@ref) accepts any number of PRAS
[result specifications](@extref PRASCore :doc:`PRAS/results`) and returns
one result object per specification, in order.

## Example

```@example request_results
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

method = SequentialMonteCarlo(; samples=20, seed=1, threaded=false, verbose=false)
shortfall, surplus, storage_energy = assess(
    sys,
    PSY.Area,
    method,
    Shortfall(),
    Surplus(),
    StorageEnergy(),
)
typeof(shortfall), typeof(surplus), typeof(storage_energy)
```

## Common specifications

| Specification                                       | Purpose                                                                                                               |
|:--------------------------------------------------- |:--------------------------------------------------------------------------------------------------------------------- |
| [`PRASCore.Results.Shortfall`](@extref)             | Unserved load by region and sample; use with [`PRASCore.Results.LOLE`](@extref) and [`PRASCore.Results.EUE`](@extref) |
| [`PRASCore.Results.Surplus`](@extref)               | Excess generation periods                                                                                             |
| [`PRASCore.Results.Flow`](@extref)                  | Interface flow utilization                                                                                            |
| [`PRASCore.Results.Utilization`](@extref)           | Generator output utilization                                                                                          |
| [`PRASCore.Results.StorageEnergy`](@extref)         | Storage energy state statistics                                                                                       |
| [`PRASCore.Results.GeneratorAvailability`](@extref) | Generator availability time series per sample                                                                         |

Append `Samples` to a specification type (for example [`PRASCore.Results.ShortfallSamples`](@extref))
to retain per-sample arrays instead of summarized regional results. See the
[PRAS results documentation](@extref PRASCore :doc:`PRAS/results`) for
full semantics.

## See also

  - [Resource adequacy workflow](@ref resource_adequacy_workflow) — tutorial steps requesting multiple specifications
  - [PRASCore API — Results](@extref PRASCore :label:`Results`)
