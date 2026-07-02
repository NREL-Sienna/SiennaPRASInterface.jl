# [Default outage values](@id default_outage_values)

```@meta
CurrentModule = SiennaPRASInterface
```

When a [`PowerSystems.System`](@extref) does not contain any
[`PowerSystems.GeometricDistributionForcedOutage`](@extref) supplemental attributes attached to
components, [`generate_pras_system`](@ref) applies a set of default outage rates defined in the
[Default Outage Rates CSV](https://github.com/Sienna-Platform/SiennaPRASInterface.jl/blob/main/src/util/descriptors/outage-rates-ERCOT-modified.csv),
based on rates from ERCOT.

For any remaining components not captured by the CSV defaults — such as lines, renewables, and
storage — the outage rates are zero and those devices will not fail stochastically.

You can also request default injection for specific device types when building a
[`RATemplate`](@ref) by setting `add_default_transition_probabilities=true` on formulations such
as [`GeneratorPRAS`](@ref). That path adds
[`PowerSystems.GeometricDistributionForcedOutage`](@extref) attributes only where they are
missing, using the same tabulated values.

## See also

  - [Prepare stochastic outage data for resource adequacy](@ref prepare_stochastic_outage_data)
  - [Outage and contingency data](@extref PowerSystems :label:`outage_and_contingency_data`)
