# [Prepare stochastic outage data for resource adequacy](@id prepare_stochastic_outage_data)

```@meta
CurrentModule = SiennaPRASInterface
DocTestSetup = quote
    using PowerSystems
    using PowerSystemCaseBuilder
end
```

Resource adequacy studies in PRAS model random equipment failures and recoveries. On a Sienna
system, attach a [`PowerSystems.GeometricDistributionForcedOutage`](@extref) supplemental
attribute to each device that should participate in stochastic outages.

This differs from the deterministic and planned outage patterns in the PowerSystems
[Model Outages](@extref PowerSystems :label:`model_outages`)
how-to, which targets security-constrained studies with [`PowerSystems.FixedForcedOutage`](@extref)
and [`PowerSystems.PlannedOutage`](@extref). For background on outage types, see
[Outage and contingency data](@extref PowerSystems :label:`outage_and_contingency_data`).

## Attach fixed outage rates

Use [`PowerSystems.add_supplemental_attribute!`](@extref PowerSystems :jl:method:`PowerSystems.add_supplemental_attribute!-Tuple{System, Component, SupplementalAttribute}`) following the general supplemental
attribute pattern in
[Attach supplemental data to components](@extref PowerSystems :label:`attach_contextual_data`).

```@example prepare_outage_data
using PowerSystems
import PowerSystemCaseBuilder
const PSCB = PowerSystemCaseBuilder
const PSY = PowerSystems
sys = PSCB.build_system(PSCB.SPISystems, "RTS_GMLC_Hourly with Static Outage Data")
PSY.set_units_base_system!(sys, PSY.UnitSystem.NATURAL_UNITS)

transition_data = GeometricDistributionForcedOutage(;
    mean_time_to_recovery=10.0,
    outage_transition_probability=0.005,
)
component = get_component(ThermalStandard, sys, "101_CT_1")
add_supplemental_attribute!(sys, component, transition_data)
transition_data
```

`mean_time_to_recovery` is in hours. `outage_transition_probability` is the per-hour probability
of transitioning into an outage state.

## Time-varying outage rates

To vary failure and recovery rates over the study horizon, attach time series named
`outage_probability` and `recovery_probability` to the supplemental attribute. The interface
reads these names when building the PRAS model. For the general time-series workflow, see
[Working with Time Series Data](@extref PowerSystems :doc:`tutorials/generated_working_with_time_series`).

## When no outage data is present

If a system has no [`PowerSystems.GeometricDistributionForcedOutage`](@extref) attributes,
[`generate_pras_system`](@ref) can apply tabulated default rates or add defaults through
[`GeneratorPRAS`](@ref) with `add_default_transition_probabilities=true`. See
[Default outage values](@ref default_outage_values) for how those defaults are chosen.

## See also

  - [Resource adequacy workflow](@ref resource_adequacy_workflow) — end-to-end assessment tutorial
  - [Configure device mappings for PRAS](@ref configure_device_mappings) — enable default outage injection per device type
  - [PRAS System Model Specification](@extref PRASCore :doc:`PRAS/sysmodelspec`) — how outage statistics appear in PRAS
