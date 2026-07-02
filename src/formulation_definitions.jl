"""
    AbstractRAFormulation

Abstract type for translating a Sienna object in PRAS. Multiple objects can
use the same formulation.

Formulations are also intended to contain information about their
configuration such as time series names.

PRAS formulation subtypes for specific PRAS types

  - [`GeneratorPRAS`](@ref)
  - [`StoragePRAS`](@ref) (abstract type)
  - [`GeneratorStoragePRAS`](@ref) (abstract type)
"""
abstract type AbstractRAFormulation end

"""
    GeneratorPRAS(;
        max_active_power = "max_active_power",
        lump_renewable_generation = false,
        add_default_transition_probabilities = false,
        outage_probability = "outage_probability",
        recovery_probability = "recovery_probability",
    )

GeneratorPRAS maps [`PowerSystems.Generator`](@extref) devices to PRAS generator entries.

# Arguments

$(TYPEDFIELDS)
"""
struct GeneratorPRAS <: AbstractRAFormulation
    """
    Name of time series to use for max active power
    """
    max_active_power::String
    """
    Whether to lump renewable generation to regions
    """
    lump_renewable_generation::Bool
    """
    Whether to add default outage data to generators
    """
    add_default_transition_probabilities::Bool
    """
    Name of time series to use for outage_probability
    """
    outage_probability::String
    """
    Name of time series to use for recovery_probability
    """
    recovery_probability::String

    function GeneratorPRAS(;
        max_active_power="max_active_power",
        lump_renewable_generation=false,
        add_default_transition_probabilities=false,
        outage_probability="outage_probability",
        recovery_probability="recovery_probability",
    )
        return new(
            max_active_power,
            lump_renewable_generation,
            add_default_transition_probabilities,
            outage_probability,
            recovery_probability,
        )
    end
end

"""
Get max active power time series name
"""
function get_max_active_power(f::GeneratorPRAS)
    return f.max_active_power
end

"""
Get whether renewable generation is lumped to regions
"""
function get_lump_renewable_generation(f::GeneratorPRAS)
    return f.lump_renewable_generation
end

"""
Get whether default forced outages needed to be added to generators
"""
function get_add_default_transition_probabilities(f::GeneratorPRAS)
    return f.add_default_transition_probabilities
end

"""
Get outage_probability time series name
"""
function get_outage_probability(f::GeneratorPRAS)
    return f.outage_probability
end

"""
Get recovery_probability time series name
"""
function get_recovery_probability(f::GeneratorPRAS)
    return f.recovery_probability
end

"""
    GeneratorStoragePRAS <: AbstractRAFormulation

Objects in Sienna that behave like generator and storage are mapped to generator-storage entries in [`PRASCore.Systems.SystemModel`](@extref).

To add a generator storage formulation, you must also add a [`assign_to_gen_stor_matrices!`](@ref) function.

  - [`HybridSystemPRAS`](@ref)
  - [`HydroEnergyReservoirPRAS`](@ref)
"""
abstract type GeneratorStoragePRAS <: AbstractRAFormulation end
"""
    HybridSystemPRAS(;
        max_active_power = "max_active_power",
        add_default_transition_probabilities = false,
        outage_probability = "outage_probability",
        recovery_probability = "recovery_probability",
    )

HybridSystemPRAS maps hybrid systems to PRAS generator-storage entries.

# Arguments

$(TYPEDFIELDS)
"""
struct HybridSystemPRAS <: GeneratorStoragePRAS
    """
    Name of time series to use for max active power
    """
    max_active_power::String
    """
    Whether to add default outage data
    """
    add_default_transition_probabilities::Bool
    """
    Name of time series to use for outage_probability
    """
    outage_probability::String
    """
    Name of time series to use for recovery_probability
    """
    recovery_probability::String

    function HybridSystemPRAS(;
        max_active_power="max_active_power",
        add_default_transition_probabilities=false,
        outage_probability="outage_probability",
        recovery_probability="recovery_probability",
    )
        return new(
            max_active_power,
            add_default_transition_probabilities,
            outage_probability,
            recovery_probability,
        )
    end
end

"""
    HydroEnergyReservoirPRAS(;
        max_active_power = "max_active_power",
        inflow = "inflow",
        storage_capacity = "storage_capacity",
        add_default_transition_probabilities = false,
        outage_probability = "outage_probability",
        recovery_probability = "recovery_probability",
    )

Maps hydro energy reservoirs to PRAS generator-storage entries.

# Arguments

$(TYPEDFIELDS)
"""
struct HydroEnergyReservoirPRAS <: GeneratorStoragePRAS
    """
    Name of time series to use for max active power
    """
    max_active_power::String
    """
    Name of time series to use for inflow
    """
    inflow::String
    """
    Name of time series to use for storage capacity
    """
    storage_capacity::String
    """
    Whether to add default outage data
    """
    add_default_transition_probabilities::Bool
    """
    Name of time series to use for outage_probability
    """
    outage_probability::String
    """
    Name of time series to use for recovery_probability
    """
    recovery_probability::String

    function HydroEnergyReservoirPRAS(;
        max_active_power="max_active_power",
        inflow="inflow",
        storage_capacity="storage_capacity",
        add_default_transition_probabilities=false,
        outage_probability="outage_probability",
        recovery_probability="recovery_probability",
    )
        return new(
            max_active_power,
            inflow,
            storage_capacity,
            add_default_transition_probabilities,
            outage_probability,
            recovery_probability,
        )
    end
end

"""
Get max active power time series name
"""
function get_max_active_power(f::GeneratorStoragePRAS)
    return f.max_active_power
end

"""
Get whether default forced outages needed to be added to generatorstorages
"""
function get_add_default_transition_probabilities(f::GeneratorStoragePRAS)
    return f.add_default_transition_probabilities
end

"""
Get outage_probability time series name
"""
function get_outage_probability(f::GeneratorStoragePRAS)
    return f.outage_probability
end

"""
Get recovery_probability time series name
"""
function get_recovery_probability(f::GeneratorStoragePRAS)
    return f.recovery_probability
end

"""
Get inflow time series name
"""
function get_inflow(f::HydroEnergyReservoirPRAS)
    return f.inflow
end

"""
Get storage capacity time series name
"""
function get_storage_capacity(f::HydroEnergyReservoirPRAS)
    return f.storage_capacity
end

"""
    StoragePRAS <: AbstractRAFormulation

Objects in Sienna that behave like storage are mapped to storage entries in [`PRASCore.Systems.SystemModel`](@extref).

Subtypes must provide [`assign_to_stor_matrices!`](@ref) function.
"""
abstract type StoragePRAS <: AbstractRAFormulation end

"""
    EnergyReservoirSoC(;
        add_default_transition_probabilities = false,
        outage_probability = "outage_probability",
        recovery_probability = "recovery_probability",
    )

Storage formulation that tracks state of charge for energy reservoir devices.

# Arguments

$(TYPEDFIELDS)
"""
struct EnergyReservoirSoC <: StoragePRAS
    """
    Whether to add default outage data
    """
    add_default_transition_probabilities::Bool
    """
    Name of time series to use for outage_probability
    """
    outage_probability::String
    """
    Name of time series to use for recovery_probability
    """
    recovery_probability::String

    function EnergyReservoirSoC(;
        add_default_transition_probabilities=false,
        outage_probability="outage_probability",
        recovery_probability="recovery_probability",
    )
        return new(
            add_default_transition_probabilities,
            outage_probability,
            recovery_probability,
        )
    end
end

"""
Get whether default forced outages needed to be added to generatorstorages
"""
function get_add_default_transition_probabilities(f::StoragePRAS)
    return f.add_default_transition_probabilities
end

"""
Get outage_probability time series name
"""
function get_outage_probability(f::StoragePRAS)
    return f.outage_probability
end

"""
Get recovery_probability time series name
"""
function get_recovery_probability(f::StoragePRAS)
    return f.recovery_probability
end

"""
    InterfacePRAS <: AbstractRAFormulation

InterfacePRAS produces interface entries in PRAS.

If not supplied, then interfaces will be generated by Lines

Subtypes must provide [`add_to_interface!`](@ref) function.
"""
abstract type InterfacePRAS <: AbstractRAFormulation end

"""
    AreaInterchangeLimit <: InterfacePRAS

AreaInterchangeLimit produces interfaces from [`PowerSystems.AreaInterchange`](@extref) objects.

Each line must have a corresponding AreaInterchange. All AreaInterchange
objects will be consolidated for each pair of directly connected regions.
"""
struct AreaInterchangeLimit <: InterfacePRAS end

"""
    LinePRAS <: AbstractRAFormulation

LinePRAS produces line entries in PRAS.

See [`assign_to_line_matrices!`](@ref) for the formulation handling. Any
subtypes must implement this.
"""
struct LinePRAS <: AbstractRAFormulation end

"""
    LoadPRAS <: AbstractRAFormulation

See [`add_to_load_matrix!`](@ref) for how the formulation is used to add load to
regions.
"""
abstract type LoadPRAS <: AbstractRAFormulation end

"""
    StaticLoadPRAS(; max_active_power = "max_active_power")

Maps static loads to PRAS regional load entries.

# Arguments

$(TYPEDFIELDS)
"""
struct StaticLoadPRAS <: LoadPRAS
    """
    Name of time series to use for max active power
    """
    max_active_power::String

    function StaticLoadPRAS(; max_active_power="max_active_power")
        return new(max_active_power)
    end
end

"""
    DeviceRAModel{D <: PSY.Device, B <: AbstractRAFormulation}

# Arguments

  - `D <: `[`PowerSystems.Device`](@extref): device type
    $(TYPEDFIELDS)

A DeviceRAModel, like a DeviceModel in PowerSimulations, assigns a [`PowerSystems.Device`](@extref)
to a specific formulation. Unlike Sienna, we put configuration information
in the formulation itself.
"""
struct DeviceRAModel{D <: PSY.Device, B <: AbstractRAFormulation}
    """
    Formulation containing configuration
    """
    formulation::B

    function DeviceRAModel(
        ::Type{D},
        formulation::B,
    ) where {D <: PSY.Device, B <: AbstractRAFormulation}
        return new{D, B}(formulation)
    end
end

"""
Get formulation from a [`DeviceRAModel`](@ref)
"""
function get_formulation(f::DeviceRAModel)
    return f.formulation
end

"""
    $(TYPEDSIGNATURES)

# Arguments

  - `::Type{D}`: [`PowerSystems.Device`](@extref) type
  - `::Type{B}`: [`AbstractRAFormulation`](@ref) type
  - `time_series_names::Dict{Symbol, String}`: Mapping of time series `Symbol` to names
  - `kwargs...`: Additional arguments to pass to the formulation constructor

Keyword arguments in [`DeviceRAModel`](@ref) are passed to the
formulation constructor.

You may also pass a `time_series_names` Dict to map time series `Symbol` to names.

# Example

```julia
DeviceRAModel(
    PSY.Generator,
    GeneratorPRAS(; max_active_power="max_active_power"),
)
```

```julia
DeviceRAModel(
    PSY.HydroEnergyReservoir,
    HydroEnergyReservoirPRAS;
    max_active_power="max_active_power",
    inflow="inflow",
    storage_capacity="storage_capacity",
)
```

```julia
DeviceRAModel(
    PSY.HybridSystem,
    HybridSystemPRAS;
    time_series_names=Dict(; :max_active_power="max_active_power"),
)
```
"""
function DeviceRAModel(
    ::Type{D},
    ::Type{B};
    time_series_names::Dict{Symbol, String}=Dict{Symbol, String}(),
    kwargs...,
) where {D <: PSY.Device, B <: AbstractRAFormulation}
    formulation = B(; time_series_names..., kwargs...)
    return DeviceRAModel(D, formulation)
end

"""
Check whether a [`DeviceRAModel`](@ref) applies to a given type
"""
function appliestodevice(::DeviceRAModel{D}, ::Type{T}) where {D, T}
    return T <: D
end

"""
Uses a PRAS device model to find all components matching it in a system.
"""
function get_available_components(
    ::DeviceRAModel{D},
    sys::PSY.System,
) where {D <: PSY.Device}
    return PSY.get_available_components(D, sys)
end

"""
    $(TYPEDSIGNATURES)

# Arguments

$(TYPEDFIELDS)

The RATemplate contains all configuration necessary for building
a PRAS simulation from a [`PowerSystems.System`](@extref).

PRAS models are processed in reverse order, with later models taking precedence.

# Example

```julia
template = RATemplate(
    PSY.Area,
    [
        DeviceRAModel(
            PSY.Generator,
            GeneratorPRAS(; max_active_power="max_active_power"),
        ),
        DeviceRAModel(
            PSY.HydroEnergyReservoir,
            HydroEnergyReservoirPRAS(;
                max_active_power="max_active_power",
                inflow="inflow",
                storage_capacity="storage_capacity",
            ),
        ),
    ],
)
```
"""
mutable struct RATemplate{T <: PSY.AggregationTopology}
    """
    [`PowerSystems.AggregationTopology`](@extref) level of aggregation to use for PRAS regions, since PRAS is an area-based model
    """
    aggregation::Type{T}
    """
    [`DeviceRAModel`](@ref)s to translate components into PRAS
    """
    device_models::Array{DeviceRAModel}

    function RATemplate(
        aggregation::Type{T}=PSY.Area,
        device_models::Vector{<:DeviceRAModel}=DeviceRAModel[],
    ) where {T <: PSY.AggregationTopology}
        return new{T}(aggregation, convert(Vector{DeviceRAModel}, device_models))
    end
end

"""
    $(TYPEDSIGNATURES)

# Arguments

  - `template::`[`RATemplate`](@ref): template to add device model to
  - `device_model::`[`DeviceRAModel`](@ref)`{D}`: device model to add, where `D` is a [`PowerSystems.Device`](@extref) type

Add a device model to a [`RATemplate`](@ref). If an existing model
already applies to the given device type, then a warning
is issued. However, newer models will take precedence.
"""
function set_device_model!(template::RATemplate, device_model::DeviceRAModel{D}) where {D}
    for existing_model in template.device_models
        if appliestodevice(existing_model, D)
            @warn "Device model $(D) already exists in template"
        end
    end
    push!(template.device_models, device_model)
end

"""
    $(TYPEDSIGNATURES)

# Arguments

  - `template::`[`RATemplate`](@ref): template to add device model to
  - `::Type{D}`: [`PowerSystems.Device`](@extref) type
  - `::Type{B}`: [`AbstractRAFormulation`](@ref) type

Adds a device model to a [`RATemplate`](@ref) by passing the type
to a constructor.
"""
function set_device_model!(
    template::RATemplate,
    ::Type{D},
    ::Type{B},
) where {D <: PSY.Device, B <: AbstractRAFormulation}
    set_device_model!(template, DeviceRAModel(D, B()))
end

"""
    $(SIGNATURES)

Constructs a dictionary from Sienna Devices to formulation objects
"""
function build_component_to_formulation(
    ::Type{B},
    sys::PSY.System,
    device_models::Array{DeviceRAModel},
)::Dict{PSY.Device, B} where {B <: AbstractRAFormulation}
    component_to_formulation = Dict{PSY.Device, B}()
    for device_model in reverse(device_models)
        if !(device_model.formulation isa B)
            continue
        end
        for component in get_available_components(device_model, sys)
            if haskey(component_to_formulation, component)
                @warn "Component $(PSY.get_name(component)) has multiple formulations. Choosing last applied"
                continue
            end
            component_to_formulation[component] = device_model.formulation
        end
    end
    return component_to_formulation
end

"""
    $(SIGNATURES)

Filter the dictionary from Sienna Devices to [`GeneratorPRAS`](@ref) formulation objects for Lumped vs. NonLumped
"""
function filter_component_to_formulation(gens_to_formula::Dict{PSY.Device, GeneratorPRAS})
    lumped_gens_to_formula = filter(
        ((k, v),) ->
            (
                get_lump_renewable_generation(v) &&
                PSY.has_supplemental_attributes(PSY.GeometricDistributionForcedOutage, k) &&
                all(
                    iszero.(
                        PSY.get_outage_transition_probability.(
                            PSY.get_supplemental_attributes(
                                PSY.GeometricDistributionForcedOutage,
                                k,
                            ),
                        ),
                    ),
                )
            ) || (
                get_lump_renewable_generation(v) &&
                !PSY.has_supplemental_attributes(PSY.GeometricDistributionForcedOutage, k)
            ),
        gens_to_formula,
    )
    nonlumped_gens_to_formula =
        filter(((k, v),) -> k ∉ keys(lumped_gens_to_formula), gens_to_formula)

    return lumped_gens_to_formula, nonlumped_gens_to_formula
end
