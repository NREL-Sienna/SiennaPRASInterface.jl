"""
    AbstractPRASFormulation

Abstract type for translating a Sienna object in PRAS. Multiple objects can
use the same formulation.

Formulations are also intended to contain information about their
configuration such as time series names.
"""
abstract type AbstractPRASFormulation end

"""
    GeneratorFormulation(; max_active_power, lump_renewable_generation) <: AbstractPRASFormulation

# Arguments
$(TYPEDFIELDS)

GeneratorFormulation produces generator entires in PRAS.
"""
struct GeneratorFormulation <: AbstractPRASFormulation
    "Name of time series to use for max active power"
    max_active_power::String
    "Whether to lump renewable generation to regions"
    lump_renewable_generation::Bool

    function GeneratorFormulation(;
        max_active_power="max_active_power",
        lump_renewable_generation=false,
    )
        return new(max_active_power, lump_renewable_generation)
    end
end

"""
    GeneratorStorageFormulation <: AbstractPRASFormulation

Objects in Sienna that behave like generator and storage are mapped to generatorstorage in PRAS.

Applies only to [`HybridSystemFormulation`](@ref)
"""
abstract type GeneratorStorageFormulation <: AbstractPRASFormulation end

"""
    HybridSystemFormulation(; max_active_power) <: GeneratorStorageFormulation

# Arguments
$(TYPEDFIELDS)

HybridSystemFormulation produces generatorstorage entries in PRAS.
"""
struct HybridSystemFormulation <: GeneratorStorageFormulation
    "Name of time series to use for max active power"
    max_active_power::String

    function HybridSystemFormulation(; max_active_power="max_active_power")
        return new(max_active_power)
    end
end

"""
    HydroEnergyReservoirFormulation <: GeneratorStorageFormulation

# Arguments
$(TYPEDFIELDS)
"""
struct HydroEnergyReservoirFormulation <: GeneratorStorageFormulation
    "Name of time series to use for max active power"
    max_active_power::String
    "Name of time series to use for inflow"
    inflow::String
    "Name of time series to use for storage capacity"
    storage_capacity::String

    function HydroEnergyReservoirFormulation(;
        max_active_power="max_active_power",
        inflow="inflow",
        storage_capacity="storage_capacity",
    )
        return new(max_active_power, inflow, storage_capacity)
    end
end

"""
    StorageFormulation <: AbstractPRASFormulation

Objects in Sienna that behave like storage are mapped to storage in PRAS.
"""
struct StorageFormulation <: AbstractPRASFormulation end

"""
    DevicePRASModel{D <: PSY.Device, B <: AbstractPRASFormulation}

# Arguments

- D <: PSY.Device: Device type
$(TYPEDFIELDS)

A DevicePRASModel, like a DeviceModel in PowerSimulations, assigns a type of Component
to a specific formulation. Unlike Sienna, we put configuration information
in the formulation itself.
"""
struct DevicePRASModel{D <: PSY.Device, B <: AbstractPRASFormulation}
    "Formulation containing configuration"
    formulation::B

    function DevicePRASModel(
        ::Type{D},
        formulation::B,
    ) where {D <: PSY.Device, B <: AbstractPRASFormulation}
        return new{D, B}(formulation)
    end
end

"""
    $(TYPEDSIGNATURES)
    
# Arguments
- `::Type{D}`: Device type
- `::Type{B}`: Formulation type
- `time_series_names::Dict{Symbol, String}`: Mapping of time series `Symbol` to names
- `kwargs...`: Additional arguments to pass to the formulation constructor

Keyword arguments in DevicePRASModel are passed to the
formulation constructor.

You may also pass a `time_series_names` Dict to map time series `Symbol` to names.

# Example

```julia
DevicePRASModel(
    PSY.Generator,
    GeneratorFormulation(max_active_power="max_active_power"),
)
```

```julia
DevicePRASModel(
    PSY.HydroEnergyReservoir,
    HydroEnergyReservoirFormulation;
    max_active_power="max_active_power",
    inflow="inflow",
    storage_capacity="storage_capacity",
)
```

```julia
DevicePRASModel(
    PSY.HybridSystem,
    HybridSystemFormulation;
    time_series_names=Dict(:max_active_power="max_active_power"),
)
```
"""
function DevicePRASModel(
    ::Type{D},
    ::Type{B};
    time_series_names::Dict{Symbol, String}=Dict{Symbol, String}(),
    kwargs...,
) where {D <: PSY.Device, B <: AbstractPRASFormulation}
    formulation = B(; time_series_names..., kwargs...)
    return DevicePRASModel(D, formulation)
end

"""
Check whether a DevicePRASModel applies to a given type
"""
function appliestodevice(::DevicePRASModel{D}, ::Type{T}) where {D, T}
    return T <: D
end

"""
Uses a PRAS device model to find all components matching it in a system.
"""
function get_available_components(::DevicePRASModel{D}, sys::PSY.System) where {D}
    return PSY.get_available_components(D, sys)
end

"""
    $(TYPEDSIGNATURES)

# Arguments
$(TYPEDFIELDS)

The PRASProblemTemplate contains all configuration necessary for building
a PRAS Simulation from a PowerSystems.jl System.

Since PRAS is an area-based model, we provide a level of aggregation to apply.

PRAS models are processed in reverse order, with later models taking precedence.

# Example

```julia
template = PRASProblemTemplate(
    PSY.Area,
    [
        DevicePRASModel(
            PSY.Generator,
            GeneratorFormulation(max_active_power="max_active_power"),
        ),
        DevicePRASModel(
            PSY.HydroEnergyReservoir,
            HydroEnergyReservoirFormulation(
                max_active_power="max_active_power",
                inflow="inflow",
                storage_capacity="storage_capacity",
            ),
        ),
    ],
)
```
"""
mutable struct PRASProblemTemplate{T <: PSY.AggregationTopology}
    "Level of aggregation to use for PRAS regions"
    aggregation::Type{T}
    "DevicePRASModels to translate components into PRAS"
    device_models::Array{DevicePRASModel}

    function PRASProblemTemplate(
        aggregation::Type{T}=PSY.Area,
        device_models::Array{DevicePRASModel}=DevicePRASModel[],
    ) where {T <: PSY.AggregationTopology}
        return new{T}(aggregation, device_models)
    end
end

"""
    $(TYPEDSIGNATURES)

# Arguments
- `template::PRASProblemTemplate`: Template to add device model to
- `device_model::DevicePRASModel{D}`: Device model to add

Add a device model to a PRASProblemTemplate. If an existing model
already applies to the given device type, then a warning
is issued. However, newer models will take precedence.
"""
function set_device_model!(
    template::PRASProblemTemplate,
    device_model::DevicePRASModel{D},
) where {D}
    for existing_model in template.device_models
        if appliestodevice(existing_model, D)
            @warn "Device model $(D) already exists in template"
        end
        return
    end
    push!(template.device_models, device_model)
end

"""
    $(TYPEDSIGNATURES)

# Arguments
- `template::PRASProblemTemplate`: Template to add device model to
- `::Type{D}`: Device type
- `::Type{B}`: Formulation type

Adds a device model to a PRASProblemTemplate by passing the type
to a constructor.
"""
function set_device_model!(
    template::PRASProblemTemplate,
    ::Type{D},
    ::Type{B},
) where {D <: DevicePRASModel, B <: AbstractPRASFormulation}
    set_device_model!(template, DevicePRASModel(D, B()))
end

"""
    $(SIGNATURES)

Constructs a dictionary from Sienna Devices to formulation objects
"""
function build_component_to_formulation(
    ::Type{B},
    sys::PSY.System,
    device_models::Array{DevicePRASModel},
)::Dict{PSY.Device, B} where {B <: AbstractPRASFormulation}
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
