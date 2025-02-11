abstract type AbstractPRASFormulation end

struct GeneratorFormulation <: AbstractPRASFormulation
    max_active_power::String
    lump_renewable_generation::Bool

    function GeneratorFormulation(;
        max_active_power="max_active_power",
        lump_renewable_generation=false,
    )
        return new(max_active_power, lump_renewable_generation)
    end
end

abstract type GeneratorStorageFormulation <: AbstractPRASFormulation end

struct HybridSystemFormulation <: GeneratorStorageFormulation
    max_active_power::String

    function HybridSystemFormulation(; max_active_power="max_active_power")
        return new(max_active_power)
    end
end

struct HydroEnergyReservoirFormulation <: GeneratorStorageFormulation
    max_active_power::String
    inflow::String
    storage_capacity::String

    function HydroEnergyReservoirFormulation(;
        max_active_power="max_active_power",
        inflow="inflow",
        storage_capacity="storage_capacity",
    )
        return new(max_active_power, inflow, storage_capacity)
    end
end

struct StorageFormulation <: AbstractPRASFormulation end

struct DevicePRASModel{D <: PSY.Device, B <: AbstractPRASFormulation}
    formulation::B

    function DevicePRASModel(
        ::Type{D},
        formulation::B,
    ) where {D <: PSY.Device, B <: AbstractPRASFormulation}
        return new{D, B}(formulation)
    end
end

function DevicePRASModel(
    ::Type{D},
    ::Type{B};
    time_series_names::Dict{Symbol, String}=Dict{Symbol, String}(),
    kwargs...,
) where {D <: PSY.Device, B <: AbstractPRASFormulation}
    formulation = B(; time_series_names..., kwargs...)
    return DevicePRASModel(D, formulation)
end

function appliestodevice(::DevicePRASModel{D}, ::Type{T}) where {D, T}
    return T <: D
end

function get_available_components(::DevicePRASModel{D}, sys::PSY.System) where {D}
    return PSY.get_available_components(D, sys)
end

mutable struct PRASProblemTemplate{T <: PSY.AggregationTopology}
    aggregation::Type{T}
    device_models::Array{DevicePRASModel}

    function PRASProblemTemplate(
        aggregation::Type{T}=PSY.Area,
        device_models::Array{DevicePRASModel}=DevicePRASModel[],
    ) where {T <: PSY.AggregationTopology}
        return new{T}(aggregation, device_models)
    end
end

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

function set_device_model!(
    template::PRASProblemTemplate,
    ::Type{D},
    ::Type{B},
) where {D <: DevicePRASModel, B <: AbstractPRASFormulation}
    set_device_model!(template, DevicePRASModel(D, B()))
end

function build_component_to_formulation(
    ::Type{B},
    sys::PSY.System,
    device_models::Array{DevicePRASModel},
) where {B <: AbstractPRASFormulation}
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
