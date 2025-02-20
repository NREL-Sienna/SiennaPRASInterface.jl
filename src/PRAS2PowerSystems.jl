"""
    $(TYPEDSIGNATURES)

Analyze resource adequacy using Monte Carlo simulation and add the asset status from the worst sample
to PSY.TimeSeriesForcedOutage of the component.

# Arguments

  - `sys::PSY.System`: PowerSystems.jl system model
  - `template::RATemplate`: PRAS problem template
  - `method::PRASCore.SequentialMonteCarlo`: Simulation method to use

# Returns

 - PSY System with PSY.TimeSeriesForcedOutage for all components for which asset status is available
"""
function generate_outage_profile!(
    sys::PSY.System,
    template::RATemplate,
    method::PRASCore.SequentialMonteCarlo,
)
    pras_system = generate_pras_system(sys, template)
    resultsspecs = get_outage_pras_resultspec(template)
    results = PRASCore.assess(pras_system, method, resultsspecs...)
    add_asset_status!(sys, results, template, resultsspecs)
    return sys
end

"""
    generate_outage_profile!(
        sys::PSY.System,
        aggregation::Type{AT},
        method::PRASCore.SequentialMonteCarlo,
    ) where {AT <: PSY.AggregationTopology, RM <: PRASCore.Results.ReliabilityMetric}

Analyze resource adequacy using Monte Carlo simulation and add the asset status from the worst sample
to PSY.TimeSeriesForcedOutage of the component.

# Arguments

  - `sys::PSY.System`: PowerSystems.jl system model
  - `aggregation::Type{AT}`: Aggregation topology to use in translating to PRAS
  - `method::PRASCore.SequentialMonteCarlo`: Simulation method to use

# Returns

  - PSY System with PSY.TimeSeriesForcedOutage for all components for which asset status is available
"""
function generate_outage_profile!(
    sys::PSY.System,
    aggregation::Type{AT},
    method::PRASCore.SequentialMonteCarlo,
) where {AT <: PSY.AggregationTopology}
    template = RATemplate(aggregation, DEFAULT_DEVICE_MODELS)
    sys = generate_outage_profile!(sys, template, method)
    return sys
end

"""
    $(TYPEDSIGNATURES)

Analyze resource adequacy using Monte Carlo simulation and add the asset status from the worst sample
to PSY.TimeSeriesForcedOutage of the component.

Uses default template with PSY.Area AggregationTopology.

# Arguments

  - `sys::PSY.System`: PowerSystems.jl system model
  - `method::PRASCore.SequentialMonteCarlo`: Simulation method to use

# Returns

    - PSY System with PSY.TimeSeriesForcedOutage for all components for which asset status is available
"""
function generate_outage_profile!(sys::PSY.System, method::PRASCore.SequentialMonteCarlo)
    sys = generate_outage_profile!(sys, DEFAULT_TEMPLATE, method)
    return sys
end

"""
    $(TYPEDSIGNATURES)

Add the asset status from the worst sample to PSY.TimeSeriesForcedOutage of the component.

# Arguments

  - `sys::PSY.System`: PowerSystems.jl system model
  - `results::T`: Tuple of results from the PRAS assess call
  - `template::RATemplate`: PRAS problem template
  - `resultsspecs::Vector{PRASCore.Results.ResultSpec}`: PRAS ResultSpecs
  
# Returns

    - PSY System with PSY.TimeSeriesForcedOutage for all components for which asset status is available
"""
function add_asset_status!(
    sys::PSY.System,
    results::T,
    template::RATemplate,
    resultsspecs::Vector{PRASCore.Results.ResultSpec},
) where {T <: Tuple{Vararg{PRASCore.Results.Result}}}
    # Time series timestamps
    all_ts = PSY.get_time_series_multiple(sys, x -> (typeof(x) <: PSY.StaticTimeSeries))
    ts_timestamps = TimeSeries.timestamp(first(all_ts).data)

    shortfall_samp_idx =
        findfirst(type_name.(typeof.(results)) .== PRASCore.Results.ShortfallSamplesResult)
    shortfall_samples = results[shortfall_samp_idx]

    sample_idx = argmax(shortfall_samples[])

    for result_spec in filter(x -> !(x == ShortfallSamples()), resultsspecs)
        device_ramodel = get_device_ramodel(typeof(result_spec))
        gens_to_formula =
            build_component_to_formulation(device_ramodel, sys, template.device_models)

        gens_avail_idx = findfirst(
            type_name.(typeof.(results)) .== get_pras_resulttype(typeof(result_spec)),
        )
        gens_avail = results[gens_avail_idx]

        for gen in keys(gens_to_formula)
            ts_forced_outage =
                PSY.TimeSeriesForcedOutage(; outage_status_scenario="WorstShortfallSample")
            PSY.add_supplemental_attribute!(sys, gen, ts_forced_outage)

            availability_data = TimeSeries.TimeArray(
                ts_timestamps,
                getindex.(gens_avail[gen.name, :], sample_idx),
            )
            availability_timeseries =
                PSY.SingleTimeSeries("availability", availability_data)

            PSY.add_time_series!(sys, ts_forced_outage, availability_timeseries)
            @info "Added availability time series to TimeSeriesForcedOutage supplemental attribute of $(gen.name)."
        end
    end
end

"""
    get_outage_pras_resultspec(
        template::RATemplate,
    )

Get ResultSpec to be used by PRAS.assess() to get asset availability based on the 
template.

# Arguments

  - `template::RATemplate`: PRAS problem template

# Returns

  - PRAS ResultSpec
"""

function get_outage_pras_resultspec(template::RATemplate)
    # We have to change this once we define a formulation for Lines
    resultsspecs = PRASCore.Results.ResultSpec[]
    push!(resultsspecs, ShortfallSamples())

    for model in template.device_models
        result_spec = get_pras_resultspec(typeof(get_formulation(model)))
        if !(result_spec in resultsspecs)
            push!(resultsspecs, result_spec)
        end
    end

    return resultsspecs
end

"""
Get PRAS ResultSpec for RAFormulation
"""
function get_pras_resultspec(::Type{GeneratorPRAS})
    return GeneratorAvailability()
end

"""
Get PRAS ResultSpec for RAFormulation
"""
function get_pras_resultspec(::Type{B}) where {B <: StoragePRAS}
    return StorageAvailability()
end

"""
Get PRAS ResultSpec for RAFormulation
"""
function get_pras_resultspec(::Type{B}) where {B <: GeneratorStoragePRAS}
    return GeneratorStorageAvailability()
end

"""
Get PRAS ResultSpec for AbstractRAFormulation
"""
function get_pras_resultspec(::Type{B}) where {B <: AbstractRAFormulation}
    return error("PRAS ResultSpec not defined for $(B)")
end

"""
Get DeviceRAModel for PRAS ResultSpec
"""
function get_device_ramodel(::Type{GeneratorAvailability})
    return GeneratorPRAS
end

"""
Get DeviceRAModel for PRAS ResultSpec
"""
function get_device_ramodel(::Type{StorageAvailability})
    return StoragePRAS
end

"""
Get DeviceRAModel for PRAS ResultSpec
"""
function get_device_ramodel(::Type{GeneratorStorageAvailability})
    return GeneratorStoragePRAS
end

"""
Get DeviceRAModel for PRAS ResultSpec
"""
function get_device_ramodel(::Type{R}) where {R <: PRASCore.Results.ResultSpec}
    return error("DeviceRAModel not defined for $(R)")
end

"""
Get AvailabilityResult Type for PRAS ResultSpec
"""
function get_pras_resulttype(::Type{GeneratorAvailability})
    return PRASCore.Results.GeneratorAvailabilityResult
end

"""
Get AvailabilityResult Type for PRAS ResultSpec
"""
function get_pras_resulttype(::Type{StorageAvailability})
    return PRASCore.Results.StorageAvailabilityResult
end

"""
Get AvailabilityResult Type for PRAS ResultSpec
"""
function get_pras_resulttype(::Type{GeneratorStorageAvailability})
    return PRASCore.Results.GeneratorStorageAvailabilityResult
end

"""
Get AvailabilityResult Type for PRAS ResultSpec
"""
function get_pras_resulttype(::Type{R}) where {R <: PRASCore.Results.ResultSpec}
    return error("PRAS AbstractAvailabilityResult not defined for $(R)")
end

"""
Get name of Type to strip out parameters for easy comparison
"""
type_name(::Type{T}) where {T} = (isempty(T.parameters) ? T : T.name.wrapper)
