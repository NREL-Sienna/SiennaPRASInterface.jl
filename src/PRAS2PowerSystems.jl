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
    add_asset_status!(sys, SPIOutageResult(results), template, resultsspecs)
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
    results::SPIOutageResult,
    template::RATemplate,
    resultsspecs::Vector{PRASCore.Results.ResultSpec},
)
    # Time series timestamps
    all_ts = PSY.get_time_series_multiple(sys, x -> (typeof(x) <: PSY.StaticTimeSeries))
    ts_timestamps = TimeSeries.timestamp(first(all_ts).data)

    shortfall_samples = results.shortfall_samples

    sample_idx = argmax(shortfall_samples[])

    for result_spec in filter(x -> !(x == ShortfallSamples()), resultsspecs)
        device_ramodel = get_device_ramodel(typeof(result_spec))
        gens_to_formula =
            build_component_to_formulation(device_ramodel, sys, template.device_models)

        gens_avail = getproperty(results, get_outage_result_fname(typeof(result_spec)))

        if !(isnothing(gens_avail))
            for gen in keys(gens_to_formula)
                ts_forced_outage = PSY.TimeSeriesForcedOutage(;
                    outage_status_scenario="WorstShortfallSample",
                )
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
end
