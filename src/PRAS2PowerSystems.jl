"""
    generate_outage_profile(pras_system,num_runs,psy_sys,num_scenarios,location)

Process the assess results to get timeseries of generator status and include
this timeseries data to the corresponding component in PSY System and exported
using to_json method (serializing the PSY System).

...

# Arguments

  - `pras_system::PRASCore.SystemModel`: PRAS System
  - `num_runs::Int64`: Number of PRAS runs
  - `psy_sys::PSY.System`: PSY System
  - `num_scenarios::Int64`: Number of scenarios of user interest.
  - `location::String`: Location to store outage profile.
    ...

# Examples

```julia-repl
julia> generate_outage_profile(results, pras_sys, psy_sys, 1)
PSY System exported using to_json method in InfrastructureSystems
```
"""
function generate_outage_profile(
    pras_system::PRASCore.SystemModel,
    psy_sys::PSY.System;
    location::Union{Nothing, String}=nothing,
    num_runs::Union{Nothing, Int}=nothing,
    num_scenarios::Union{Nothing, Int}=nothing,
)

    #kwargs handling
    if (location === nothing)
        location = dirname(dirname(@__DIR__))
        @warn "Location to save the exported PSY System not specified. Using the data folder of the module."
    end

    if (num_runs === nothing)
        num_runs = 10000
        @warn "Number of samples of PRAS Assess not specified. Using 10,000 samples."
    end

    if (num_scenarios === nothing)
        num_scenarios = 1
        @warn "Number of scenarios to be exported not specified. Only exporting one scenario."
    end

    # Run PRAS Analysis
    @info "Running PRAS SequentialMonteCarlo Resource Adequacy analysis for $(num_runs) runs ..."
    shortfall_samples, gens_avail, stors_avail, gen_stors_avail, lines_avail =
        PRASCore.assess(
            pras_system,
            PRASCore.SequentialMonteCarlo(samples=num_runs, threaded=true, verbose=false),
            PRASCore.ShortfallSamples(),
            PRASCore.GeneratorAvailability(),
            PRASCore.StorageAvailability(),
            PRASCore.GeneratorStorageAvailability(),
            PRASCore.LineAvailability(),
        )
    @info "Successfully completed PRAS Runs. Exporting outage profiles now..."
    # Setup to save the new PSY System with scenario timeseries data
    #working_dir = pwd();
    working_dir = location
    dt_now = Dates.format(Dates.now(), "dd-u-yy-H-M-S")
    sys_name = string(psy_sys.internal.uuid)

    dir_name =
        working_dir * "/data/Generated-Outage-Profile-JSON/" * sys_name * "/" * dt_now
    mkpath(dir_name)
    #  **TODO: IF generating systems for multiple scenarios, remember to delete the timeseries data from previous scenario.
    asset_dict = Dict([
        (:generators, (gens_avail, PSY.Generator)),
        (:storages, (stors_avail, PSY.Storage)),
        (:generatorstorages, (gen_stors_avail, PSY.StaticInjection)),
        (:lines, (lines_avail, PSY.Branch)),
    ])
    asset_keys = []

    for k in keys(asset_dict)
        if (length(getfield(pras_system, k)) != 0)
            push!(asset_keys, k)
        end
    end
    @info "Scenarios of interest will be sorted according to sample unserved energy. The PSY systems for all the scenarios will be exported using to_json() method in InfrastructureSystems."

    sample_idx = sortperm(shortfall_samples[], rev=true)
    resolution = Dates.Hour(1)
    for i in 1:num_scenarios
        for k in asset_keys
            asset_status_all = asset_dict[k][1].available[:, :, sample_idx[i]]
            asset_names = getfield(asset_dict[k][1], k)
            for (j, asset_name) in enumerate(asset_names)
                # Creating TimeSeries Data
                timestamps = range(
                    Dates.DateTime(pras_system.timestamps.start),
                    step=resolution,
                    length=length(pras_system.timestamps),
                )
                availability_data = TimeSeries.TimeArray(timestamps, asset_status_all[j, :])
                availability_timeseries =
                    PSY.SingleTimeSeries("availability", availability_data)
                #Removing TimeSeries Data associated with the component if it already exists
                if (i != 1)
                    PSY.remove_time_series!(
                        psy_sys,
                        PSY.SingleTimeSeries,
                        PSY.get_component(asset_dict[k][2], psy_sys, asset_name),
                        "availability",
                    )
                end
                # Adding TimeSeries Data to PSY System
                PSY.add_time_series!(
                    psy_sys,
                    PSY.get_component(asset_dict[k][2], psy_sys, asset_name),
                    availability_timeseries,
                )
            end
        end

        file_name = dir_name * "/$i.json"
        PSY.IS.to_json(psy_sys, file_name, force=false)
        @info "Succesfully exported PSY System with outage profile for scenario $(i)."
    end
    @info "Succesfully exported PSY System(s). The .json files and corresponding data .h5 files for all scenarios are available here : $(dir_name)."
end

"""
    generate_outage_profile(pras_system,num_runs,psy_sys,num_scenarios,location)

Process the assess results to get timeseries of generator status and include
this timeseries data to the corresponding component in PSY System and exported
using to_json method (serializing the PSY System).

...

# Arguments

  - `pras_system::PRASCore.SystemModel`: PRAS System
  - `num_runs::Int64`: Number of PRAS runs
  - `psy_sys::PSY.System`: PSY System
  - `num_scenarios::Int64`: Number of scenarios of user interest.
  - `location::String`: Location to store outage profile.
    ...

# Examples

```julia-repl
julia> generate_outage_profile(results, pras_sys, psy_sys, 1)
PSY System exported using to_json method in InfrastructureSystems
```
"""
function generate_csv_outage_profile(
    pras_system::PRASCore.SystemModel;
    location::Union{Nothing, String}=nothing,
    num_runs::Union{Nothing, Int}=nothing,
    num_scenarios::Union{Nothing, Int}=nothing,
)

    #kwargs handling
    if (location === nothing)
        location = dirname(dirname(@__DIR__))
        @warn "Location to save the exported PSY System not specified. Using the data folder of the module."
    end

    if (num_runs === nothing)
        num_runs = 10000
        @warn "Number of samples of PRAS Assess not specified. Using 10,000 samples."
    end

    if (num_scenarios === nothing)
        num_scenarios = 1
        @warn "Number of scenarios to be exported not specified. Only exporting one scenario."
    end

    # Run PRAS Analysis
    @info "Running PRAS SequentialMonteCarlo Resource Adequacy analysis for $(num_runs) runs..."
    shortfall_samples, gens_avail, stors_avail, gen_stors_avail, lines_avail =
        PRASCore.assess(
            pras_system,
            PRASCore.SequentialMonteCarlo(samples=num_runs, threaded=true, verbose=false),
            PRASCore.ShortfallSamples(),
            PRASCore.GeneratorAvailability(),
            PRASCore.StorageAvailability(),
            PRASCore.GeneratorStorageAvailability(),
            PRASCore.LineAvailability(),
        )
    @info "Successfully completed PRAS Runs. Exporting outage profiles now..."
    # Setup to save the new PSY System with scenario timeseries data
    #working_dir = pwd();
    working_dir = location
    dt_now = Dates.format(Dates.now(), "dd-u-yy-H-M-S")

    dir_name = joinpath(
        working_dir,
        "data",
        "Generated-Outage-Profile-JSON",
        string(UUIDs.uuid4()),
        dt_now,
    )

    #  **TODO: IF generating systems for multiple scenarios, remember to delete the timeseries data from previous scenario.
    asset_dict = Dict([
        (:generators, (gens_avail, PSY.Generator)),
        (:storages, (stors_avail, PSY.Storage)),
        (:generatorstorages, (gen_stors_avail, PSY.StaticInjection)),
        (:lines, (lines_avail, PSY.Branch)),
    ])
    asset_keys = []

    for k in keys(asset_dict)
        if (length(getfield(pras_system, k)) != 0)
            push!(asset_keys, k)
        end
    end
    @info "Scenarios of interest will be sorted according to sample unserved energy. The availability for individual asset types will be exported to CSV sheets."

    sample_idx = sortperm(shortfall_samples[], rev=true)

    for i in 1:num_scenarios
        mkpath(joinpath(dir_name, string(i)))

        for k in asset_keys
            asset_status_all = asset_dict[k][1].available[:, :, sample_idx[i]]
            asset_names = getfield(asset_dict[k][1], k)

            df_outage = DataFrames.DataFrame()

            for (j, asset_name) in enumerate(asset_names)
                df_outage[!, asset_name] = asset_status_all[j, :]
            end

            csv_path = joinpath(dir_name, string(i), string(asset_dict[k][2], ".csv"))
            CSV.write(csv_path, df_outage, writeheader=true)
        end
        @info "Succesfully exported the outage profile for scenario $(i)."
    end
    @info "Succesfully exported the outage profiles. The CSV files for all asset types for all scenarios are available here : $(dir_name)."
end

"""
    generate_outage_profile(
        sys::PSY.System,
        aggregation::Type{AT},
        method::PRASCore.SequentialMonteCarlo,
        metric::Type{RM}
    ) where {AT <: PSY.AggregationTopology, RM <: PRASCore.Results.ReliabilityMetric}

Analyze resource adequacy using Monte Carlo simulation and add the asset status from the worst sample
to PSY.TimeSeriesForcedOutage of the component.

# Arguments

  - `sys::PSY.System`: PowerSystems.jl system model
  - `aggregation::Type{AT}`: Aggregation topology to use in translating to PRAS
  - `method::PRASCore.SequentialMonteCarlo`: Simulation method to use
  - `metric::Type{RM}` : ReliabilityMetric to use for sorting

# Returns

  - PSY System with PSY.TimeSeriesForcedOutage for all components for which asset status is available
"""
function generate_outage_profile(
    sys::PSY.System,
    aggregation::Type{AT},
    method::PRASCore.SequentialMonteCarlo,
    metric::Type{RM}
) where {AT <: PSY.AggregationTopology, RM <: PRASCore.Results.ReliabilityMetric}
    pras_system = generate_pras_system(sys, aggregation)
    return PRASCore.assess(pras_system, method, resultsspecs...)
end

"""
    $(TYPEDSIGNATURES)

Analyze resource adequacy using Monte Carlo simulation and add the asset status from the worst sample
to PSY.TimeSeriesForcedOutage of the component.

# Arguments

  - `sys::PSY.System`: PowerSystems.jl system model
  - `template::RATemplate`: PRAS problem template
  - `method::PRASCore.SequentialMonteCarlo`: Simulation method to use
  - `metric::Type{RM}` : ReliabilityMetric to use for sorting

# Returns

 - PSY System with PSY.TimeSeriesForcedOutage for all components for which asset status is available
"""
function generate_outage_profile(
    sys::PSY.System,
    template::RATemplate,
    method::PRASCore.SequentialMonteCarlo,
    metric::Type{RM}
) where {RM <: PRASCore.Results.ReliabilityMetric}
    pras_system = generate_pras_system(sys, template)
    return PRASCore.assess(pras_system, method, resultsspecs...)
end

"""
    $(TYPEDSIGNATURES)

Analyze resource adequacy using Monte Carlo simulation and add the asset status from the worst sample
to PSY.TimeSeriesForcedOutage of the component.

Uses default template with PSY.Area AggregationTopology.

# Arguments

  - `sys::PSY.System`: PowerSystems.jl system model
  - `method::PRASCore.SequentialMonteCarlo`: Simulation method to use
  - `metric::Type{RM}` : ReliabilityMetric to use for sorting

# Returns

    - PSY System with PSY.TimeSeriesForcedOutage for all components for which asset status is available
"""
function generate_outage_profile(
    sys::PSY.System,
    method::PRASCore.SequentialMonteCarlo,
    metric::Type{RM}
) where {RM <: PRASCore.Results.ReliabilityMetric}
    pras_system = generate_pras_system(sys, DEFAULT_TEMPLATE)
    return PRASCore.assess(pras_system, method, resultsspecs...)
end