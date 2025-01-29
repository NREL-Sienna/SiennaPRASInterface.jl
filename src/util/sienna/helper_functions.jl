"""
    get_available_components_in_aggregation_topology(
        type::Type{<:PSY.StaticInjection},
        sys::PSY.System,
        region::PSY.AggregationTopology,
    )

Get available components in the AggregationTopology region of the given type.
"""
function get_available_components_in_aggregation_topology(
    type::Type{<:PSY.StaticInjection},
    sys::PSY.System,
    region::PSY.AggregationTopology,
)
    avail_comps = filter(
        x -> (PSY.get_available(x)),
        collect(PSY.get_components_in_aggregation_topology(type, sys, region)),
    )
    return avail_comps
end

"""
    get_generator_category(gen::StaticInjection)

Get the category of the generator.

# Arguments

  - `gen::StaticInjection`: Generator

# Returns

  - `String`: Category of the generator
"""
function get_generator_category(gen::PSY.StaticInjection)
    error("get_generator_category isn't defined for $(typeof(gen))")
end

function get_generator_category(gen::GEN) where {GEN <: PSY.RenewableGen}
    return string(PSY.get_prime_mover_type(gen))
end

function get_generator_category(gen::GEN) where {GEN <: PSY.ThermalGen}
    return string(PSY.get_fuel(gen))
end

function get_generator_category(gen::GEN) where {GEN <: PSY.HydroGen}
    return "Hydro"
end

function get_generator_category(stor::GEN) where {GEN <: PSY.Storage}
    if (occursin("Distributed", PSY.get_name(stor)))
        return "Distributed_Storage"
    elseif (occursin("Battery", PSY.get_name(stor)))
        return "Battery_Storage"
    else
        return "Battery"
    end
end

function get_generator_category(stor::GEN) where {GEN <: PSY.HybridSystem}
    return "Hybrid-System"
end

"""
    line_rating(line::Branch)

Get the line rating.

# Arguments

  - `line::Branch`: Line

# Returns

  - `Tuple{forward_capacity::Float64, backward_capacity::Float64}`: Line rating
"""
function line_rating(line::PSY.Branch)
    error("line_rating isn't defined for $(typeof(line))")
end

function line_rating(line::PSY.Line)
    rate = PSY.get_rating(line)
    return (forward_capacity=abs(rate), backward_capacity=abs(rate))
end

function line_rating(line::PSY.MonitoredLine)
    rate = PSY.get_flow_limits(line)
    return (
        forward_capacity=abs(getfield(rate, :from_to)),
        backward_capacity=abs(getfield(rate, :to_from)),
    )
end

function line_rating(line::PSY.TwoTerminalHVDCLine)
    forward_capacity = getfield(PSY.get_active_power_limits_from(line), :max)
    backward_capacity = getfield(PSY.get_active_power_limits_to(line), :max)
    return (
        forward_capacity=abs(forward_capacity),
        backward_capacity=abs(backward_capacity),
    )
end

function line_rating(line::DCLine) where {DCLine <: HVDCLineTypes}
    error("line_rating isn't defined for $(typeof(line))")
end

"""
    line_type(line::Branch)

Get the line type.

# Arguments

  - `line::Branch`: Line

# Returns

  - `String`: Line type
"""
function line_type(line::PSY.Branch)
    error("line_type isn't defined for $(typeof(line))")
end

function line_type(line::Union{PSY.Line, PSY.MonitoredLine})
    return "ACLine"
end

function line_type(line::PSY.TwoTerminalHVDCLine)
    return "HVDCLine"
end

function line_type(line::DCLine) where {DCLine <: HVDCLineTypes}
    error("line_type isn't defined for $(typeof(line))")
end

"""
    get_ts_values(ts::PSY.TimeSeriesData)

Get the time series values.

# Arguments

  - `ts::TimeSeries`: Time series

# Returns

  - `Array{Float64}`: Time series values
"""
function get_ts_values(ts::PSY.TimeSeriesData)
    error("get_ts_values isn't defined for $(typeof(ts))")
end

function get_ts_values(ts::TS) where {TS <: PSY.AbstractDeterministic}
    if (typeof(ts) == PSY.DeterministicSingleTimeSeries)
        forecast_vals = get_ts_values(ts.single_time_series)
    else
        forecast_vals = []
        for it in collect(keys(PSY.get_data(ts)))
            append!(
                forecast_vals,
                collect(values(PSY.get_window(ts, it; len=interval_len))),
            )
        end
    end
    return forecast_vals
end

function get_ts_values(ts::TS) where {TS <: PSY.StaticTimeSeries}
    forecast_vals = values(PSY.get_data(ts))
    return forecast_vals
end

"""
    get_first_ts(ts::TS)

Get time series and handle components with no time series (such as unavailable)
"""
function get_first_ts(ts::TS) where {TS <: Channel{Any}}
    if isempty(ts)
        return nothing
    else
        return first(ts)
    end
end
function get_forecast_values(ts::Nothing)
    return zeros(length(period_of_interest))
end

function get_outage_time_series_data(
    gen::Union{PSY.StaticInjection, PSY.Branch},
    s2p_meta::S2P_metadata,
)
    # Get GeometricForcedOutage SupplementalAttribute of the generator g
    outage_sup_attrs =
        PSY.get_supplemental_attributes(PSY.GeometricDistributionForcedOutage, gen)
    λ_gen, μ_gen = zeros(Float64, 1, s2p_meta.N), ones(Float64, 1, s2p_meta.N)
    if (length(outage_sup_attrs) > 0)
        transition_data = first(outage_sup_attrs)
        λ = PSY.get_outage_transition_probability(transition_data)
        μ = if (iszero(PSY.get_mean_time_to_recovery(transition_data)))
            1.0
        else
            1 / PSY.get_mean_time_to_recovery(transition_data)
        end

        λ_gen, μ_gen = if (PSY.has_time_series(transition_data, PSY.SingleTimeSeries))
            PSY.get_time_series_values(
                PSY.SingleTimeSeries,
                transition_data,
                "outage_probability",
            ),
            PSY.get_time_series_values(
                PSY.SingleTimeSeries,
                transition_data,
                "recovery_probability",
            )
        else
            fill.(λ, 1, s2p_meta.N), fill.(μ, 1, s2p_meta.N)
        end
    else
        @warn "No GeometricForcedOutage SupplementalAttribute available for $(PSY.get_name(gen)) of $(typeof(gen)). Using nominal outage and recovery probabilities for this component."
    end

    return λ_gen, μ_gen
end
