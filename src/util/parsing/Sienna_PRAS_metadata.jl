"""
    S2P_metadata

Struct to store metadata for the Sienna to PRAS conversion
"""
mutable struct S2P_metadata
    has_static_timeseries::Bool
    has_forecasts::Bool
    filter_func::Function
    N::Int64
    first_timestamp::Union{Nothing, Dates.DateTime}
    first_timeseries::Union{Nothing, Union{<:PSY.Forecast, <:PSY.StaticTimeSeries}}
    hs_uuids::Vector{Base.UUID}
    pras_resolution::Type{T} where {T <: Dates.Period}
    pras_timestep::Int64

    S2P_metadata(
        has_st_timeseries=false,
        has_forecasts=false,
        filter_func=x -> (typeof(x) <: PSY.StaticTimeSeries),
        N=0,
        first_timestamp=nothing,
        first_ts=nothing,
        hs_uuids=Vector{Base.UUID}[],
        pras_res=Dates.Hour,
        pras_ts=1,
    ) = new(
        has_st_timeseries,
        has_forecasts,
        filter_func,
        N,
        first_timestamp,
        first_ts,
        hs_uuids,
        pras_res,
        pras_ts,
    )
end

function add_N!(s2p_meta::S2P_metadata)
    if (s2p_meta.has_static_timeseries)
        s2p_meta.N = length(s2p_meta.first_timeseries.data)
    end
    if (s2p_meta.has_forecasts)
        N = if ~(isnothing(s2p_meta.first_timeseries.single_time_series))
            length(s2p_meta.first_timeseries.single_time_series.data)
        else
            s2p_meta.first_timeseries.count.interval_len
        end
        s2p_meta.N = N
    end
end
