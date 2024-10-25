#######################################################
# Structs to parse and store the outage information
#######################################################
mutable struct S2P_metadata
    has_static_timeseries::Bool
    has_forecasts::Bool
    N::Int64
    first_timestamp::Union{Nothing, Dates.DateTime}
    first_timeseries::Union{Nothing, Union{<:PSY.Forecast, <:PSY.StaticTimeSeries}}
    
    S2P_metadata(has_st_timeseries = true, has_forecasts = false, N = 0, first_timestamp = nothing, first_ts = nothing) =new(has_st_timeseries,has_forecasts,N,first_timestamp, first_ts)
end 

function add_N!(s2p_meta::S2P_metadata)
    if (s2p_meta.has_static_timeseries)
        s2p_meta.N = length(s2p_meta.first_timeseries.data)
    end
    if (s2p_meta.has_forecasts)
        N = 
        if ~(isnothing(s2p_meta.first_timeseries.single_time_series))
            length(s2p_meta.first_timeseries.single_time_series.data)
        else
            s2p_meta.first_timeseries.count.interval_len
        end
        s2p_meta.N = N
    end
end