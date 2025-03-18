"""
    S2P_metadata

Struct to store metadata for the Sienna to PRAS conversion
"""
mutable struct S2P_metadata
    N::Int64
    first_timestamp::Dates.DateTime
    hs_uuids::Vector{Base.UUID}
    pras_resolution::Type{T} where {T <: Dates.Period}
    pras_timestep::Int64
end

function S2P_metadata(df::DataFrames.DataFrame)
    ts_period = df[1, "resolution"].periods[1]
    return S2P_metadata(
        df[1, "time_step_count"],
        Dates.DateTime(df[1, "initial_timestamp"]),
        Vector{Base.UUID}[],
        typeof(ts_period),
        ts_period.value,
    )
end
