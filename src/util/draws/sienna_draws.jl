# TODO Try Threefry and others (https://juliarandom.github.io/RandomNumbers.jl/stable/man/random123/)
rng = Random123.Philox4x((0, 0), 10)

"""
    $(TYPEDSIGNATURES)

Draw stochastic generator availability profiles and attach them as time series.

For each [`PowerSystems.Generator`](@extref) with a positive
[`PowerSystems.GeometricDistributionForcedOutage`](@extref) transition probability,
simulates sequential outage and recovery states and writes an `availability` time series
to the supplemental attribute.

# Arguments

  - `sys`: [`PowerSystems.System`](@extref) to update in place
  - `initial_time`: start time for generated time series (defaults to system time series metadata when omitted)
  - `resolution`: time step for generated time series
  - `steps`: number of time steps per horizon segment
  - `horizon`: number of horizon segments to simulate

# Returns

  - `nothing`; availability time series are added to `sys`
"""
function make_generator_outage_draws!(
    sys::PSY.System,
    initial_time::Dates.DateTime=nothing,
    resolution::TIMEPERIOD=nothing,
    steps::Int=nothing,
    horizon::Int=nothing,
) where {TIMEPERIOD <: Dates.TimePeriod}
    if ~(resolution in PSY.get_time_series_resolutions(sys))
        error("Cannot generate outage time series with $(resolution) resolution ...")
    end

    outage_gens = PSY.get_components(
        x -> (
            PSY.has_supplemental_attributes(x, PSY.GeometricDistributionForcedOutage) &&
            PSY.get_outage_transition_probability(
                first(
                    PSY.get_supplemental_attributes(
                        PSY.GeometricDistributionForcedOutage,
                        x,
                    ),
                ),
            ) > 0.0
        ),
        PSY.Generator,
        sys,
    )
    ngens = length(collect(outage_gens))

    gens_available = Vector{Bool}(undef, ngens)
    gens_nexttransition = Vector{Int}(undef, ngens)

    initial_availability, next_transition = initialize_availability!(
        rng,
        gens_available,
        gens_nexttransition,
        collect(outage_gens),
        (steps * horizon),
    )
    gen_availability, gen_next_transition = initial_availability, next_transition

    all_gen_availability = Vector{Int64}[]
    for i in 1:(steps * horizon)
        gen_availability, gen_next_transition = update_availability!(
            rng,
            gen_availability,
            gen_next_transition,
            collect(outage_gens),
            i,
            (steps * horizon),
        )
        push!(all_gen_availability, gen_availability)
    end

    for (idx, gen) in enumerate(outage_gens)
        @info "Adding avaialability time series to $(gen.name) ..."
        PSY.add_time_series!(
            sys,
            first(
                PSY.get_supplemental_attributes(PSY.GeometricDistributionForcedOutage, gen),
            ),
            PSY.SingleTimeSeries(
                "availability",
                TimeSeries.TimeArray(
                    get_timestamps(initial_time, resolution, steps, horizon),
                    getindex.(all_gen_availability, idx),
                ),
            ),
        )
    end
end
