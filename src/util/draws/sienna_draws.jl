###########################################
# Add Required Packages
###########################################
import PowerSystems
import PowerSystemCaseBuilder
import CSV
import DataFrames
import TimeSeries: TimeArray
import Dates: DateTime, Hour, TimePeriod

const PSY = PowerSystems
const PSCB = PowerSystemCaseBuilder

include("draw_helper_functions.jl")

###########################################
# Build RTS-GMLC System using PSCB
###########################################
rts_da_sys = PSCB.build_system(PSCB.PSISystems, "RTS_GMLC_DA_sys");
PSY.set_units_base_system!(rts_da_sys, "natural_units")

###########################################
# Parse the gen.csv and add OutageData 
# SupplementalAttribute to components for 
# which we have this data 
###########################################
gen_for_data = CSV.read("gen.csv", DataFrames.DataFrame);

for row in DataFrames.eachrow(gen_for_data)
    λ, μ = rate_to_probability(row.FOR, row["MTTR Hr"])
    transition_data = PSY.GeometricDistributionForcedOutage(;
        mean_time_to_recovery=row["MTTR Hr"],
        outage_transition_probability=λ,
    )
    comp = PSY.get_component(PSY.Generator, rts_da_sys, row["GEN UID"])

    if ~(isnothing(comp))
        PSY.add_supplemental_attribute!(rts_da_sys, comp, transition_data)
        @info "Added outage data supplemental attribute to $(row["GEN UID"]) generator"
    else
        @warn "$(row["GEN UID"]) generator doesn't exist in the System."
    end
end

##############################################
# Main function to make generator outage draws
##############################################
import Random123

# TODO Try Threefry and others (https://juliarandom.github.io/RandomNumbers.jl/stable/man/random123/)
rng = Random123.Philox4x((0, 0), 10)

function make_generator_outage_draws!(
    sys::PSY.System,
    initial_time::DateTime=nothing,
    resolution::TIMEPERIOD=nothing,
    steps::Int=nothing,
    horizon::Int=nothing,
) where {TIMEPERIOD <: TimePeriod}
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
                TimeArray(
                    get_timestamps(initial_time, resolution, steps, horizon),
                    getindex.(all_gen_availability, idx),
                ),
            ),
        )
    end
end

##############################################
# Make generator outage draws
##############################################
initial_time = DateTime("2020-01-01") # Initial timestamp similar to DecisionModel etc.
resolution = Hour(1) # Because now you can have mutiple resolutions of time series in a System
horizon = 24 # Horizon similar to DecisionModel etc.
steps = 2 # Number of steps similar to a simulation

make_generator_outage_draws!(rts_da_sys, initial_time, resolution, steps, horizon)
