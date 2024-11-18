module PRASInterface
#################################################################################
# Exports
#################################################################################
export generate_pras_system
export generate_outage_profile
export generate_csv_outage_profile
export add_csv_time_series!
export add_csv_time_series_single_stage!
export make_generator_outage_draws!
export PRAS
export estimate_resource_adequacy
#################################################################################
# Imports
#################################################################################
using PowerSystems: PowerSystems
using Dates: Dates
using TimeZones: TimeZones
using DataFrames: DataFrames
using CSV: CSV
using JSON: JSON
using UUIDs: UUIDs
using TimeSeries: TimeSeries
using Random123: Random123
using Random: Random

const PSY = PowerSystems
#################################################################################
# Includes
#################################################################################

module PRAS
using Reexport
const PRAS_VERSION = "v0.6.0"
include("MinCostFlows/MinCostFlows.jl")
include("PRASBase/PRASBase.jl")
include("ResourceAdequacy/ResourceAdequacy.jl")
include("CapacityCredit/CapacityCredit.jl")
end

include("util/definitions.jl")
include("util/runchecks.jl")

include("util/parsing/Sienna_PRAS_metadata.jl")
include("util/parsing/lines_and_interfaces.jl")
include("util/parsing/outage_data_helper_functions.jl")
include("util/parsing/PRAS_export.jl")

include("util/sienna/helper_functions.jl")
include("util/sienna/add_csv_time_series_data.jl")

include("util/draws/draw_helper_functions.jl")
include("util/draws/sienna_draws.jl")

include("PSY2PRAS.jl")
include("PRAS2PSY.jl")

"""
    estimate_resource_adequacy(
        sys::PSY.System;
        samples::Int=10_000,
        resultsspecs::Vector{PRAS.ResultSpec}=[Shortfall],
        seed::Int=0,
        verbose::Bool=false,
        threaded::Bool=true,
        aggregation::Type{AT}=PSY.Area,
        availability::Bool=true,
        lump_region_renewable_gens::Bool=false,
        export_location::Union{Nothing, String}=nothing,
    ) where {AT <: PSY.AggregationTopology}

Estimate resource adequacy using Monte Carlo simulation.

# Arguments

  - `sys::PSY.System`: PowerSystems.jl system model
  - `samples::Int`: Number of samples to use in the Monte Carlo simulation
  - `resultsspecs::Vector{PRAS.ResourceAdequacy.ResultSpec}`: Vector of result specifications to use in PRAS simulation
  - `seed::Int`: Random seed for the Monte Carlo simulation
  - `verbose::Bool`: Print progress information in PRAS simulation
  - `threaded::Bool`: Use threading in PRAS simulation
  - `aggregation::Type{AT}`: Aggregation topology to use in translating to PRAS
  - `availability::Bool`: Use available components
  - `lump_region_renewable_gens::Bool`: Lump renewable generators in PRAS model
  - `export_location::Union{Nothing, String}`: Location to export PRAS system. Default is no export.

# Returns

  - Tuple of results from `resultsspec`: default is ([`ShortfallResult`](@ref),)
"""
function estimate_resource_adequacy(
    sys::PSY.System;
    samples::Int=10_000,
    resultsspecs::Vector{PRAS.ResourceAdequacy.ResultSpec}=PRAS.ResourceAdequacy.ResultSpec[PRAS.Shortfall()],
    seed::Int=0,
    verbose::Bool=false,
    threaded::Bool=true,
    aggregation::Type{AT}=PSY.Area,
    availability::Bool=true,
    lump_region_renewable_gens::Bool=false,
    export_location::Union{Nothing, String}=nothing,
) where {AT <: PSY.AggregationTopology}
    pras_system = generate_pras_system(
        sys,
        aggregation;
        availability,
        lump_region_renewable_gens,
        export_location,
    )
    method = PRAS.SequentialMonteCarlo(; samples, seed, verbose, threaded)
    return PRAS.assess(pras_system, method, resultsspecs...)
end

end
