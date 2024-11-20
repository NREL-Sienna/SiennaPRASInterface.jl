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
    assess(
        sys::PSY.System,
        aggregation::Type{AT},
        method::PRAS.SimulationSpec,
        resultsspecs::PRAS.ResourceAdequacy.ResultSpec...,
    ) where {AT <: PSY.AggregationTopology}

Estimate resource adequacy using Monte Carlo simulation.

# Arguments

  - `sys::PSY.System`: PowerSystems.jl system model
  - `aggregation::Type{AT}`: Aggregation topology to use in translating to PRAS
  - `method::PRAS.SimulationSpec`: Simulation method to use
  - `resultsspec::PRAS.ResourceAdequacy.ResultSpec...`: Results to compute

# Returns

  - Tuple of results from `resultsspec`: default is ([`ShortfallResult`](@ref),)
"""
function PRAS.assess(
    sys::PSY.System,
    aggregation::Type{AT},
    method::PRAS.ResourceAdequacy.SimulationSpec,
    resultsspecs::PRAS.ResourceAdequacy.ResultSpec...,
) where {AT <: PSY.AggregationTopology}
    pras_system = generate_pras_system(sys, aggregation)
    return PRAS.assess(pras_system, method, resultsspecs...)
end

end
