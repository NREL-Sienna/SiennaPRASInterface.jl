"""
PowerSystems Interface for Probabilistic Resource Adequacy Studies (PRAS)

# Key Functions

  - [`generate_pras_system`](@ref): convert PSY to PRAS model
  - [`assess`](@ref): assess PRAS model

# Key PRAS Types

  - [`SystemModel`](@ref): PRAS data structure
  - [`SequentialMonteCarlo`](@ref): method for PRAS analysis
  - [`Shortfall`](@ref): PRAS metric for missing generation
  - [`LOLE`](@ref): PRAS metric for loss of load expectation
  - [`EUE`](@ref): PRAS metric for energy unserved expectation
"""
module SiennaPRASInterface
#################################################################################
# Exports
#################################################################################
export generate_pras_system
export SystemModel
export assess
export SequentialMonteCarlo
# ResultSpecs
export Shortfall
export ShortfallSamples
export Surplus
export SurplusSamples
export Flow
export FlowSamples
export Utilization
export UtilizationSamples
export StorageEnergy
export StorageEnergySamples
export GeneratorStorageEnergy
export GeneratorStorageEnergySamples
export GeneratorAvailability
export StorageAvailability
export GeneratorStorageAvailability
export LineAvailability

export LOLE
export EUE
export val
export stderror
export generate_outage_profile
export generate_csv_outage_profile
export add_csv_time_series!
export add_csv_time_series_single_stage!
export make_generator_outage_draws!

#################################################################################
# Imports
#################################################################################
import PowerSystems
import Dates
import TimeZones
import DataFrames
import CSV
import JSON
import UUIDs
import TimeSeries
import Random123
import Random

const PSY = PowerSystems
#################################################################################
# Includes
#################################################################################

import PRASCore

import PRASCore:
    assess,
    LOLE,
    EUE,
    val,
    stderror,
    SequentialMonteCarlo,
    Shortfall,
    ShortfallSamples,
    Surplus,
    SurplusSamples,
    Flow,
    FlowSamples,
    Utilization,
    UtilizationSamples,
    StorageEnergy,
    StorageEnergySamples,
    GeneratorStorageEnergy,
    GeneratorStorageEnergySamples,
    GeneratorAvailability,
    StorageAvailability,
    GeneratorStorageAvailability,
    LineAvailability,
    SystemModel

import PRASFiles

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

include("PowerSystems2PRAS.jl")
include("PRAS2PowerSystems.jl")

"""
    assess(
        sys::PSY.System,
        aggregation::Type{AT},
        method::PRASCore.SequentialMonteCarlo,
        resultsspecs::PRASCore.Results.ResultSpec...,
    ) where {AT <: PSY.AggregationTopology}

Estimate resource adequacy using Monte Carlo simulation.

# Arguments

  - `sys::PSY.System`: PowerSystems.jl system model
  - `aggregation::Type{AT}`: Aggregation topology to use in translating to PRAS
  - `method::PRASCore.SequentialMonteCarlo`: Simulation method to use
  - `resultsspec::PRASCore.Results.ResultSpec...`: Results to compute

# Returns

  - Tuple of results from `resultsspec`: default is ([`ShortfallResult`](@ref),)
"""
function PRASCore.assess(
    sys::PSY.System,
    aggregation::Type{AT},
    method::PRASCore.SequentialMonteCarlo,
    resultsspecs::PRASCore.Results.ResultSpec...,
) where {AT <: PSY.AggregationTopology}
    pras_system = generate_pras_system(sys, aggregation)
    return PRASCore.assess(pras_system, method, resultsspecs...)
end

end
