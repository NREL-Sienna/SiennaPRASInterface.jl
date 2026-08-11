"""
PowerSystems Interface for Probabilistic Resource Adequacy Studies (PRAS)

# Key Functions

  - [`generate_pras_system`](@ref): convert PSY to PRAS model
  - [`assess`](@ref): assess PRAS model

# Key PRAS Types

  - [`PRASCore.Systems.SystemModel`](@extref): PRAS data structure
  - [`PRASCore.Simulations.SequentialMonteCarlo`](@extref): method for PRAS analysis
  - [`PRASCore.Results.Shortfall`](@extref): PRAS metric for missing generation
  - [`PRASCore.Results.LOLE`](@extref): PRAS metric for loss of load expectation
  - [`PRASCore.Results.EUE`](@extref): PRAS metric for energy unserved expectation
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
export generate_outage_profile!
export make_generator_outage_draws!

export GeneratorPRAS
export HybridSystemPRAS
export HydroEnergyReservoirPRAS
export DeviceRAModel
export GeneratorStoragePRAS
export EnergyReservoirSoC
export StoragePRAS
export LinePRAS
export AreaInterchangeLimit
export StaticLoadPRAS
export set_device_model!
export RATemplate

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
using DocStringExtensions

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

include("util/runchecks.jl")

include("util/parsing/Sienna_PRAS_metadata.jl")
include("util/parsing/lines_and_interfaces.jl")
include("util/parsing/outage_data_helper_functions.jl")
include("util/parsing/PRAS_export.jl")

include("formulation_definitions.jl")

include("util/sienna/helper_functions.jl")

include("util/draws/draw_helper_functions.jl")
include("util/draws/sienna_draws.jl")

include("util/definitions.jl")
include("PowerSystems2PRAS.jl")

include("util/parsing/result_export_helper_functions.jl")
include("PRAS2PowerSystems.jl")

"""
    $(TYPEDSIGNATURES)

Analyze resource adequacy using Monte Carlo simulation.

# Arguments

  - `sys`: [`PowerSystems.System`](@extref) to translate and assess
  - `aggregation`: [`PowerSystems.AggregationTopology`](@extref) type used for PRAS region aggregation
  - `method`: [`PRASCore.Simulations.SequentialMonteCarlo`](@extref) simulation method
  - `resultsspecs`: [PRAS result specifications](@extref PRASCore :doc:`PRAS/results`) to compute (for example [`PRASCore.Results.Shortfall`](@extref))

# Returns

  - Tuple of result objects, one per requested result specification (for example [`PRASCore.Results.Shortfall`](@extref) when that is the only specification)
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

"""
    $(TYPEDSIGNATURES)

Analyze resource adequacy using Monte Carlo simulation.

# Arguments

  - `sys`: [`PowerSystems.System`](@extref) to translate and assess
  - `template`: [`RATemplate`](@ref) defining aggregation topology and device mappings
  - `method`: [`PRASCore.Simulations.SequentialMonteCarlo`](@extref) simulation method
  - `resultsspecs`: [PRAS result specifications](@extref PRASCore :doc:`PRAS/results`) to compute (for example [`PRASCore.Results.Shortfall`](@extref))

# Returns

  - Tuple of result objects, one per requested result specification (for example [`PRASCore.Results.Shortfall`](@extref) when that is the only specification)
"""
function PRASCore.assess(
    sys::PSY.System,
    template::RATemplate,
    method::PRASCore.SequentialMonteCarlo,
    resultsspecs::PRASCore.Results.ResultSpec...,
)
    pras_system = generate_pras_system(sys, template)
    return PRASCore.assess(pras_system, method, resultsspecs...)
end

"""
    $(TYPEDSIGNATURES)

Analyze resource adequacy using Monte Carlo simulation.

Uses default template with [`PowerSystems.Area`](@extref) level aggregation.

# Arguments

  - `sys`: [`PowerSystems.System`](@extref) to translate and assess
  - `method`: [`PRASCore.Simulations.SequentialMonteCarlo`](@extref) simulation method
  - `resultsspecs`: [PRAS result specifications](@extref PRASCore :doc:`PRAS/results`) to compute (for example [`PRASCore.Results.Shortfall`](@extref))

# Returns

  - Tuple of result objects, one per requested result specification (for example [`PRASCore.Results.Shortfall`](@extref) when that is the only specification)
"""
function PRASCore.assess(
    sys::PSY.System,
    method::PRASCore.SequentialMonteCarlo,
    resultsspecs::PRASCore.Results.ResultSpec...,
)
    pras_system = generate_pras_system(sys, DEFAULT_TEMPLATE)
    return PRASCore.assess(pras_system, method, resultsspecs...)
end

end
