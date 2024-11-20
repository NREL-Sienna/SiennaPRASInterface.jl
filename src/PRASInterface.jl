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
module PRASInterface
#################################################################################
# Exports
#################################################################################
export generate_pras_system
export PRAS
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

"""
    PRAS

Module for Probabilistic Resource Adequacy Studies (PRAS).

Re-exported in PRASInterface

# Source

https://github.com/NREL/PRAS.jl
"""
module PRAS
using Reexport
const PRAS_VERSION = "v0.6.0"
include("MinCostFlows/MinCostFlows.jl")
include("PRASBase/PRASBase.jl")
include("ResourceAdequacy/ResourceAdequacy.jl")
include("CapacityCredit/CapacityCredit.jl")
end

import .PRAS.assess
import .PRAS.LOLE
import .PRAS.EUE
import .PRAS.val
import .PRAS.stderror
import .PRAS.SequentialMonteCarlo

import .PRAS.Shortfall
import .PRAS.ShortfallSamples
import .PRAS.Surplus
import .PRAS.SurplusSamples
import .PRAS.Flow
import .PRAS.FlowSamples
import .PRAS.Utilization
import .PRAS.UtilizationSamples
import .PRAS.StorageEnergy
import .PRAS.StorageEnergySamples
import .PRAS.GeneratorStorageEnergy
import .PRAS.GeneratorStorageEnergySamples
import .PRAS.GeneratorAvailability
import .PRAS.StorageAvailability
import .PRAS.GeneratorStorageAvailability
import .PRAS.LineAvailability

import .PRAS.SystemModel

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

end
