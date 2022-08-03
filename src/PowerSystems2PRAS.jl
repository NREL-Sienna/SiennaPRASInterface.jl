#######################################################
# Surya
# NREL
# January 2021
# SIIP --> PRAS Linkage Package
#######################################################
module PowerSystems2PRAS
#################################################################################
# Exports
#################################################################################
export make_pras_system
export make_ercot_pras_system
export make_RT_ercot_pras_system

export generate_outage_profile
export generate_csv_outage_profile
export extract_csv_outage_profile

export export_pras_system
export export_pros_system
export load_pros_system
export reduce_pras_system_size
export collapse_pras_regions
export extract_regions
export add_csv_time_series!
export add_csv_time_series_single_stage!
#################################################################################
# Imports
#################################################################################
import PowerSystems
import PRAS
import PrOS

import Dates
import TimeSeries
import HDF5
import TimeZones
import DataFrames
import CSV
import YAML
import JLD
import JSON

import UUIDs
import Statistics
#################################################################################
# Includes
#################################################################################
include("parsers/power_system_table_data.jl") # Over-writes some PSY functions.
include("main/PSY2PRAS.jl")
include("main/PRAS2PSY.jl")
include("main/PSY2PrOS.jl") 
include("main/ERCOT-PSY2PrOS.jl")
include("main/ERCOT-RT-PSY2PrOS.jl")
include("main/PrOS2PSY.jl")
include("util/save-pras-system.jl")
include("util/save-pros-system.jl")
include("util/pras_system_size_reducer.jl")
include("util/collapse_pras_regions.jl")
include("util/extract_regions.jl")
include("util/add_csv_time_series_data.jl")
end