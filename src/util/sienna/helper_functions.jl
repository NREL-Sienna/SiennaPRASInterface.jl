#######################################################
# Function to get available components in AggregationTopology
#######################################################
function get_available_components_in_aggregation_topology(type::Type{<:PSY.StaticInjection}, sys::PSY.System, region::PSY.AggregationTopology)
    avail_comps =  filter(x ->(PSY.get_available(x)),collect(PSY.get_components_in_aggregation_topology(type, sys, region)))
    return avail_comps
end

#######################################################
# Helper Functions
# Generators
#######################################################
#######################################################
# Functions to get generator category
#######################################################
function get_generator_category(gen::GEN) where {GEN <: PSY.RenewableGen}
    return string(PSY.get_prime_mover_type(gen))
end

function get_generator_category(gen::GEN) where {GEN <: PSY.ThermalGen}
    return string(PSY.get_fuel(gen))
end

function get_generator_category(gen::GEN) where {GEN <: PSY.HydroGen}
    return "Hydro"
end

function get_generator_category(stor::GEN) where {GEN <: PSY.Storage}
    if (occursin("Distributed",PSY.get_name(stor)))
        return "Distributed_Storage"
    elseif (occursin("Battery",PSY.get_name(stor)))
        return "Battery_Storage"
    else
        return "Battery"
    end
end

function get_generator_category(stor::GEN) where {GEN <: PSY.HybridSystem}
    return "Hybrid-System"
end
#######################################################
# Helper Functions
# Lines
#######################################################
#######################################################
# Line Rating
#######################################################
function line_rating(line::Union{PSY.Line,PSY.MonitoredLine})
    rate = PSY.get_rate(line);
    return(forward_capacity = abs(rate) , backward_capacity = abs(rate))
end

function line_rating(line::PSY.TwoTerminalHVDCLine)
    forward_capacity = getfield(PSY.get_active_power_limits_from(line), :max)
    backward_capacity = getfield(PSY.get_active_power_limits_to(line), :max)
    return(forward_capacity = abs(forward_capacity), backward_capacity = abs(backward_capacity))
end

function line_rating(line::DCLine) where {DCLine<:HVDCLineTypes}
    error("line_rating isn't defined for $(typeof(line))")
end


