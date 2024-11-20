"""
    outage_data

Struct to store the outage information for a component.
"""
mutable struct outage_data
    prime_mover::PSY.PrimeMovers
    fuel::Union{PSY.ThermalFuels, Missing}
    max_capacity::Int64
    FOR::Float64
    MTTR::Int64
    min_capacity::Int64

    outage_data(
        prime_mover=PSY.PrimeMovers.OT,
        fuel=PSY.ThermalFuels.OTHER,
        max_capacity=100,
        FOR=0.5,
        MTTR=50,
        min_capacity=0,
    ) = new(prime_mover, fuel, max_capacity, FOR, MTTR, min_capacity)
end

# Defined on outage_data for Base.sort! to work
Base.isless(data_1::outage_data, data_2::outage_data) =
    data_1.max_capacity < data_2.max_capacity

function add_min_capacity!(outage_values::Vector{outage_data})
    for (p_m, fuel) in
        unique(zip(getfield.(outage_values, :prime_mover), getfield.(outage_values, :fuel)))
        filtered_data = if (ismissing(fuel))
            filter(x -> x.prime_mover == p_m, outage_values)
        else
            filter(x -> (x.prime_mover == p_m && x.fuel == fuel), outage_values)
        end

        sort!(filtered_data)
        for (i, data) in enumerate(filtered_data)
            if (i == 1)
                continue
            else
                data.min_capacity = filtered_data[i - 1].max_capacity
            end
        end
    end
end
