#######################################################
# Structs to parse and store the outage information
#######################################################
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

##############################################
# Converting FOR and MTTR to λ and μ
##############################################
function rate_to_probability(for_gen::Float64, mttr::Int64)

    # Check to make sure FOR data is not in % 
    if (for_gen > 1.0)
        for_gen = for_gen / 100
    end

    if (for_gen == 1.0)
        λ = 1.0
        μ = 0.0
    else
        if ~(mttr == 0)
            μ = 1 / mttr
        else
            μ = 1.0
        end
        λ = (μ * for_gen) / (1 - for_gen)
    end

    return (λ=λ, μ=μ)
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
