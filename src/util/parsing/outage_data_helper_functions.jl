#######################################################
# Structs to parse and store the outage information
#######################################################
struct outage_data
    prime_mover::String
    thermal_fuel::String
    capacity::Int64
    FOR::Float64
    MTTR::Int64

    outage_data(prime_mover  = "PrimeMovers.Default", thermal_fuel ="ThermalFuels.Default", capacity = 100, FOR=0.5,MTTR = 50) =new(prime_mover,thermal_fuel,capacity,FOR,MTTR)
end 

##############################################
# Converting FOR and MTTR to λ and μ
##############################################
function outage_to_rate(for_gen::Float64, mttr::Int64)

    # Check to make sure FOR data is not in % 
    if (for_gen > 1.0)
        for_gen = for_gen/100
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

    return (λ = λ, μ = μ)
end