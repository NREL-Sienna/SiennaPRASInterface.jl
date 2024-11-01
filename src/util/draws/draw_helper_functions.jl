##############################################
# Helpful Reference ; https://core.ac.uk/download/pdf/13643059.pdf
# Converting FOR and MTTR to λ and μ
##############################################
function rate_to_probability(for_gen::Float64, mttr::Int64)
    if (for_gen > 1.0)
        for_gen = for_gen / 100
    end

    if for_gen == 1.0
        return (λ=1.0, μ=0.0)  # can we error here instead?
    end
    if mttr != 0
        μ = 1 / mttr
    else # MTTR of 0.0 doesn't make much sense.
        μ = 1.0
    end
    return (λ=(μ * for_gen) / (1 - for_gen), μ=μ)
end

##############################################
# Initializing Availability
##############################################
function initialize_availability!(
    rng::Random.AbstractRNG,
    availability::Vector{Bool},
    nexttransition::Vector{Int},
    devices::Vector{PSY.Generator},
    t_last::Int,
)

    # TODO: When time series is available, we should get the time series
    all_device_λ = reduce(
        vcat,
        fill.(
            PSY.get_outage_transition_probability.(
                first.(
                    PSY.get_supplemental_attributes.(
                        PSY.GeometricDistributionForcedOutage,
                        devices,
                    )
                )
            ),
            1,
            t_last,
        ),
    )
    all_device_μ = reduce(
        vcat,
        fill.(
            1.0 ./
            PSY.get_mean_time_to_recovery.(
                first.(
                    PSY.get_supplemental_attributes.(
                        PSY.GeometricDistributionForcedOutage,
                        devices,
                    )
                )
            ),
            1,
            t_last,
        ),
    )

    for (i, device) in enumerate(devices)

        # TODO: When time series is available, we should get the time series
        transition_data = first(
            PSY.get_supplemental_attributes(PSY.GeometricDistributionForcedOutage, device),
        )
        λ = PSY.get_outage_transition_probability(transition_data) # TODO: When time series is available, we should get the first element of time series
        μ = 1 / PSY.get_mean_time_to_recovery(transition_data) # TODO: When time series is available, we should get the first element of time series

        online = rand(rng) < μ / (λ + μ)

        availability[i] = online

        transitionprobs = online ? all_device_λ : all_device_μ

        nexttransition[i] = randtransitiontime(rng, transitionprobs, i, 1, t_last)
    end

    return availability, nexttransition
end
##############################################
# Update Availability
##############################################

function update_availability!(
    rng::Random.AbstractRNG,
    availability::Vector{Bool},
    nexttransition::Vector{Int},
    devices::Vector{PSY.Generator},
    t_now::Int,
    t_last::Int,
)

    # TODO: When time series is available, we should get the time series
    all_device_λ = reduce(
        vcat,
        fill.(
            PSY.get_outage_transition_probability.(
                first.(
                    PSY.get_supplemental_attributes.(
                        PSY.GeometricDistributionForcedOutage,
                        devices,
                    )
                )
            ),
            1,
            t_last,
        ),
    )
    all_device_μ = reduce(
        vcat,
        fill.(
            1.0 ./
            PSY.get_mean_time_to_recovery.(
                first.(
                    PSY.get_supplemental_attributes.(
                        PSY.GeometricDistributionForcedOutage,
                        devices,
                    )
                )
            ),
            1,
            t_last,
        ),
    )

    for i in 1:length(devices)
        if nexttransition[i] == t_now # Unit switches states
            transitionprobs = (availability[i] ⊻= true) ? all_device_λ : all_device_μ
            nexttransition[i] = randtransitiontime(rng, transitionprobs, i, t_now, t_last)
        end
    end

    return availability, nexttransition
end

function randtransitiontime(
    rng::Random.AbstractRNG,
    p::Matrix{Float64},
    i::Int,
    t_now::Int,
    t_last::Int,
)
    cdf = 0.0
    p_noprevtransition = 1.0

    x = rand(rng)
    t = t_now + 1

    while t <= t_last
        p_it = p[i, t]
        cdf += p_noprevtransition * p_it
        x < cdf && return t
        p_noprevtransition *= (1.0 - p_it)
        t += 1
    end

    return t_last + 1
end

# Initial time, resolution and length
function get_timestamps(
    initial_time::Dates.DateTime,
    resolution::TIMEPERIOD,
    steps::Int,
    horizon::Int,
) where {TIMEPERIOD <: Dates.TimePeriod}
    finish_time = initial_time + (resolution * ((steps * horizon) - 1))
    return collect(StepRange(initial_time, resolution, finish_time))
end
