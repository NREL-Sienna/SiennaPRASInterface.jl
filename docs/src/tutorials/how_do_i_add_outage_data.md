# How do I add outage data to Sienna?

You can attach outage data to `PowerSystems` `Components` by using the
supplemental attribute [`GeometricDistributionForcedOutage`](https://nrel-sienna.github.io/PowerSystems.jl/stable/api/public/#PowerSystems.GeometricDistributionForcedOutage).

## Step 1 : Parse your outage data into Sienna

`SiennaPRASInterface.jl` uses outage information in the form of independent `mean_time_to_recovery`
in units of hours and `outage_transition_probability` in probability of outage per hour.
A simple Markov model models the transitions between out and active using these parameters.

We support data either being fixed and specified in the `GeometricDistributionForcedOutage` object
or attached as time-series to the `GeometricDistributionForcedOutage` struct.

### Creating a `GeometricDistributionForcedOutage` from fixed data

```julia
using PowerSystems
transition_data = GeometricDistributionForcedOutage(;
    mean_time_to_recovery=10,  # Units of hours
    outage_transition_probability=0.005,  # Probability for outage per hour
)
```

### Creating a `GeometricDistributionForcedOutage` from time series data

Time series should be attached to a `GeometricDistributionForcedOutage` object
under the keys `recovery_probability` (1/`mean_time_to_recovery`) and `outage_probability`.

See the [Sienna time-series documentation on working with time-series](https://nrel-sienna.github.io/PowerSystems.jl/stable/tutorials/working_with_time_series/).

```julia
using PowerSystems
using Dates
using TimeSeries

transition_data = GeometricDistributionForcedOutage(;
    mean_time_to_recovery=10,  # Units of hours
    outage_transition_probability=0.005,  # Probability for outage per hour
)

outage_probability = [0.1, 0.1, 0.2, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4]
recovery_probability = [0.1, 0.1, 0.2, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4]

# Your resolution and length must match the other SingleTimeSeries in your System.
resolution = Dates.Minute(5)
timestamps = range(DateTime("2020-01-01T08:00:00"); step = resolution, length = 24)
outage_timearray = TimeArray(timestamps, outage_probability)
outage_time_series = SingleTimeSeries(;
    name = "outage_probability",
    data = outage_timearray,
)

recovery_timearray = TimeArray(timestamps, recovery_probability)
recovery_time_series = SingleTimeSeries(;
    name = "recovery_probability",
    data = recovery_timearray,
)

# Here we assume you have a system named sys
PSY.add_time_series!(sys, transition_data, outage_time_series)
PSY.add_time_series!(sys, transition_data, recovery_time_series)
```

## Step 2 : Attaching Data to Components

Once you have a `GeometricDistributionForcedOutage` object, then you can add it to
any components with that data:

```julia
component = get_component(Generator, sys, "test_generator")
add_supplemental_attribute!(sys, component, transition_data)
```

## Step 3 : Run simulations and verify result

```julia
using SiennaPRASInterface
method = SequentialMonteCarlo(samples=10_000, seed=1)
shortfalls, = assess(sys, PowerSystems.Area, sequential_monte_carlo, Shortfall())
eue = EUE(shortfalls)
```