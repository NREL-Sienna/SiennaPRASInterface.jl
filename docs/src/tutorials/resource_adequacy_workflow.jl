# # [Resource adequacy workflow](@id resource_adequacy_workflow)
#
# In this tutorial we will run a probabilistic resource adequacy study on a
# [`PowerSystems.System`](@extref) using `SiennaPRASInterface.jl`. We load an RTS-GMLC test case,
# inspect outage data, run a [`PRASCore.Simulations.SequentialMonteCarlo`](@extref) simulation,
# request multiple result specifications, inspect the translated PRAS
# [`PRASCore.Systems.SystemModel`](@extref), and customize
# device mappings with a [`RATemplate`](@ref).
#
# ## Prerequisites
#
# - Julia 1.8 or later
# - `SiennaPRASInterface`, `PowerSystems`, and `PowerSystemCaseBuilder` installed in your environment

using SiennaPRASInterface
using PowerSystems
using PowerSystemCaseBuilder

# ## Step 1: Load a test system
#
# We load the RTS-GMLC hourly data set with static outage rates already attached to generators
# from [`PowerSystemCaseBuilder`](https://sienna-platform.github.io/PowerSystemCaseBuilder.jl/stable/),
# using [`PowerSystemCaseBuilder.build_system`](@extref). We then call
# [`PowerSystems.set_units_base_system!`](@extref) to work in natural units (MW).

sys = build_system(SPISystems, "RTS_GMLC_Hourly with Static Outage Data")
set_units_base_system!(sys, UnitSystem.NATURAL_UNITS)
sys

# Notice the system summary lists generators, storage, and time series counts. Resource adequacy
# analysis requires at least one static time series on the system.

# ## Step 2: Inspect regions and outage data
#
# PRAS aggregates the network into regions. By default, [`assess`](@ref) uses
# [`PowerSystems.Area`](@extref) as the aggregation topology. We call
# [`PowerSystems.show_components`](@extref) to list the areas in the system:

show_components(Area, sys)

# Generators participating in the simulation should carry a
# [`PowerSystems.GeometricDistributionForcedOutage`](@extref) supplemental attribute.
# See [Outage and contingency data](@extref PowerSystems :label:`outage_and_contingency_data`)
# for background on stochastic forced outages. We use [`PowerSystems.get_component`](@extref PowerSystems :jl:method:`InfrastructureSystems.get_component-Union{Tuple{T}, Tuple{Type{T}, System, AbstractString}} where T<:Component`) and
# [`PowerSystems.get_supplemental_attributes`](@extref PowerSystems :jl:method:`InfrastructureSystems.get_supplemental_attributes-Union{Tuple{T}, Tuple{Function, Type{T}, System}} where T<:SupplementalAttribute`) to inspect a thermal generator:

gen = get_component(ThermalStandard, sys, "101_CT_1")

outage = get_supplemental_attributes(GeometricDistributionForcedOutage, gen)[1]
outage

# See that this generator already has mean time to recovery and outage transition probability
# values attached. If your system lacks outage data, see
# [Prepare stochastic outage data for resource adequacy](@ref prepare_stochastic_outage_data).

# ## Step 3: Configure the Monte Carlo simulation
#
# [`PRASCore.Simulations.SequentialMonteCarlo`](@extref) draws sequential outage states for each
# sample. We construct a simulation method with a modest number of samples so the tutorial runs
# quickly; increase `samples` for production studies:

sequential_monte_carlo =
    SequentialMonteCarlo(; samples=50, seed=1, threaded=false, verbose=false)
sequential_monte_carlo

# ## Step 4: Run a shortfall assessment
#
# [`assess`](@ref) translates the Sienna system to PRAS, runs the simulation, and returns one
# result object per requested [`PRASCore.Results.Shortfall`](@extref) specification. We pass
# [`PowerSystems.Area`](@extref) as the aggregation topology:

shortfalls, = assess(sys, Area, sequential_monte_carlo, Shortfall())
shortfalls

# Notice the return type is a PRAS shortfall result object that stores regional time series of
# unserved energy across Monte Carlo samples.

# ## Step 5: Compute reliability metrics
#
# [`PRASCore.Results.LOLE`](@extref) and [`PRASCore.Results.EUE`](@extref) summarize shortfall
# results into standard adequacy metrics. We extract numeric values with `val`:

lole = LOLE(shortfalls)
eue = EUE(shortfalls)
val(lole), val(eue)

# LOLE is reported in hours per study horizon;
# EUE is reported in megawatt-hours of expected unserved energy.

# ## Step 6: Request multiple result specifications
#
# A single [`assess`](@ref) call can request several
# [PRAS result specifications](@extref PRASCore :doc:`PRAS/results`) at
# once. We add [`PRASCore.Results.Surplus`](@extref) and
# [`PRASCore.Results.StorageEnergy`](@extref) alongside [`PRASCore.Results.Shortfall`](@extref):

shortfall, surplus, storage_energy = assess(
    sys,
    Area,
    sequential_monte_carlo,
    Shortfall(),
    Surplus(),
    StorageEnergy(),
)
typeof(shortfall), typeof(surplus), typeof(storage_energy)

# See that each specification returns its own result container. Use surplus results to study
# periods with excess generation and storage energy results to track state of charge statistics.

# ## Step 7: Inspect the translated PRAS system
#
# Before or instead of running a simulation, you can build the intermediate PRAS
# [`PRASCore.Systems.SystemModel`](@extref) with
# [`generate_pras_system`](@ref):

pras_system = generate_pras_system(sys, Area)
pras_system

# Notice the printed summary lists PRAS regions, generators, and interfaces derived from the
# Sienna system. To persist the model for use with native PRAS workflows, pass an `export_location`
# ending in `.pras` — see [Inspect or export a PRAS SystemModel](@ref inspect_pras_system).

# ## Step 8: Customize device mappings with a RATemplate
#
# The default translation lumps devices by type. To change how renewables are aggregated, we build
# a [`RATemplate`](@ref) and register a custom [`DeviceRAModel`](@ref) with [`set_device_model!`](@ref):

template = RATemplate(Area, deepcopy(SiennaPRASInterface.DEFAULT_DEVICE_MODELS))
set_device_model!(
    template,
    DeviceRAModel(
        RenewableDispatch,
        GeneratorPRAS;
        lump_renewable_generation=true,
    ),
)
template

# We re-run [`assess`](@ref) with the customized template and compute [`PRASCore.Results.LOLE`](@extref) from the
# shortfall result:

shortfalls_custom, = assess(sys, template, sequential_monte_carlo, Shortfall())
val(LOLE(shortfalls_custom))

# See that this LOLE value is lower than the default-area assessment from Step 5. Lumping renewable
# generation within a region can change adequacy statistics when intermittent resources are
# modeled at aggregate capacity rather than unit level.

# ## Next steps
#
# - [Prepare stochastic outage data for resource adequacy](@ref prepare_stochastic_outage_data) —
#   attach or update [`PowerSystems.GeometricDistributionForcedOutage`](@extref) data
# - [Configure device mappings for PRAS](@ref configure_device_mappings) — further customize
#   [`RATemplate`](@ref) translations
# - [PRAS 101 Walkthrough](@extref PRASCore :doc:`examples/pras_walkthrough`) —
#   native PRAS concepts and file formats
# - [Working with Time Series Data](@extref PowerSystems :doc:`tutorials/generated_working_with_time_series`) —
#   time-varying outage probabilities on a Sienna system
