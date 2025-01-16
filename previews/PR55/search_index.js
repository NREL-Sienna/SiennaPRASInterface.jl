var documenterSearchIndex = {"docs":
[{"location":"tutorials/intro_page/#Examples","page":"Tutorials","title":"Examples","text":"","category":"section"},{"location":"tutorials/intro_page/","page":"Tutorials","title":"Tutorials","text":"Tutorials to be created soon.","category":"page"},{"location":"api/internal/#Internal-API","page":"Internal API Reference","title":"Internal API","text":"","category":"section"},{"location":"api/internal/","page":"Internal API Reference","title":"Internal API Reference","text":"Modules = [SiennaPRASInterface, PRAS]\nPublic = false","category":"page"},{"location":"#SiennaPRASInterface.jl","page":"Welcome Page","title":"SiennaPRASInterface.jl","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"CurrentModule = SiennaPRASInterface","category":"page"},{"location":"#About","page":"Welcome Page","title":"About","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"SiennaPRASInterface.jl is a Julia package that provides an interface to PRAS.jl from Sienna's PowerSystem.jl's System data model.","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"The Probabilistic Resource Adequacy Suite (PRAS) analyzes the resource adequacy of a bulk power system using Monte Carlo methods.","category":"page"},{"location":"#Getting-Started","page":"Welcome Page","title":"Getting Started","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"To use SiennaPRASInterface.jl, you first need a System from PowerSystems.jl","category":"page"},{"location":"#.-Install","page":"Welcome Page","title":"1. Install","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"] add SiennaPRASInterface","category":"page"},{"location":"#.-Add-Data","page":"Welcome Page","title":"2. Add Data","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"Add outage information to generators using the supplemental attribute GeometricDistributionForcedOutage.","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"using PowerSystems\ntransition_data = GeometricDistributionForcedOutage(;\n    mean_time_to_recovery=10,  # Units of hours\n    outage_transition_probability=0.005,  # Probability for outage per hour\n)\ncomponent = get_component(Generator, sys, \"test_generator\")\nadd_supplemental_attribute!(sys, component, transition_data)","category":"page"},{"location":"#.-Calculate-Shortfalls-and-Expected-Unserved-Energy-on-System","page":"Welcome Page","title":"3. Calculate Shortfalls and Expected Unserved Energy on System","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"using SiennaPRASInterface\nmethod = SequentialMonteCarlo(samples=10_000, seed=1)\nshortfalls, = assess(sys, PowerSystems.Area, sequential_monte_carlo, Shortfall())\neue = EUE(shortfalls)","category":"page"},{"location":"#Documentation","page":"Welcome Page","title":"Documentation","text":"","category":"section"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"PRAS Documentation","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"Pages = [\"api/public.md\", \"tutorials\"]\nDepth = 2","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"","category":"page"},{"location":"","page":"Welcome Page","title":"Welcome Page","text":"SiennaPRASInterface has been developed as part of the Transmission Planning Tools Maintenance project at the U.S. Department of Energy's National Renewable Energy Laboratory (NREL) funded by DOE Grid Deployment Office (GDO).","category":"page"},{"location":"tutorials/how_do_i_add_outage_data/#How-do-I-add-outage-data-to-Sienna?","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"","category":"section"},{"location":"tutorials/how_do_i_add_outage_data/","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"You can attach outage data to PowerSystems Components by using the supplemental attribute GeometricDistributionForcedOutage.","category":"page"},{"location":"tutorials/how_do_i_add_outage_data/#Step-1-:-Parse-your-outage-data-into-Sienna","page":"How do I add outage data to Sienna?","title":"Step 1 : Parse your outage data into Sienna","text":"","category":"section"},{"location":"tutorials/how_do_i_add_outage_data/","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"SiennaPRASInterface.jl uses outage information in the form of independent mean_time_to_recovery in units of hours and outage_transition_probability in probability of outage per hour. A simple Markov model models the transitions between out and active using these parameters.","category":"page"},{"location":"tutorials/how_do_i_add_outage_data/","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"We support data either being fixed and specified in the GeometricDistributionForcedOutage object or attached as time-series to the GeometricDistributionForcedOutage struct.","category":"page"},{"location":"tutorials/how_do_i_add_outage_data/#Creating-a-GeometricDistributionForcedOutage-from-fixed-data","page":"How do I add outage data to Sienna?","title":"Creating a GeometricDistributionForcedOutage from fixed data","text":"","category":"section"},{"location":"tutorials/how_do_i_add_outage_data/","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"using PowerSystems\ntransition_data = GeometricDistributionForcedOutage(;\n    mean_time_to_recovery=10,  # Units of hours\n    outage_transition_probability=0.005,  # Probability for outage per hour\n)","category":"page"},{"location":"tutorials/how_do_i_add_outage_data/#Creating-a-GeometricDistributionForcedOutage-from-time-series-data","page":"How do I add outage data to Sienna?","title":"Creating a GeometricDistributionForcedOutage from time series data","text":"","category":"section"},{"location":"tutorials/how_do_i_add_outage_data/","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"Time series should be attached to a GeometricDistributionForcedOutage object under the keys recovery_probability (1/mean_time_to_recovery) and outage_probability.","category":"page"},{"location":"tutorials/how_do_i_add_outage_data/","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"See the Sienna time-series documentation on working with time-series.","category":"page"},{"location":"tutorials/how_do_i_add_outage_data/","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"using PowerSystems\nusing Dates\nusing TimeSeries\n\ntransition_data = GeometricDistributionForcedOutage(;\n    mean_time_to_recovery=10,  # Units of hours\n    outage_transition_probability=0.005,  # Probability for outage per hour\n)\n\noutage_probability = [0.1, 0.1, 0.2, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4]\nrecovery_probability = [0.1, 0.1, 0.2, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.4]\n\n# Your resolution and length must match the other SingleTimeSeries in your System.\nresolution = Dates.Minute(5)\ntimestamps = range(DateTime(\"2020-01-01T08:00:00\"); step=resolution, length=24)\noutage_timearray = TimeArray(timestamps, outage_probability)\noutage_time_series = SingleTimeSeries(; name=\"outage_probability\", data=outage_timearray)\n\nrecovery_timearray = TimeArray(timestamps, recovery_probability)\nrecovery_time_series =\n    SingleTimeSeries(; name=\"recovery_probability\", data=recovery_timearray)\n\n# Here we assume you have a system named sys\nPSY.add_time_series!(sys, transition_data, outage_time_series)\nPSY.add_time_series!(sys, transition_data, recovery_time_series)","category":"page"},{"location":"tutorials/how_do_i_add_outage_data/#Step-2-:-Attaching-Data-to-Components","page":"How do I add outage data to Sienna?","title":"Step 2 : Attaching Data to Components","text":"","category":"section"},{"location":"tutorials/how_do_i_add_outage_data/","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"Once you have a GeometricDistributionForcedOutage object, then you can add it to any components with that data:","category":"page"},{"location":"tutorials/how_do_i_add_outage_data/","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"component = get_component(Generator, sys, \"test_generator\")\nadd_supplemental_attribute!(sys, component, transition_data)","category":"page"},{"location":"tutorials/how_do_i_add_outage_data/#Step-3-:-Run-simulations-and-verify-result","page":"How do I add outage data to Sienna?","title":"Step 3 : Run simulations and verify result","text":"","category":"section"},{"location":"tutorials/how_do_i_add_outage_data/","page":"How do I add outage data to Sienna?","title":"How do I add outage data to Sienna?","text":"using SiennaPRASInterface\nmethod = SequentialMonteCarlo(samples=10_000, seed=1)\nshortfalls, = assess(sys, PowerSystems.Area, sequential_monte_carlo, Shortfall())\neue = EUE(shortfalls)","category":"page"},{"location":"api/public/#Public-API-Reference","page":"Public API Reference","title":"Public API Reference","text":"","category":"section"},{"location":"api/public/","page":"Public API Reference","title":"Public API Reference","text":"SiennaPRASInterface\ngenerate_pras_system\nSystemModel\nassess\nSequentialMonteCarlo\nShortfall\nSurplus\nFlow\nUtilization\nStorageEnergy\nGeneratorStorageEnergy\nLOLE\nEUE\nval\nstderror\ngenerate_outage_profile\ngenerate_csv_outage_profile\nadd_csv_time_series!\nadd_csv_time_series_single_stage!\nmake_generator_outage_draws!\nShortfallSamples\nSurplusSamples\nFlowSamples\nUtilizationSamples\nStorageEnergySamples\nGeneratorStorageEnergySamples\nGeneratorAvailability\nGeneratorStorageAvailability\nLineAvailability","category":"page"},{"location":"api/public/#SiennaPRASInterface","page":"Public API Reference","title":"SiennaPRASInterface","text":"PowerSystems Interface for Probabilistic Resource Adequacy Studies (PRAS)\n\nKey Functions\n\ngenerate_pras_system: convert PSY to PRAS model\nassess: assess PRAS model\n\nKey PRAS Types\n\nSystemModel: PRAS data structure\nSequentialMonteCarlo: method for PRAS analysis\nShortfall: PRAS metric for missing generation\nLOLE: PRAS metric for loss of load expectation\nEUE: PRAS metric for energy unserved expectation\n\n\n\n\n\n","category":"module"},{"location":"api/public/#SiennaPRASInterface.generate_pras_system","page":"Public API Reference","title":"SiennaPRASInterface.generate_pras_system","text":"generate_pras_system(sys::PSY.System, aggregation; kwargs...)\n\nSienna/Data PowerSystems.jl System is the input and an object of PRAS SystemModel is returned. ...\n\nArguments\n\nsys::PSY.System: Sienna/Data PowerSystems.jl System\naggregation<:PSY.AggregationTopology: \"PSY.Area\" (or) \"PSY.LoadZone\" {Optional}\nlump_region_renewable_gens::Bool: Whether to lumps PV and Wind generators in a region because usually these generators don't have FOR data {Optional}\nexport_location::String: Export location of the .pras file ...\n\nReturns\n\n- `PRASCore.SystemModel`: PRAS SystemModel object\n\nExamples\n\njulia> generate_pras_system(psy_sys)\nPRAS SystemModel\n\n\n\n\n\ngenerate_pras_system(sys_location::String, aggregation; kwargs...)\n\nGenerate a PRAS SystemModel from a Sienna/Data PowerSystems System JSON file.\n\nArguments\n\nsys_location::String: Location of the Sienna/Data PowerSystems System JSON file\naggregation::Type{AT}: Aggregation topology type\nlump_region_renewable_gens::Bool: Lumping of region renewable generators\nexport_location::Union{Nothing, String}: Export location of the .pras file\n\nReturns\n\nPRASCore.SystemModel: PRAS SystemModel\n\n\n\n\n\n","category":"function"},{"location":"api/public/#PRASCore.Systems.SystemModel","page":"Public API Reference","title":"PRASCore.Systems.SystemModel","text":"SystemModel\n\nA SystemModel contains a representation of a power system to be studied with PRAS.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Simulations.assess","page":"Public API Reference","title":"PRASCore.Simulations.assess","text":"assess(system::SystemModel, method::SequentialMonteCarlo, resultspecs::ResultSpec...)\n\nRun a Sequential Monte Carlo simulation on a system using the method data and return resultspecs.\n\nArguments\n\nsystem::SystemModel: PRAS data structure\nmethod::SequentialMonteCarlo: method for PRAS analysis\nresultspecs::ResultSpec...: PRAS metric for metrics like Shortfall missing generation\n\nReturns\n\nresults::Tuple{Vararg{ResultAccumulator{SequentialMonteCarlo}}}: PRAS metric results\n\n\n\n\n\nassess(\n    sys::PSY.System,\n    aggregation::Type{AT},\n    method::PRASCore.SequentialMonteCarlo,\n    resultsspecs::PRASCore.Results.ResultSpec...,\n) where {AT <: PSY.AggregationTopology}\n\nEstimate resource adequacy using Monte Carlo simulation.\n\nArguments\n\nsys::PSY.System: PowerSystems.jl system model\naggregation::Type{AT}: Aggregation topology to use in translating to PRAS\nmethod::PRASCore.SequentialMonteCarlo: Simulation method to use\nresultsspec::PRASCore.Results.ResultSpec...: Results to compute\n\nReturns\n\nTuple of results from resultsspec: default is (ShortfallResult,)\n\n\n\n\n\n","category":"function"},{"location":"api/public/#PRASCore.Simulations.SequentialMonteCarlo","page":"Public API Reference","title":"PRASCore.Simulations.SequentialMonteCarlo","text":"SequentialMonteCarlo(;\n    samples::Int=10_000,\n    seed::Integer=rand(UInt64),\n    verbose::Bool=false,\n    threaded::Bool=true\n)\n\nSequential Monte Carlo simulation parameters for PRAS analysis\n\nIt it recommended that you fix the random seed for reproducibility.\n\nArguments\n\nsamples::Int=10_000: Number of samples\nseed::Integer=rand(UInt64): Random seed\nverbose::Bool=false: Print progress\nthreaded::Bool=true: Use multi-threading\n\nReturns\n\nSequentialMonteCarlo: PRAS simulation specification\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.Shortfall","page":"Public API Reference","title":"PRASCore.Results.Shortfall","text":"Shortfall\n\nThe Shortfall result specification reports expectation-based resource adequacy risk metrics such as EUE and LOLE, producing a ShortfallResult.\n\nA ShortfallResult can be directly indexed by a region name and a timestamp to retrieve a tuple of sample mean and standard deviation, estimating  the average unserved energy in that region and timestep. However, in most cases it's simpler to use EUE and LOLE constructors to directly retrieve standard risk metrics.\n\nExample:\n\nshortfall, =\n    assess(sys, SequentialMonteCarlo(samples=1000), Shortfall())\n\nperiod = ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")\n\n# Unserved energy mean and standard deviation\nsf_mean, sf_std = shortfall[\"Region A\", period]\n\n# System-wide risk metrics\neue = EUE(shortfall)\nlole = LOLE(shortfall)\n\n# Regional risk metrics\nregional_eue = EUE(shortfall, \"Region A\")\nregional_lole = LOLE(shortfall, \"Region A\")\n\n# Period-specific risk metrics\nperiod_eue = EUE(shortfall, period)\nperiod_lolp = LOLE(shortfall, period)\n\n# Region- and period-specific risk metrics\nperiod_eue = EUE(shortfall, \"Region A\", period)\nperiod_lolp = LOLE(shortfall, \"Region A\", period)\n\nSee ShortfallSamples for recording sample-level shortfall results.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.Surplus","page":"Public API Reference","title":"PRASCore.Results.Surplus","text":"Surplus\n\nThe Surplus result specification reports unused generation and storage discharge capability of Regions, producing a SurplusResult.\n\nA SurplusResult can be indexed by region name and timestamp to retrieve a tuple of sample mean and standard deviation, estimating the average unused capacity in that region and timestep.\n\nExample:\n\nsurplus, =\n    assess(sys, SequentialMonteCarlo(samples=1000), Surplus())\n\nsurplus_mean, surplus_std =\n    surplus[\"Region A\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\nSee SurplusSamples for sample-level surplus results.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.Flow","page":"Public API Reference","title":"PRASCore.Results.Flow","text":"Flow\n\nThe Flow result specification reports the estimated average flow across transmission Interfaces, producing a FlowResult.\n\nA FlowResult can be indexed by a directional Pair of region names and a timestamp to retrieve a tuple of sample mean and standard deviation, estimating the average net flow magnitude and direction relative to the given directed interface in that timestep. For a query of \"Region A\" => \"Region B\", if estimated average flow was from A to B, the reported value would be positive, while if average flow was in the reverse direction, from B to A, the value would be negative.\n\nExample:\n\nflows, =\n    assess(sys, SequentialMonteCarlo(samples=1000), Flow())\n\nflow_mean, flow_std =\n    flows[\"Region A\" => \"Region B\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\nflow2_mean, flow2_std =\n    flows[\"Region B\" => \"Region A\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n@assert flow_mean == -flow2_mean\n\nSee FlowSamples for sample-level flow results.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.Utilization","page":"Public API Reference","title":"PRASCore.Results.Utilization","text":"Utilization\n\nThe Utilization result specification reports the estimated average absolute utilization of Interfaces, producing a UtilizationResult.\n\nWhereas Flow reports the average directional power transfer across an interface, Utilization reports the absolute value of flow relative to the interface's transfer capability (counting the effects of line outages). For example, a symmetrically-constrained interface which is fully congested with max power flowing in one direction in half of the samples, and the other direction in the remaining samples, would have an average flow of 0 MW, but an average utilization of 100%.\n\nA UtilizationResult can be indexed by a Pair of region names and a timestamp to retrieve a tuple of sample mean and standard deviation, estimating the average utilization of the interface. Given the absolute value nature of the outcome, results are independent of direction. Querying \"Region A\" => \"Region B\" will yield the same result as \"Region B\" => \"Region A\".\n\nExample:\n\nutils, =\n    assess(sys, SequentialMonteCarlo(samples=1000), Utilization())\n\nutil_mean, util_std =\n    utils[\"Region A\" => \"Region B\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\nutil2_mean, util2_std =\n    utils[\"Region B\" => \"Region A\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert util_mean == util2_mean\n\nSee UtilizationSamples for sample-level utilization results.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.StorageEnergy","page":"Public API Reference","title":"PRASCore.Results.StorageEnergy","text":"StorageEnergy\n\nThe StorageEnergy result specification reports the average state of charge of Storages, producing a StorageEnergyResult.\n\nA StorageEnergyResult can be indexed by storage device name and a timestamp to retrieve a tuple of sample mean and standard deviation, estimating the average energy level for the given storage device in that timestep.\n\nExample:\n\nstorenergy, =\n    assess(sys, SequentialMonteCarlo(samples=1000), StorageEnergy())\n\nsoc_mean, soc_std =\n    storenergy[\"MyStorage123\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\nSee StorageEnergySamples for sample-level storage states of charge.\n\nSee GeneratorStorageEnergy for average generator-storage states of charge.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.GeneratorStorageEnergy","page":"Public API Reference","title":"PRASCore.Results.GeneratorStorageEnergy","text":"GeneratorStorageEnergy\n\nThe GeneratorStorageEnergy result specification reports the average state of charge of GeneratorStorages, producing a GeneratorStorageEnergyResult.\n\nA GeneratorStorageEnergyResult can be indexed by generator-storage device name and a timestamp to retrieve a tuple of sample mean and standard deviation, estimating the average energy level for the given generator-storage device in that timestep.\n\nExample:\n\ngenstorenergy, =\n    assess(sys, SequentialMonteCarlo(samples=1000), GeneratorStorageEnergy())\n\nsoc_mean, soc_std =\n    genstorenergy[\"MyGeneratorStorage123\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\nSee GeneratorStorageEnergySamples for sample-level generator-storage states of charge.\n\nSee StorageEnergy for average storage states of charge.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.LOLE","page":"Public API Reference","title":"PRASCore.Results.LOLE","text":"LOLE\n\nLOLE reports loss of load expectation over a particular time period and regional extent. When the reporting period is a single simulation timestep, the metric is equivalent to loss of load probability (LOLP).\n\nContains both the estimated value itself as well as the standard error of that estimate, which can be extracted with val and stderror, respectively.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.EUE","page":"Public API Reference","title":"PRASCore.Results.EUE","text":"EUE\n\nEUE reports expected unserved energy over a particular time period and regional extent.\n\nContains both the estimated value itself as well as the standard error of that estimate, which can be extracted with val and stderror, respectively.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#SiennaPRASInterface.generate_outage_profile","page":"Public API Reference","title":"SiennaPRASInterface.generate_outage_profile","text":"generate_outage_profile(pras_system,num_runs,psy_sys,num_scenarios,location)\n\nProcess the assess results to get timeseries of generator status and include this timeseries data to the corresponding component in PSY System and exported using to_json method (serializing the PSY System).\n\n...\n\nArguments\n\npras_system::PRASCore.SystemModel: PRAS System\nnum_runs::Int64: Number of PRAS runs\npsy_sys::PSY.System: PSY System\nnum_scenarios::Int64: Number of scenarios of user interest.\nlocation::String: Location to store outage profile. ...\n\nExamples\n\njulia> generate_outage_profile(results, pras_sys, psy_sys, 1)\nPSY System exported using to_json method in InfrastructureSystems\n\n\n\n\n\n","category":"function"},{"location":"api/public/#SiennaPRASInterface.generate_csv_outage_profile","page":"Public API Reference","title":"SiennaPRASInterface.generate_csv_outage_profile","text":"generate_outage_profile(pras_system,num_runs,psy_sys,num_scenarios,location)\n\nProcess the assess results to get timeseries of generator status and include this timeseries data to the corresponding component in PSY System and exported using to_json method (serializing the PSY System).\n\n...\n\nArguments\n\npras_system::PRASCore.SystemModel: PRAS System\nnum_runs::Int64: Number of PRAS runs\npsy_sys::PSY.System: PSY System\nnum_scenarios::Int64: Number of scenarios of user interest.\nlocation::String: Location to store outage profile. ...\n\nExamples\n\njulia> generate_outage_profile(results, pras_sys, psy_sys, 1)\nPSY System exported using to_json method in InfrastructureSystems\n\n\n\n\n\n","category":"function"},{"location":"api/public/#SiennaPRASInterface.add_csv_time_series!","page":"Public API Reference","title":"SiennaPRASInterface.add_csv_time_series!","text":"add_csv_time_series!(\n    sys_DA,\n    sys_RT,\n    outage_csv_location::String;\n    days_of_interest::Union{Nothing, UnitRange}=nothing,\n    add_scenario::Union{Nothing, Int}=nothing,\n)\n\nGenerates outage profile for two stage PowerSimulation and adds availability time series data to Generators in PSY System from CSV files.\n\n\n\n\n\n","category":"function"},{"location":"api/public/#SiennaPRASInterface.add_csv_time_series_single_stage!","page":"Public API Reference","title":"SiennaPRASInterface.add_csv_time_series_single_stage!","text":"add_csv_time_series_single_stage!(\n    sys_DA,\n    outage_csv_location::String;\n    days_of_interest::Union{Nothing, UnitRange}=nothing,\n    add_scenario::Union{Nothing, Int}=nothing,\n)\n\nGenerates outage profile for single stage PowerSimulation and adds availability time series.\n\n\n\n\n\n","category":"function"},{"location":"api/public/#SiennaPRASInterface.make_generator_outage_draws!","page":"Public API Reference","title":"SiennaPRASInterface.make_generator_outage_draws!","text":"make_generator_outage_draws!(\n    sys,\n    initial_time::Dates.DateTime=nothing,\n    resolution::TIMEPERIOD=nothing,\n    steps::Int=nothing,\n    horizon::Int=nothing,\n) where {TIMEPERIOD <: Dates.TimePeriod}\n\nAdds availability time series to the generators in the system.\n\nMain function to make generator outage draws.\n\n\n\n\n\n","category":"function"},{"location":"api/public/#PRASCore.Results.ShortfallSamples","page":"Public API Reference","title":"PRASCore.Results.ShortfallSamples","text":"ShortfallSamples\n\nThe ShortfallSamples result specification reports sample-level unserved energy outcomes, producing a ShortfallSamplesResult.\n\nA ShortfallSamplesResult can be directly indexed by a region name and a timestamp to retrieve a vector of sample-level unserved energy results in that region and timestep. EUE and LOLE constructors can also be used to retrieve standard risk metrics.\n\nExample:\n\nshortfall, =\n    assess(sys, SequentialMonteCarlo(samples=10), ShortfallSamples())\n\nperiod = ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")\n\nsamples = shortfall[\"Region A\", period]\n\n@assert samples isa Vector{Float64}\n@assert length(samples) == 10\n\n# System-wide risk metrics\neue = EUE(shortfall)\nlole = LOLE(shortfall)\n\n# Regional risk metrics\nregional_eue = EUE(shortfall, \"Region A\")\nregional_lole = LOLE(shortfall, \"Region A\")\n\n# Period-specific risk metrics\nperiod_eue = EUE(shortfall, period)\nperiod_lolp = LOLE(shortfall, period)\n\n# Region- and period-specific risk metrics\nperiod_eue = EUE(shortfall, \"Region A\", period)\nperiod_lolp = LOLE(shortfall, \"Region A\", period)\n\nNote that this result specification requires large amounts of memory for larger sample sizes. See Shortfall for average shortfall outcomes when sample-level granularity isn't required.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.SurplusSamples","page":"Public API Reference","title":"PRASCore.Results.SurplusSamples","text":"SurplusSamples\n\nThe SurplusSamples result specification reports sample-level unused generation and storage discharge capability of Regions, producing a SurplusSamplesResult.\n\nA SurplusSamplesResult can be indexed by region name and timestamp to retrieve a vector of sample-level surplus values in that region and timestep.\n\nExample:\n\nsurplus, =\n    assess(sys, SequentialMonteCarlo(samples=10), SurplusSamples())\n\nsamples = surplus[\"Region A\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert samples isa Vector{Float64}\n@assert length(samples) == 10\n\nNote that this result specification requires large amounts of memory for larger sample sizes. See Surplus for estimated average surplus values when sample-level granularity isn't required.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.FlowSamples","page":"Public API Reference","title":"PRASCore.Results.FlowSamples","text":"FlowSamples\n\nThe FlowSamples result specification reports the sample-level magnitude and direction of power flows across Interfaces, producing a FlowSamplesResult.\n\nA FlowSamplesResult can be indexed by a directional Pair of region names and a timestamp to retrieve a vector of sample-level net flow magnitudes and directions relative to the given directed interface in that timestep. For a query of \"Region A\" => \"Region B\", if flow in one sample was from A to B, the reported value would be positive, while if flow was in the reverse direction, from B to A, the value would be negative.\n\nExample:\n\nflows, =\n    assess(sys, SequentialMonteCarlo(samples=10), FlowSamples())\n\nsamples = flows[\"Region A\" => \"Region B\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert samples isa Vector{Float64}\n@assert length(samples) == 10\n\nsamples2 = flows[\"Region B\" => \"Region A\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert samples == -samples2\n\nNote that this result specification requires large amounts of memory for larger sample sizes. See Flow for estimated average flow results when sample-level granularity isn't required.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.UtilizationSamples","page":"Public API Reference","title":"PRASCore.Results.UtilizationSamples","text":"UtilizationSamples\n\nThe UtilizationSamples result specification reports the sample-level absolute utilization of Interfaces, producing a UtilizationSamplesResult.\n\nWhereas FlowSamples reports the directional power transfer across an interface, UtilizationSamples reports the absolute value of flow relative to the interface's transfer capability (counting the effects of line outages). For example, a 100 MW symmetrically-constrained interface which is fully congested may have a flow of +100 or -100 MW, but in both cases the utilization will be 100%. If a 50 MW line in the interface went on outage, flow may drop to +50 or -50 MW, but utilization would remain at 100%.\n\nA UtilizationSamplesResult can be indexed by a Pair of region names and a timestamp to retrieve a vector of sample-level utilizations of the interface in that timestep. Given the absolute value nature of the outcome, results are independent of direction. Querying \"Region A\" => \"Region B\" will yield the same result as \"Region B\" => \"Region A\".\n\nExample:\n\nutils, =\n    assess(sys, SequentialMonteCarlo(samples=10), UtilizationSamples())\n\nsamples =\n    utils[\"Region A\" => \"Region B\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert samples isa Vector{Float64}\n@assert length(samples) == 10\n\nsamples2 =\n    utils[\"Region B\" => \"Region A\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert samples == samples2\n\nNote that this result specification requires large amounts of memory for larger sample sizes. See Utilization for sample-averaged utilization results when sample-level granularity isn't required.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.StorageEnergySamples","page":"Public API Reference","title":"PRASCore.Results.StorageEnergySamples","text":"StorageEnergySamples\n\nThe StorageEnergySamples result specification reports the sample-level state of charge of Storages, producing a StorageEnergySamplesResult.\n\nA StorageEnergySamplesResult can be indexed by storage device name and a timestamp to retrieve a vector of sample-level charge states for the device in the given timestep.\n\nExample:\n\nstorenergy, =\n    assess(sys, SequentialMonteCarlo(samples=10), StorageEnergySamples())\n\nsamples = storenergy[\"MyStorage123\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert samples isa Vector{Float64}\n@assert length(samples) == 10\n\nNote that this result specification requires large amounts of memory for larger sample sizes. See StorageEnergy for estimated average storage state of charge when sample-level granularity isn't required.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.GeneratorStorageEnergySamples","page":"Public API Reference","title":"PRASCore.Results.GeneratorStorageEnergySamples","text":"GeneratorStorageEnergySamples\n\nThe GeneratorStorageEnergySamples result specification reports the sample-level state of charge of GeneratorStorages, producing a GeneratorStorageEnergySamplesResult.\n\nA GeneratorStorageEnergySamplesResult can be indexed by generator-storage device name and a timestamp to retrieve a vector of sample-level charge states for the device in the given timestep.\n\nExample:\n\ngenstorenergy, =\n    assess(sys, SequentialMonteCarlo(samples=10), GeneratorStorageEnergySamples())\n\nsamples = genstorenergy[\"MyGeneratorStorage123\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert samples isa Vector{Float64}\n@assert length(samples) == 10\n\nNote that this result specification requires large amounts of memory for larger sample sizes. See GeneratorStorageEnergy for estimated average generator-storage state of charge when sample-level granularity isn't required.\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.GeneratorAvailability","page":"Public API Reference","title":"PRASCore.Results.GeneratorAvailability","text":"GeneratorAvailability\n\nThe GeneratorAvailability result specification reports the sample-level discrete availability of Generators, producing a GeneratorAvailabilityResult.\n\nA GeneratorAvailabilityResult can be indexed by generator name and timestamp to retrieve a vector of sample-level availability states for the unit in the given timestep. States are provided as a boolean with true indicating that the unit is available and false indicating that it's unavailable.\n\nExample:\n\ngenavail, =\n    assess(sys, SequentialMonteCarlo(samples=10), GeneratorAvailability())\n\nsamples = genavail[\"MyGenerator123\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert samples isa Vector{Bool}\n@assert length(samples) == 10\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.GeneratorStorageAvailability","page":"Public API Reference","title":"PRASCore.Results.GeneratorStorageAvailability","text":"GeneratorStorageAvailability\n\nThe GeneratorStorageAvailability result specification reports the sample-level discrete availability of GeneratorStorages, producing a GeneratorStorageAvailabilityResult.\n\nA GeneratorStorageAvailabilityResult can be indexed by generator-storage name and timestamp to retrieve a vector of sample-level availability states for the unit in the given timestep. States are provided as a boolean with true indicating that the unit is available and false indicating that it's unavailable.\n\nExample:\n\ngenstoravail, =\n    assess(sys, SequentialMonteCarlo(samples=10), GeneratorStorageAvailability())\n\nsamples = genstoravail[\"MyGenerator123\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert samples isa Vector{Bool}\n@assert length(samples) == 10\n\n\n\n\n\n","category":"type"},{"location":"api/public/#PRASCore.Results.LineAvailability","page":"Public API Reference","title":"PRASCore.Results.LineAvailability","text":"LineAvailability\n\nThe LineAvailability result specification reports the sample-level discrete availability of Lines, producing a LineAvailabilityResult.\n\nA LineAvailabilityResult can be indexed by line name and timestamp to retrieve a vector of sample-level availability states for the unit in the given timestep. States are provided as a boolean with true indicating that the unit is available and false indicating that it's unavailable.\n\nExample:\n\nlineavail, =\n    assess(sys, SequentialMonteCarlo(samples=10), LineAvailability())\n\nsamples = lineavail[\"MyLine123\", ZonedDateTime(2020, 1, 1, 0, tz\"UTC\")]\n\n@assert samples isa Vector{Bool}\n@assert length(samples) == 10\n\n\n\n\n\n","category":"type"}]
}
