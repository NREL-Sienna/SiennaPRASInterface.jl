@testset "RTS GMLC DA" begin
    rts_da_sys = get_rts_gmlc_outage("DA")
    area_names = PSY.get_name.(PSY.get_components(PSY.Area, rts_da_sys))
    generator_names =
        PSY.get_name.(
            PSY.get_components(
                x -> PSY.get_available(x) && PSY.get_rating(x) > 0,
                Union{PSY.HydroDispatch, PSY.RenewableGen, PSY.ThermalGen},
                rts_da_sys,
            )
        )
    storage_names = PSY.get_name.(PSY.get_components(PSY.Storage, rts_da_sys))
    generatorstorage_names =
        PSY.get_name.(
            PSY.get_components(
                x -> PSY.get_available(x) && PSY.get_rating(x) > 0,
                Union{PSY.HydroEnergyReservoir, PSY.HybridSystem},
                rts_da_sys,
            )
        )
    line_names =
        PSY.get_name.(
            PSY.get_components(PSY.Branch, rts_da_sys) do c
                PSY.get_available(c) &&
                    !(any(c isa T for T in SiennaPRASInterface.TransformerTypes))  # From definitions.jl
                PSY.get_area(PSY.get_from_bus(c)) != PSY.get_area(PSY.get_to_bus(c))
            end
        )

    # Make a PRAS System from PSY-4.X System
    rts_pras_sys = generate_pras_system(rts_da_sys, PSY.Area)
    @test rts_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel

    @test test_names_equal(rts_pras_sys.regions.names, area_names)
    @test test_names_equal(rts_pras_sys.generators.names, generator_names)
    @test test_names_equal(rts_pras_sys.storages.names, storage_names)
    @test test_names_equal(rts_pras_sys.generatorstorages.names, generatorstorage_names)
    @test test_names_equal(rts_pras_sys.lines.names, line_names)

    # Test that timestamps look right
    # get time series length
    psy_ts = first(PSY.get_time_series_multiple(rts_da_sys, type=PSY.SingleTimeSeries))
    @test all(
        TimeSeries.timestamp(psy_ts.data) .== collect(DateTime.(rts_pras_sys.timestamps)),
    )

    # Test capacities
    # 201_HYDRO_4 should have max active power from time series
    idx = findfirst(x -> x == "201_HYDRO_4", rts_pras_sys.generators.names)
    hydro_component = PSY.get_component(PSY.HydroDispatch, rts_da_sys, "201_HYDRO_4")
    max_power_ts =
        floor.(
            PSY.get_time_series_values(
                PSY.SingleTimeSeries,
                hydro_component,
                "max_active_power",
            )
        )
    @test rts_pras_sys.generators.capacity[idx, 1] == max_power_ts[1]
    @test all(rts_pras_sys.generators.capacity[idx, :] .== max_power_ts)

    thermal_component = PSY.get_component(PSY.ThermalStandard, rts_da_sys, "322_CT_6")
    @test array_all_equal(
        rts_pras_sys.generators.capacity[
            findfirst(x -> x == "322_CT_6", rts_pras_sys.generators.names),
            :,
        ],
        Int(floor(PSY.get_max_active_power(thermal_component))),
    )

    load_values = zeros(Float64, length(rts_pras_sys.regions.names), length(psy_ts))
    for load in PSY.get_components(PSY.StaticLoad, rts_da_sys)
        region = PSY.get_area(PSY.get_bus(load))
        # Fast enough for # areas < 10
        idx = findfirst(x -> x == PSY.get_name(region), rts_pras_sys.regions.names)
        max_active_power =
            PSY.get_time_series_values(PSY.SingleTimeSeries, load, "max_active_power")
        load_values[idx, :] += max_active_power
    end
    @test all(rts_pras_sys.regions.load .== Int.(floor.(load_values)))

    @testset "Lumped Renewable Generators" begin
        rts_pras_sys = generate_pras_system(rts_da_sys, PSY.Area, true)
        @test rts_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel
        @test test_names_equal(rts_pras_sys.regions.names, area_names)

        rts_pras_sys = generate_pras_system(rts_da_sys, PSY.Area, true)
        @test rts_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel
        @test test_names_equal(rts_pras_sys.regions.names, area_names)

        rts_pras_sys =
            generate_pras_system(rts_da_sys, PSY.Area, true, joinpath(@__DIR__, "rts.pras"))
        @test rts_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel
        @test test_names_equal(rts_pras_sys.regions.names, area_names)
        @test isfile(joinpath(@__DIR__, "rts.pras"))
        rts_pras_sys2 =
            SiennaPRASInterface.PRASCore.SystemModel(joinpath(@__DIR__, "rts.pras"))
    end
end

@testset "RTS GMLC DA with RATemplate" begin
    rts_da_sys = get_rts_gmlc_outage("RT")
    area_names = PSY.get_name.(PSY.get_components(PSY.Area, rts_da_sys))
    load_names = PSY.get_name.(PSY.get_components(PSY.StaticLoad, rts_da_sys))
    generator_names =
        PSY.get_name.(
            PSY.get_components(
                x -> PSY.get_available(x) && PSY.get_rating(x) > 0,
                Union{PSY.HydroDispatch, PSY.RenewableGen, PSY.ThermalGen},
                rts_da_sys,
            )
        )
    storage_names = PSY.get_name.(PSY.get_components(PSY.Storage, rts_da_sys))
    generatorstorage_names =
        PSY.get_name.(
            PSY.get_components(
                x -> PSY.get_available(x) && PSY.get_rating(x) > 0,
                Union{PSY.HydroEnergyReservoir, PSY.HybridSystem},
                rts_da_sys,
            )
        )
    line_names =
        PSY.get_name.(
            PSY.get_components(PSY.Branch, rts_da_sys) do c
                PSY.get_available(c) &&
                    !(any(c isa T for T in SiennaPRASInterface.TransformerTypes))  # From definitions.jl
                PSY.get_area(PSY.get_from_bus(c)) != PSY.get_area(PSY.get_to_bus(c))
            end
        )
    for (type, names) in zip(
        [
            PSY.StaticLoad,  # Bajer has two ElectricLoads
            PSY.Generator,
            PSY.Storage,
            Union{PSY.HydroEnergyReservoir, PSY.HybridSystem},
        ],
        [load_names, generator_names, storage_names, generatorstorage_names],
    )
        for name in names
            comp = PSY.get_component(type, rts_da_sys, name)
            if PSY.has_time_series(PSY.SingleTimeSeries, comp, "max_active_power")
                sts = PSY.get_time_series(PSY.SingleTimeSeries, comp, "max_active_power")
                PSY.remove_time_series!(
                    rts_da_sys,
                    PSY.SingleTimeSeries,
                    comp,
                    "max_active_power",
                )
                sts.name = "max_active_POWER"
                PSY.add_time_series!(rts_da_sys, comp, sts)
            end
        end
    end

    problem_template = SiennaPRASInterface.RATemplate(
        PSY.Area,
        [
            SiennaPRASInterface.DeviceRAModel(PSY.Line, SiennaPRASInterface.LinePRAS()),
            SiennaPRASInterface.DeviceRAModel(
                PSY.MonitoredLine,
                SiennaPRASInterface.LinePRAS(),
            ),
            SiennaPRASInterface.DeviceRAModel(
                PSY.TwoTerminalHVDCLine,
                SiennaPRASInterface.LinePRAS(),
            ),
            SiennaPRASInterface.DeviceRAModel(
                PSY.StaticLoad,
                SiennaPRASInterface.StaticLoadPRAS(max_active_power="max_active_POWER"),
            ),
            SiennaPRASInterface.DeviceRAModel(
                PSY.ThermalGen,
                SiennaPRASInterface.GeneratorPRAS(max_active_power="max_active_POWER"),
            ),
            SiennaPRASInterface.DeviceRAModel(
                PSY.HydroDispatch,
                SiennaPRASInterface.GeneratorPRAS(max_active_power="max_active_POWER"),
            ),
            SiennaPRASInterface.DeviceRAModel(
                PSY.RenewableGen,
                SiennaPRASInterface.GeneratorPRAS(max_active_power="max_active_POWER"),
            ),
            SiennaPRASInterface.DeviceRAModel(
                PSY.EnergyReservoirStorage,
                SiennaPRASInterface.EnergyReservoirLossless(),
            ),
            SiennaPRASInterface.DeviceRAModel(
                PSY.HydroEnergyReservoir,
                SiennaPRASInterface.HydroEnergyReservoirPRAS(
                    max_active_power="max_active_POWER",
                ),
            ),
        ],
    )
    # Make a PRAS System from PSY-4.X System
    rts_pras_sys = generate_pras_system(rts_da_sys, problem_template)
    @test rts_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel

    @test test_names_equal(rts_pras_sys.regions.names, area_names)
    @test test_names_equal(rts_pras_sys.generators.names, generator_names)
    @test test_names_equal(rts_pras_sys.storages.names, storage_names)
    @test test_names_equal(rts_pras_sys.generatorstorages.names, generatorstorage_names)
    @test test_names_equal(rts_pras_sys.lines.names, line_names)

    # Test that timestamps look right
    # get time series length
    psy_ts = first(PSY.get_time_series_multiple(rts_da_sys, type=PSY.SingleTimeSeries))
    @test all(
        TimeSeries.timestamp(psy_ts.data) .== collect(DateTime.(rts_pras_sys.timestamps)),
    )

    load_values = zeros(Float64, length(rts_pras_sys.regions.names), length(psy_ts))
    for load_name in load_names
        load = PSY.get_component(PSY.StaticLoad, rts_da_sys, load_name)
        region = PSY.get_area(PSY.get_bus(load))
        # Fast enough for # areas < 10
        idx = findfirst(x -> x == PSY.get_name(region), rts_pras_sys.regions.names)
        max_active_power =
            PSY.get_time_series_values(PSY.SingleTimeSeries, load, "max_active_POWER")
        load_values[idx, :] += max_active_power
    end
    # The rounding here gets weird because there's a 1199.9999999999995 load value
    @test all(rts_pras_sys.regions.load .<= floor.(Int64, load_values .+ 0.01))
    @test all(rts_pras_sys.regions.load .>= floor.(Int64, load_values .- 0.01))
    # Test capacities
    # 201_HYDRO_4 should have max active power from time series
    idx = findfirst(x -> x == "201_HYDRO_4", rts_pras_sys.generators.names)
    hydro_component = PSY.get_component(PSY.HydroDispatch, rts_da_sys, "201_HYDRO_4")
    max_power_ts =
        floor.(
            PSY.get_time_series_values(
                PSY.SingleTimeSeries,
                hydro_component,
                "max_active_POWER",
            )
        )
    @test rts_pras_sys.generators.capacity[idx, 1] == max_power_ts[1]
    @test all(rts_pras_sys.generators.capacity[idx, :] .== max_power_ts)

    thermal_component = PSY.get_component(PSY.ThermalStandard, rts_da_sys, "322_CT_6")
    @test array_all_equal(
        rts_pras_sys.generators.capacity[
            findfirst(x -> x == "322_CT_6", rts_pras_sys.generators.names),
            :,
        ],
        Int(floor(PSY.get_max_active_power(thermal_component))),
    )
end

@testset "RTS GMLC DA with default data" begin
    rts_da_sys = PSCB.build_system(PSCB.PSISystems, "RTS_GMLC_DA_sys")
    area_names = PSY.get_name.(PSY.get_components(PSY.Area, rts_da_sys))
    generator_names =
        PSY.get_name.(
            PSY.get_components(
                x -> PSY.get_available(x) && PSY.get_rating(x) > 0,
                Union{PSY.HydroDispatch, PSY.RenewableGen, PSY.ThermalGen},
                rts_da_sys,
            )
        )
    storage_names = PSY.get_name.(PSY.get_components(PSY.Storage, rts_da_sys))

    rts_pras_sys = generate_pras_system(rts_da_sys, PSY.Area)
    @test rts_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel
    @test test_names_equal(rts_pras_sys.regions.names, area_names)
    @test test_names_equal(rts_pras_sys.generators.names, generator_names)
    @test test_names_equal(rts_pras_sys.storages.names, storage_names)
    # 201_HYDRO_4 should have 5.15 FOR, 11.6 POR, and 22 MTTR
    idx = findfirst(x -> x == "201_HYDRO_4", rts_pras_sys.generators.names)
    λ, μ = SiennaPRASInterface.rate_to_probability(5.15, 22)
    @test array_all_equal(rts_pras_sys.generators.λ[idx, :], λ)
    @test array_all_equal(rts_pras_sys.generators.μ[idx, :], μ)
end

@testset "RTS GMLC RT with default data" begin
    rts_rt_sys = PSCB.build_system(PSCB.PSISystems, "RTS_GMLC_RT_sys")

    rts_pras_sys = generate_pras_system(rts_rt_sys, PSY.Area)
    @test rts_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel

    # Test that timestamps look right
    # get time series length
    psy_ts = first(PSY.get_time_series_multiple(rts_rt_sys, type=PSY.SingleTimeSeries))
    @test all(
        TimeSeries.timestamp(psy_ts.data) .== collect(DateTime.(rts_pras_sys.timestamps)),
    )
end

@testset "Two Area PJM with default data AreaInterchange" begin
    pjm_sys = PSCB.build_system(PSCB.PSISystems, "two_area_pjm_DA")
    template = RATemplate(PSY.Area, copy(SiennaPRASInterface.DEFAULT_DEVICE_MODELS))
    set_device_model!(
        template,
        PSY.AreaInterchange,
        SiennaPRASInterface.AreaInterchangeLimit,
    )
    pjm_pras_sys = generate_pras_system(pjm_sys, template)
    @test pjm_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel

    # Test that PRAS Interfaces limit_forward and limit_backward look right
    area_interchange = only(PSY.get_components(PSY.AreaInterchange, pjm_sys))

    @test size(pjm_pras_sys.interfaces.limit_forward) == (1, 168)
    @test pjm_pras_sys.interfaces.regions_from == [1]
    @test pjm_pras_sys.interfaces.regions_to == [2]
    if pjm_pras_sys.regions.names[1] == PSY.get_name(PSY.get_from_area(area_interchange))
        @test pjm_pras_sys.regions.names == [
            PSY.get_name(PSY.get_from_area(area_interchange)),
            PSY.get_name(PSY.get_to_area(area_interchange)),
        ]
        @test all(
            pjm_pras_sys.interfaces.limit_forward .==
            Int(floor(PSY.get_flow_limits(area_interchange).from_to)),
        )
        @test all(
            pjm_pras_sys.interfaces.limit_backward .==
            Int(floor(PSY.get_flow_limits(area_interchange).to_from)),
        )
    else
        @test pjm_pras_sys.regions.names == [
            PSY.get_name(PSY.get_to_area(area_interchange)),
            PSY.get_name(PSY.get_from_area(area_interchange)),
        ]
        @test all(
            pjm_pras_sys.interfaces.limit_forward .==
            Int(floor(PSY.get_flow_limits(area_interchange).to_from)),
        )
        @test all(
            pjm_pras_sys.interfaces.limit_backward .==
            Int(floor(PSY.get_flow_limits(area_interchange).from_to)),
        )
    end

    line_names =
        PSY.get_name.(
            PSY.get_components(Union{PSY.MonitoredLine, PSY.Line}, pjm_sys) do c
                PSY.get_available(c) &&
                    PSY.get_area(PSY.get_from_bus(c)) != PSY.get_area(PSY.get_to_bus(c))
            end
        )
    @test test_names_equal(pjm_pras_sys.lines.names, line_names)
end

# TODO: We want to test time-series λ, μ
# TODO: test HybridSystems
# TODO: Unit test line_and_interfaces.jl
# TODO: Unit test util/sienna/helper_functions.jl
