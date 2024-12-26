function test_names_equal(
    x::Vector{String},
    y::Vector{String};
    left="Left-side",
    right="Right-side",
)
    x = sort(x)
    y = sort(y)
    idx = 1
    while idx <= length(x) && idx <= length(y)
        if x[idx] == y[idx]
            idx += 1
        elseif x[idx] != y[idx]
            prev_idx = max(idx - 2, 1)
            next_idx = min(idx + 2, length(x))
            @error "$left $(x[prev_idx:next_idx]) != $right $(y[prev_idx:next_idx])"
            return false
        end
    end
    if idx > length(x) && idx > length(y)
        return true
    end
    prev_idx = max(idx - 2, 1)
    next_idx = min(idx + 2, length(x))
    comparison = idx <= length(x) ? "more" : "less"
    @error "$left has $comparison elements $(x[prev_idx:next_idx]) than $right $(y[prev_idx:next_idx])"
    return false
end

function array_all_equal(x::AbstractVector{T}, v::T) where {T}
    for (i, xi) in enumerate(x)
        if !isapprox(xi, v)
            @error "array[$i]: $xi != $v"
            return false
        end
    end
    return true
end

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
        rts_pras_sys =
            generate_pras_system(rts_da_sys, PSY.Area, lump_region_renewable_gens=true)
        @test rts_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel
        @test test_names_equal(rts_pras_sys.regions.names, area_names)

        rts_pras_sys =
            generate_pras_system(rts_da_sys, PSY.Area, lump_region_renewable_gens=true)
        @test rts_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel
        @test test_names_equal(rts_pras_sys.regions.names, area_names)

        rts_pras_sys = generate_pras_system(
            rts_da_sys,
            PSY.Area,
            lump_region_renewable_gens=true,
            export_location=joinpath(@__DIR__, "rts.pras"),
        )
        @test rts_pras_sys isa SiennaPRASInterface.PRASCore.SystemModel
        @test test_names_equal(rts_pras_sys.regions.names, area_names)
        @test isfile(joinpath(@__DIR__, "rts.pras"))
        rts_pras_sys2 =
            SiennaPRASInterface.PRASCore.SystemModel(joinpath(@__DIR__, "rts.pras"))
    end
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

# TODO: We want to test time-series λ, μ
# TODO: test HybridSystems
# TODO: Unit test line_and_interfaces.jl
# TODO: Unit test util/sienna/helper_functions.jl
