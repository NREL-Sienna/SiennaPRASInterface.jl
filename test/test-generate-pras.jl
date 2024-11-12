import PowerSystems
import PowerSystemCaseBuilder
import CSV
import DataFrames
import Dates
import Test
import TimeSeries
using Dates: DateTime

const PSY = PowerSystems
const PSCB = PowerSystemCaseBuilder

"""
Use PSCB to build RTS-GMLC System and add outage data
as a Supplmental Attribute
"""
function get_rts_gmlc_outage()
    rts_da_sys = PSCB.build_system(PSCB.PSISystems, "RTS_GMLC_DA_sys")

    ###########################################
    # Parse the gen.csv and add OutageData
    # SupplementalAttribute to components for
    # which we have this data
    ###########################################
    gen_for_data = CSV.read(joinpath(@__DIR__, "descriptors/gen.csv"), DataFrames.DataFrame)

    for row in DataFrames.eachrow(gen_for_data)
        λ, μ = PRASInterface.rate_to_probability(row.FOR, row["MTTR Hr"])
        transition_data = PSY.GeometricDistributionForcedOutage(;
            mean_time_to_recovery=row["MTTR Hr"],
            outage_transition_probability=λ,
        )
        comp = PSY.get_component(PSY.Generator, rts_da_sys, row["GEN UID"])

        if ~(isnothing(comp))
            PSY.add_supplemental_attribute!(rts_da_sys, comp, transition_data)
            @info "Added outage data supplemental attribute to $(row["GEN UID"]) generator"
        else
            @warn "$(row["GEN UID"]) generator doesn't exist in the System."
        end
    end
    return rts_da_sys
end

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

@testset "RTS GMLC" begin
    rts_da_sys = get_rts_gmlc_outage()
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

    # Make a PRAS System from PSY-4.X System
    rts_pras_sys = generate_pras_system(rts_da_sys, PSY.Area)
    @test rts_pras_sys isa PRASInterface.PRAS.SystemModel

    @test test_names_equal(rts_pras_sys.regions.names, area_names)
    @test test_names_equal(rts_pras_sys.generators.names, generator_names)
    @test test_names_equal(rts_pras_sys.storages.names, storage_names)

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
    max_capacity =
        PSY.get_time_series_values(
            PSY.SingleTimeSeries,
            hydro_component,
            "max_active_power",
        ) .* PSY.get_max_active_power(hydro_component)
    @test isapprox(rts_pras_sys.generators.capacity[idx, 1], max_capacity[1]) broken = true
    @test all(isapprox.(rts_pras_sys.generators.capacity[idx, :], max_capacity)) broken =
        true

    # Test Assess Run
    sequential_monte_carlo = PRASInterface.PRAS.SequentialMonteCarlo(samples=2)
    shortfalls, = PRASInterface.PRAS.assess(
        rts_pras_sys,
        sequential_monte_carlo,
        PRASInterface.PRAS.Shortfall(),
    )
    lole = PRASInterface.PRAS.LOLE(shortfalls)
    eue = PRASInterface.PRAS.EUE(shortfalls)
    @test lole isa PRASInterface.PRAS.ReliabilityMetric
    @test eue isa PRASInterface.PRAS.ReliabilityMetric
    @test PRASInterface.PRAS.val(lole) >= 0 && PRASInterface.PRAS.val(lole) <= 10
    @test PRASInterface.PRAS.stderror(lole) >= 0 && PRASInterface.PRAS.stderror(lole) <= 10
    @test PRASInterface.PRAS.val(eue) >= 0 && PRASInterface.PRAS.val(eue) <= 10
    @test PRASInterface.PRAS.stderror(eue) >= 0 && PRASInterface.PRAS.stderror(eue) <= 10

    @testset "Lumped Renewable Generators" begin
        rts_pras_sys =
            generate_pras_system(rts_da_sys, PSY.Area, lump_region_renewable_gens=true)
        @test rts_pras_sys isa PRASInterface.PRAS.SystemModel
        @test test_names_equal(rts_pras_sys.regions.names, area_names)

        rts_pras_sys = generate_pras_system(
            rts_da_sys,
            PSY.Area,
            lump_region_renewable_gens=true,
            availability=false,
        )
        @test rts_pras_sys isa PRASInterface.PRAS.SystemModel
        @test test_names_equal(rts_pras_sys.regions.names, area_names)

        rts_pras_sys = generate_pras_system(
            rts_da_sys,
            PSY.Area,
            lump_region_renewable_gens=true,
            availability=false,
            export_location=joinpath(@__DIR__, "rts.pras"),
        )
        @test rts_pras_sys isa PRASInterface.PRAS.SystemModel
        @test test_names_equal(rts_pras_sys.regions.names, area_names)
        @test isfile(joinpath(@__DIR__, "rts.pras"))
        rts_pras_sys2 = PRASInterface.PRAS.SystemModel(joinpath(@__DIR__, "rts.pras"))
    end
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

@testset "RTS GMLC with default data" begin
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
    @test rts_pras_sys isa PRASInterface.PRAS.SystemModel
    @test test_names_equal(rts_pras_sys.regions.names, area_names)
    @test test_names_equal(rts_pras_sys.generators.names, generator_names)
    @test test_names_equal(rts_pras_sys.storages.names, storage_names)
    # 201_HYDRO_4 should have 5.15 FOR, 11.6 POR, and 22 MTTR
    idx = findfirst(x -> x == "201_HYDRO_4", rts_pras_sys.generators.names)
    λ, μ = PRASInterface.rate_to_probability(5.15, 22)
    @test array_all_equal(rts_pras_sys.generators.λ[idx, :], λ)
    @test array_all_equal(rts_pras_sys.generators.μ[idx, :], μ)
end