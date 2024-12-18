@testset "test assess(::PSY.System, ::Area, ...)" begin
    rts_da_sys = get_rts_gmlc_outage()

    sequential_monte_carlo = SiennaPRASInterface.SequentialMonteCarlo(samples=2, seed=1)
    shortfalls, = SiennaPRASInterface.assess(
        rts_da_sys,
        PSY.Area,
        sequential_monte_carlo,
        SiennaPRASInterface.Shortfall(),
    )
    @test shortfalls isa SiennaPRASInterface.PRASCore.Results.ShortfallResult

    lole = SiennaPRASInterface.LOLE(shortfalls)
    eue = SiennaPRASInterface.EUE(shortfalls)
    @test lole isa SiennaPRASInterface.PRASCore.ReliabilityMetric
    @test eue isa SiennaPRASInterface.PRASCore.ReliabilityMetric
    @test SiennaPRASInterface.val(lole) >= 0 && SiennaPRASInterface.val(lole) <= 10
    @test SiennaPRASInterface.stderror(lole) >= 0 &&
          SiennaPRASInterface.stderror(lole) <= 10
    @test SiennaPRASInterface.val(eue) >= 0 && SiennaPRASInterface.val(eue) <= 10
    @test SiennaPRASInterface.stderror(eue) >= 0 && SiennaPRASInterface.stderror(eue) <= 10
end