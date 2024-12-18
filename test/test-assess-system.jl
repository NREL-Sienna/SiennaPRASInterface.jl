@testset "test assess(::PSY.System, ::Area, ...)" begin
    rts_da_sys = get_rts_gmlc_outage()

    sequential_monte_carlo = PRASCore.SequentialMonteCarlo(samples=2, seed=1)
    shortfalls, =
        PRASCore.assess(rts_da_sys, PSY.Area, sequential_monte_carlo, PRASCore.Shortfall())
    @test shortfalls isa PRASCore.ResourceAdequacy.ShortfallResult

    lole = PRASCore.LOLE(shortfalls)
    eue = PRASCore.EUE(shortfalls)
    @test lole isa PRASCore.ReliabilityMetric
    @test eue isa PRASCore.ReliabilityMetric
    @test PRASCore.val(lole) >= 0 && PRASCore.val(lole) <= 10
    @test PRASCore.stderror(lole) >= 0 && PRASCore.stderror(lole) <= 10
    @test PRASCore.val(eue) >= 0 && PRASCore.val(eue) <= 10
    @test PRASCore.stderror(eue) >= 0 && PRASCore.stderror(eue) <= 10
end
