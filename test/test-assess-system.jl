@testset "test assess(::PSY.System, ::Area, ...)" begin
    rts_da_sys = get_rts_gmlc_outage()

    sequential_monte_carlo = PRAS.SequentialMonteCarlo(samples=2, seed=1)
    shortfalls, =
        PRAS.assess(rts_da_sys, PSY.Area, sequential_monte_carlo, PRAS.Shortfall())
    @test shortfalls isa PRAS.ResourceAdequacy.ShortfallResult

    lole = PRAS.LOLE(shortfalls)
    eue = PRAS.EUE(shortfalls)
    @test lole isa PRAS.ReliabilityMetric
    @test eue isa PRAS.ReliabilityMetric
    @test PRAS.val(lole) >= 0 && PRAS.val(lole) <= 10
    @test PRAS.stderror(lole) >= 0 && PRAS.stderror(lole) <= 10
    @test PRAS.val(eue) >= 0 && PRAS.val(eue) <= 10
    @test PRAS.stderror(eue) >= 0 && PRAS.stderror(eue) <= 10
end
