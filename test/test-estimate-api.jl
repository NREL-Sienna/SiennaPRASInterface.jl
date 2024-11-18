@testset "estimate_resource_adequacy" begin
    rts_da_sys = get_rts_gmlc_outage()

    shortfalls, = estimate_resource_adequacy(rts_da_sys; samples=2, seed=1)
    @test shortfalls isa PRASInterface.PRAS.ResourceAdequacy.ShortfallResult

    lole = PRASInterface.PRAS.LOLE(shortfalls)
    eue = PRASInterface.PRAS.EUE(shortfalls)
    @test lole isa PRASInterface.PRAS.ReliabilityMetric
    @test eue isa PRASInterface.PRAS.ReliabilityMetric
    @test PRASInterface.PRAS.val(lole) >= 0 && PRASInterface.PRAS.val(lole) <= 10
    @test PRASInterface.PRAS.stderror(lole) >= 0 && PRASInterface.PRAS.stderror(lole) <= 10
    @test PRASInterface.PRAS.val(eue) >= 0 && PRASInterface.PRAS.val(eue) <= 10
    @test PRASInterface.PRAS.stderror(eue) >= 0 && PRASInterface.PRAS.stderror(eue) <= 10
end