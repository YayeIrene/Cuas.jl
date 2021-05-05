using Cuas
using Test
using ExternalBallistics




@testset "Cuas.jl" begin
    # Write your tests here.
    drone = createTargetRect(15.0,25.0,800.0)
    ϵ = RectErrorB(25.0,100.0,3.0,-10.0)
    @test pHit(drone,ϵ)== 0.0881314035058921

end
