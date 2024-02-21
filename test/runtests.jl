using Tsunami
using Test

@testset "Tsunami.jl" begin
    # Write your tests here.
    #fault=faultparam(110,0.0,300.0,150.0,0.0,192.0,12.0,90.0,0.0,0.0,0.0);
    Tsunami.testokada();

end
