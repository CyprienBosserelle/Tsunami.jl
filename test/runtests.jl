using TsunamiTools
using Test

# Okada (1985) Table 2 check. Okada's reference point is the fault's bottom
# corner with strike 90; okada85 expects the centroid, so shift accordingly.
function okada_ref(x,y,d,dip,L,W,rake,slip,u3)
    ue,un,uz,uze,uzn,unn,une,uen,uee = okada85(x-L/2,y-cosd(dip)*W/2,d-sind(dip)*W/2,90,dip,L,W,rake,slip,u3; nargout=9)
    return [ue,un,uz,uee,uen,une,unn,uze,uzn]
end

@testset "TsunamiTools" begin

    @testset "Okada 1985 reference values" begin
        @testset "Case 2" begin
            x=2; y=3; d=4; dip=70; L=3; W=2
            @test okada_ref(x,y,d,dip,L,W,0,1,0) ≈ [-8.689E-3,-4.298E-3,-2.747E-3,-1.220E-3,+2.470E-4,-8.191E-3,-5.814E-4,-5.175E-3,+2.945E-4] atol=1e-6  # strike
            @test okada_ref(x,y,d,dip,L,W,90,1,0) ≈ [-4.682E-3,-3.527E-2,-3.564E-2,-8.867E-3,-1.519E-4,+4.057E-3,-1.035E-2,+4.088E-3,+2.626E-3] atol=1e-5  # dip
            @test okada_ref(x,y,d,dip,L,W,0,0,1) ≈ [-2.660E-4,+1.056E-2,+3.214E-3,-5.655E-4,+1.993E-3,-1.066E-3,+1.230E-2,-3.730E-4,+1.040E-2] atol=1e-5  # tensile
        end
        @testset "Case 3" begin
            x=0; y=0; d=4; dip=90; L=3; W=2
            @test okada_ref(x,y,d,dip,L,W,0,1,0) ≈ [0,+5.253E-3,0,0,-1.864E-2,-2.325E-3,0,0,+2.289E-2] atol=1e-5
            @test okada_ref(x,y,d,dip,L,W,90,1,0) ≈ [0,0,0,0,+2.748E-2,0,0,0,-7.166E-2] atol=1e-5
            @test okada_ref(x,y,d,dip,L,W,0,0,1) ≈ [+1.223E-2,0,-1.606E-2,-4.182E-3,0,0,-2.325E-3,-9.146E-3,0] atol=1e-5
        end
        @testset "Case 4" begin
            x=0; y=0; d=6; dip=90; L=3; W=2
            @test okada_ref(x,y,d,dip,L,W,180,1,0) ≈ [0,-1.303E-3,0,0,+2.726E-3,+7.345E-4,0,0,-4.422E-3] atol=1e-5
            @test okada_ref(x,y,d,dip,L,W,90,1,0) ≈ [0,0,0,0,+5.157E-3,0,0,0,-1.901E-2] atol=1e-5
            @test okada_ref(x,y,d,dip,L,W,0,0,1) ≈ [+3.507E-3,0,-7.740E-3,-1.770E-3,0,0,-7.345E-4,-1.843E-3,0] atol=1e-5
        end
    end

    @testset "okada85 array interface" begin
        args = (20e3, 30.0, 15.0, 100e3, 50e3, 90.0, 5.0, 0.0)
        E = [-50e3 0.0 50e3; -40e3 10e3 60e3]
        N = [-30e3 -20e3 5e3; 40e3 70e3 -80e3]

        ue,un,uz = okada85(E,N,args...; nargout=3)
        @test size(ue) == size(un) == size(uz) == size(E)
        for I in eachindex(E)
            @test [ue[I],un[I],uz[I]] ≈ okada85(E[I],N[I],args...; nargout=3)
        end
        @test okada85(E,N,args...; nargout=1) ≈ uz
        @test length(okada85(E,N,args...; nargout=9)) == 9
        @test length(okada85(E,N,args...; nargout=12)) == 12

        # vectors work too
        ue_v,un_v,uz_v = okada85(vec(E),vec(N),args...)
        @test uz_v ≈ vec(uz)

        @test_throws ErrorException okada85(E,N,args...; nargout=2)
        @test_throws ErrorException okada85(E,N[:,1:2],args...)
    end

    @testset "Magnitude and slip" begin
        fault = faultparam(0,0,300e3,150e3,20e3,0,12,90,0,0,0)
        Calcslip!(fault, 9.1)
        @test fault.slip > 0
        @test CalcMw(fault) ≈ 9.1
        @test Mw2slip(8.0, 100e3, 50e3) ≈ fault.slip * 10^(1.5*(8.0-9.1)) * (300e3*150e3)/(100e3*50e3)
    end

    @testset "Unit conversion" begin
        fault = faultparam(0,0,300,150,10,0,12,90,0,0,0)
        faultkm2m!(fault)
        @test (fault.length, fault.width, fault.depth) == (300e3, 150e3, 10e3)
    end

    @testset "Rotations" begin
        x, y, xo, yo = 3.0, -2.0, 1.0, 4.0
        for theta in (0.0, 0.3, pi/2, -2.0)
            @test collect(unrotatexy(rotatexy(x,y,xo,yo,theta)...,xo,yo,theta)) ≈ [x,y]
        end
        for alpha in (0.0, 45.0, 192.0)
            @test collect(unrotatexyCompass(rotatexyCompass(x,y,xo,yo,alpha)...,xo,yo,alpha)) ≈ [x,y]
        end
        @test collect(rotatexy(1.0,0.0,0.0,0.0,pi/2)) ≈ [0.0,1.0] atol=1e-12
    end

    @testset "Spherical geometry" begin
        R = 6356750.52
        @test sphericDist(0.0,0.0,0.0,1.0) ≈ R*deg2rad(1.0)
        @test sphericDist(0.0,0.0,1.0,0.0) ≈ R*deg2rad(1.0)
        @test sphericDist(170.0,-40.0,170.0,-40.0) ≈ 0.0 atol=1e-6
        # crossing the antimeridian
        @test sphericDist(179.5,0.0,-179.5,0.0) ≈ R*deg2rad(1.0)

        lon, lat = sphericOffset(10.0, 0.0, 0.0, 1000.0)
        @test lon ≈ 10.0
        @test lat ≈ rad2deg(1000.0/R)
        lon, lat = sphericOffset(10.0, 0.0, 1000.0, 0.0)
        @test lon ≈ 10.0 + rad2deg(1000.0/R)
        @test lat ≈ 0.0 atol=1e-12
        lon, lat = sphericOffset(144.0, 38.0, 3000.0, -4000.0)
        @test sphericDist(144.0, 38.0, lon, lat) ≈ 5000.0 rtol=1e-6
    end

    @testset "Fault reference point" begin
        fault = faultparam(144.0,38.0,100e3,50e3,0.0,0.0,30.0,90.0,1.0,0,0)
        mvBLref2centroid!(fault)
        @test fault.depth ≈ sind(30)*50e3/2
        # strike 0: centroid is half a width east and half a length north
        @test fault.lat > 38.0
        @test fault.lon > 144.0
        @test sphericDist(144.0,38.0,fault.lon,fault.lat) ≈ hypot(25e3,50e3) rtol=1e-6
    end

    @testset "Grids" begin
        x = 1.0:1.0:3.0
        y = 10.0:1.0:11.0
        E, N = emptygrid(x, y)
        @test size(E) == size(N) == (3, 2)
        @test E[:,1] == collect(x) && E[:,2] == collect(x)
        @test N[1,:] == collect(y) && N[3,:] == collect(y)
        @test emptygrid(collect(x), collect(y)) == (E, N)
        @test emptygrid(1.0, 3.0, 10.0, 11.0, 1.0) == (E, N)
        @test emptygrid(Float32[1,2,3], Float32[10,11]) == (E, N)

        ef, nf = cartdistance2eq(E, N, 2.0, 10.5)
        @test ef == E .- 2.0
        @test nf == N .- 10.5

        ef, nf = cartsphdist2eq(E, N, 2.0, 10.5)
        @test sign.(ef) == sign.(E .- 2.0)
        @test sign.(nf) == sign.(N .- 10.5)
        @test abs(nf[1,2]) ≈ sphericDist(2.0, 11.0, 2.0, 10.5)
    end

    @testset "InitTsunami (projected)" begin
        x = collect(-200e3:10e3:200e3) .+ 5e3   # avoid sampling exactly on fault edges
        y = collect(-200e3:10e3:200e3) .+ 5e3
        fault = faultparam(0.0,0.0,100e3,50e3,20e3,0.0,15.0,90.0,5.0,0,0)

        flat = InitTsunami(x, y, zeros(length(x), length(y)), fault)
        E, N = emptygrid(x, y)
        _, _, uz = okada85(E, N, fault.depth, fault.strike, fault.dip, fault.length, fault.width, fault.rake, fault.slip, 0.0)
        @test size(flat) == (length(x), length(y))
        @test all(isfinite, flat)
        @test flat ≈ uz              # no bathymetry gradient -> only vertical displacement
        @test maximum(flat) > 0      # thrust fault produces uplift

        H = [4000.0 + 0.01*xi for xi in x, yj in y]   # sloping seafloor
        sloped = InitTsunami(x, y, H, fault)
        @test all(isfinite, sloped)
        @test !(sloped ≈ flat)
    end

    include("test_kajiura.jl")

    @testset "tsunamidemo / write2nc" begin
        mktempdir() do dir
            dz, uz = tsunamidemo(dir)
            f1 = joinpath(dir, "Vertical_displacement.nc")
            f2 = joinpath(dir, "Vertical_displacement_Flat_Bathymetry.nc")
            @test isfile(f1) && isfile(f2)
            @test size(dz) == size(uz) == (151, 101)
            @test maximum(uz) > 0
        end
    end
end
