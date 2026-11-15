# Tests for the Kajiura (1963) filter.
# Run on its own with:  julia --project -e 'using TsunamiTools, Test; include("test/test_kajiura.jl")'
@testset "Kajiura filter" begin
    @testset "Green's function" begin
        # reference values from 1/(2π) ∫ k J0(k r)/cosh(k) dk (high-precision quadrature)
        @test kajiuraG(0) ≈ 0.2915609040308188 rtol=1e-10
        @test kajiuraG(1) ≈ 0.09037715249460897 rtol=1e-10
        @test kajiuraG(3) ≈ 0.002530922141285943 rtol=1e-8
        @test kajiuraG(-2) == kajiuraG(2)
        # ∫∫ G dA over r ≤ 6 (rectangle rule; the r = 0 term is zero)
        dr = 1e-3
        r = 0:dr:6
        @test sum(2pi .* r .* kajiuraG.(r)) * dr ≈ 0.9995894036248473 rtol=1e-4
        @test_throws ArgumentError kajiuraG(1.0; nterms=5, naverage=5)
    end

    h = 1000.0

    @testset "Matches 1/cosh(kh) for constant depth" begin
        x = collect(0.0:200.0:40e3)    # dx = h/5
        y = collect(0.0:200.0:16e3)
        H = fill(h, length(x), length(y))
        ic, jc = 101, 41               # node (20 km, 8 km), > 6h from the edges
        for kh in (0.5, 1.0, 2.0)
            dz = [cos(kh/h*(xi - 20e3)) for xi in x, yj in y]
            @test kajiurafilter(x, y, dz, H)[ic, jc] ≈ 1/cosh(kh) atol=2e-3
        end
    end

    @testset "Smoothing a bump conserves volume" begin
        x = 0.0:250.0:60e3
        y = 0.0:250.0:60e3
        H = fill(h, length(x), length(y))
        dz = [exp(-((xi-30e3)^2 + (yj-30e3)^2)/(2*(3e3)^2)) for xi in x, yj in y]
        out = kajiurafilter(x, y, dz, H)
        @test size(out) == size(dz)
        @test sum(out) ≈ sum(dz) rtol=1e-6
        @test 0.8*maximum(dz) < maximum(out) < maximum(dz)
    end

    @testset "Land and very shallow water are left unchanged" begin
        x = 0.0:200.0:4e3
        y = 0.0:200.0:4e3
        dz = [sin(xi/700) * cos(yj/900) for xi in x, yj in y]
        H = fill(h, size(dz))
        H[1:5, :] .= -5.0
        out = kajiurafilter(x, y, dz, H)
        @test out[1:5, :] == dz[1:5, :]
        @test !(out[10:end, :] ≈ dz[10:end, :])
        # depth much smaller than grid spacing: kernel narrower than one cell
        @test kajiurafilter(x, y, dz, fill(1.0, size(dz))) ≈ dz
    end

    @testset "Geographic grid matches projected grid near the equator" begin
        m = 6356750.52 * pi / 180
        lon = 0.0:0.002:0.3
        lat = -0.1:0.002:0.1
        dz = [exp(-((lo-0.15)^2 + la^2)/(2*0.02^2)) for lo in lon, la in lat]
        H = fill(h, size(dz))
        @test kajiurafilterGeo(lon, lat, dz, H) ≈ kajiurafilter(m .* lon, m .* lat, dz, H) rtol=1e-4
    end

    @testset "Input checks" begin
        x = 0.0:100.0:1000.0
        y = 0.0:100.0:500.0
        dz = zeros(length(x), length(y))
        @test_throws DimensionMismatch kajiurafilter(x, y, dz', fill(h, size(dz)))
        @test_throws DimensionMismatch kajiurafilter(x, y, dz, fill(h, 3, 3))
        @test_throws ArgumentError kajiurafilter(x, y, dz, fill(h, size(dz)); rmax=0)
        @test_throws ArgumentError kajiurafilter([0.0, 2.0, 1.0], y, zeros(3, length(y)), fill(h, 3, length(y)))
    end

    @testset "InitTsunami / InitTsunamiGeo with kajiura=true" begin
        x = collect(-200e3:10e3:200e3) .+ 5e3
        y = collect(-200e3:10e3:200e3) .+ 5e3
        H = fill(4000.0, length(x), length(y))
        fault = faultparam(0.0,0.0,100e3,50e3,20e3,0.0,15.0,90.0,5.0,0,0)
        raw = InitTsunami(x, y, H, fault)
        smooth = InitTsunami(x, y, H, fault; kajiura=true)
        @test smooth ≈ kajiurafilter(x, y, raw, H)
        @test all(isfinite, smooth)
        @test maximum(smooth) < maximum(raw)

        lon = collect(140.0:0.1:142.0) .+ 0.05
        lat = collect(37.0:0.1:39.0) .+ 0.05
        Hg = fill(4000.0, length(lon), length(lat))
        gfault = faultparam(141.0,38.0,100e3,50e3,20e3,0.0,15.0,90.0,5.0,0,0)
        rawg = InitTsunamiGeo(lon, lat, Hg, gfault)
        smoothg = InitTsunamiGeo(lon, lat, Hg, gfault; kajiura=true)
        @test smoothg ≈ kajiurafilterGeo(lon, lat, rawg, Hg)
        @test maximum(smoothg) < maximum(rawg)
    end
end
