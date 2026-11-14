# Kajiura (1963) filter: smoothing of co-seismic seafloor deformation by the
# water column before it is used as an initial sea-surface elevation.
#
# The sea-surface response to an instantaneous seafloor uplift dz in water of
# constant depth h is dz filtered by 1/cosh(k h) in wavenumber space.
# Equivalently, in physical space, it is a convolution with Kajiura's Green's function
#
#     eta(x,y) = ∫∫ G(r/h) dz(x',y') dx' dy' / h^2
#
#     G(r) = 1/π Σ_{n≥0} (-1)^n (2n+1) / ((2n+1)^2 + r^2)^(3/2)
#
# where r is the horizontal distance between (x,y) and (x',y').
# The spatial form allows the local water depth to be used at each point,
# which handles variable bathymetry. See e.g. Kajiura (1963), Bull. Earthq.
# Res. Inst. 41, 535-571; Glimsdal et al. (2013), J. Geophys. Res. 118, 4993-5007.

"""
    kajiuraG(r; nterms=60, naverage=10)

Kajiura (1963) Green's function for the sea-surface response to an
instantaneous point uplift of the seafloor. `r` is the horizontal distance
divided by the water depth. The function is normalised so that
`∫∫ G(r) dA = 1` (with `dA` in units of depth²).

The alternating image series is summed with `nterms` terms and accelerated by
`naverage` rounds of averaging of partial sums (Euler transform). The defaults
are accurate to near machine precision for `r ≲ 8`. Beyond that, `G` is below
`1e-6` and loses relative precision.
"""
function kajiuraG(r::Real; nterms::Integer=60, naverage::Integer=10)
    nterms > naverage >= 0 || throw(ArgumentError("need nterms > naverage >= 0"))
    r2 = float(r)^2
    S = Vector{Float64}(undef, nterms)
    s = 0.0
    for n in 0:nterms-1
        a = 2n + 1
        s += (iseven(n) ? a : -a) / (a^2 + r2)^1.5
        S[n+1] = s
    end
    for m in 1:naverage
        for i in 1:nterms-m
            S[i] = 0.5 * (S[i] + S[i+1])
        end
    end
    return S[nterms-naverage] / pi
end

# Lookup table of kajiuraG on [0, rmax] with linear interpolation (the series is
# too expensive to evaluate for every pair of grid points).
struct KajiuraTable
    dr::Float64
    g::Vector{Float64}
end
KajiuraTable(rmax::Real; dr::Real=0.01) = KajiuraTable(dr, [kajiuraG(r) for r in 0:dr:rmax+2dr])

@inline function (t::KajiuraTable)(r::Real)
    u = r / t.dr
    i = floor(Int, u)
    f = u - i
    @inbounds return (1 - f) * t.g[i+1] + f * t.g[i+2]
end

# Width associated with each node of a (possibly non-uniform) monotonic grid
function _cellwidths(x::AbstractVector)
    n = length(x)
    n == 1 && return [1.0]
    w = similar(x, Float64)
    w[1] = abs(x[2] - x[1])
    w[n] = abs(x[n] - x[n-1])
    for i in 2:n-1
        w[i] = abs(x[i+1] - x[i-1]) / 2
    end
    return w
end

function _checkfilterinput(x, y, dz, H, rmax)
    size(dz) == (length(x), length(y)) || throw(DimensionMismatch("size(dz) = $(size(dz)) must be (length(x), length(y)) = $((length(x), length(y)))"))
    size(H) == size(dz) || throw(DimensionMismatch("size(H) = $(size(H)) must match size(dz) = $(size(dz))"))
    rmax > 0 || throw(ArgumentError("rmax must be positive"))
    (length(x) > 1 && length(y) > 1) || throw(ArgumentError("x and y need at least two points"))
    (issorted(x) || issorted(x, rev=true)) && (issorted(y) || issorted(y, rev=true)) || throw(ArgumentError("x and y must be monotonic"))
end

# Generic kernel: `offset(i,j,k,l)` returns the (east, north) distance in metres
# from target node (i,j) to source node (k,l); `area[k,l]` is the relative cell
# area; `di`, `dj` are the minimum node spacings in metres used to size the search window.
function _kajiura(dz, H, offset, area, di, dj, rmax, Hmin)
    nx, ny = size(dz)
    G = KajiuraTable(rmax)
    out = Matrix{Float64}(undef, nx, ny)
    Threads.@threads for j in 1:ny
        for i in 1:nx
            h = H[i, j]
            if !(h > Hmin)
                out[i, j] = dz[i, j]   # land / very shallow: no filtering
                continue
            end
            R = rmax * h
            wi = min(nx, ceil(Int, R / di))
            wj = min(ny, ceil(Int, R / dj))
            num = 0.0
            den = 0.0
            for l in max(1, j - wj):min(ny, j + wj), k in max(1, i - wi):min(nx, i + wi)
                oe, on = offset(i, j, k, l)
                r = hypot(oe, on) / h
                r <= rmax || continue
                w = G(r) * area[k, l]
                num += w * dz[k, l]
                den += w
            end
            out[i, j] = num / den
        end
    end
    return out
end

"""
    kajiurafilter(x, y, dz, H; rmax=6.0, Hmin=0.0)

Smooth the co-seismic seafloor deformation `dz` with the Kajiura (1963) filter,
for a grid in projected coordinates (`x`, `y` in metres, e.g. Easting and Northing).

At each wet node the deformation is convolved with Kajiura's Green's function
[`kajiuraG`](@ref) scaled by the local water depth `H[i,j]`. This is
equivalent to the `1/cosh(kH)` filter of linear wave theory, but also works
with variable depth. The kernel is truncated at `rmax` depths (about 99.96% of
its weight lies within 6 depths) and renormalised over the available
neighbours, which also handles the domain edges.

- `x`, `y`: monotonic vectors (or ranges) of node coordinates [m]
- `dz`: seafloor deformation, `size(dz) == (length(x), length(y))`
- `H`: water depth (positive in the ocean), same size as `dz` [m]
- Nodes with `H <= Hmin` (land) are returned unchanged.

The cost scales with `(rmax * H / dx)^2` per node, so very fine grids in deep
water are slow. Julia threads (`julia -t auto`) are used when available.

Returns a new array the size of `dz`.
"""
function kajiurafilter(x::AbstractVector, y::AbstractVector, dz::AbstractMatrix, H::AbstractMatrix; rmax::Real=6.0, Hmin::Real=0.0)
    _checkfilterinput(x, y, dz, H, rmax)
    ax = _cellwidths(x)
    ay = _cellwidths(y)
    area = ax * ay'
    di = minimum(abs, diff(x))
    dyv = minimum(abs, diff(y))
    offset(i, j, k, l) = (x[k] - x[i], y[l] - y[j])
    return _kajiura(dz, H, offset, area, di, dyv, float(rmax), float(Hmin))
end

"""
    kajiurafilterGeo(lon, lat, dz, H; rmax=6.0, Hmin=0.0)

Same as [`kajiurafilter`](@ref) for a geographic grid (`lon`, `lat` in decimal
degrees). Distances are computed in metres using a local flat-earth
approximation around each node, which is accurate at the scale of the kernel
(a few water depths).
"""
function kajiurafilterGeo(lon::AbstractVector, lat::AbstractVector, dz::AbstractMatrix, H::AbstractMatrix; rmax::Real=6.0, Hmin::Real=0.0)
    _checkfilterinput(lon, lat, dz, H, rmax)
    earthradius = 6356750.52   # same as sphericDist
    m = earthradius * pi / 180 # metres per degree of latitude
    ax = _cellwidths(lon)
    ay = _cellwidths(lat)
    area = [ax[k] * ay[l] * cosd(lat[l]) for k in eachindex(lon), l in eachindex(lat)]
    dlon = minimum(abs, diff(lon))
    dlat = minimum(abs, diff(lat))
    di = m * dlon * cosd(min(89.0, maximum(abs, lat)))  # conservative (smallest) E-W spacing
    offset(i, j, k, l) = (m * cosd(0.5 * (lat[j] + lat[l])) * (lon[k] - lon[i]), m * (lat[l] - lat[j]))
    return _kajiura(dz, H, offset, area, di, m * dlat, float(rmax), float(Hmin))
end
