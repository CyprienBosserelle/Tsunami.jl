
"""
 Tsunami module
 Collection of function useful for tsunami simulations.

     Cyprien Bosserelle and Takuya Miyashita 2020


"""
module Tsunami

    include("Okada.jl")
    include("TsunamiUtil.jl")
    include("TsunamiData.jl")
    import NetCDF
    using DelimitedFiles,Wavelets
    export InitTsunamiGeo,InitTsunami,faultparam,tsunamidemo,faultkm2m!,rotatexy,unrotatexy,unrotatexyCompass,rotatexyCompass,sphericDist,sphericOffset,mvBLref2centroid!,emptygrid,cartsphdist2eq,cartdistance2eq,CalcMw,Mw2slip,Calcslip!

end
