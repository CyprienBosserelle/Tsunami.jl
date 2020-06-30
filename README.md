# Tsunami.jl
Julia module with tools useful tsunami science.

## Demo for generating a tsunami initial wave
        
        using Tsunami.jl
        using NetCDF
        
        # Read bathy grid for nearfield simulations
        x=NetCDF.ncread("C:\\Users\\bosserellec\\Documents\\Work\\Tohoku_cut.nc","lon")
        y=NetCDF.ncread("C:\\Users\\bosserellec\\Documents\\Work\\Tohoku_cut.nc","lat")
        zb=NetCDF.ncread("C:\\Users\\bosserellec\\Documents\\Work\\Tohoku_cut.nc","z")
        
        # Convert topography to actual water dpth assuming 0.0 mean sea level
        zs=0.0;
        H=max.(0.0,zs.-zb)

        # Define fault parameters here relative to the bottom left reference point
        # Note: bottom left corner of fault with a strike of 192 means north east corner
        fault=faultparam(144.33,39.6,300.0,150.0,0.0,192.0,12.0,90.0,0.0,0.0,0.0);

        #convert km width lengtn and depth to m
        # (It is easier to think in km but safer to keep all in a standard unit [m])
        faultkm2m!(fault)

        #convert from bottom left reference point to centroid
        mvBLref2centroid!(fault)

        # Calculate the slip for a given earthquake magnitude
        Calcslip!(fault,9.1)

        # Calculate the initial tsunami wave
        dz=InitTsunamiGeo(x,y,H,fault)

        write2nc(x,y,dz,"C:\\Users\\bosserellec\\Documents\\Work\\Tohoku_fulldz.nc")

        uz=InitTsunamiGeo(x,y,zeros(size(H)),fault)
        write2nc(x,y,uz,"C:\\Users\\bosserellec\\Documents\\Work\\Tohoku_uz_only.nc")
