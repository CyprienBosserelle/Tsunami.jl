# Tsunami.jl
Julia module with tools useful for tsunami science. At this stage the module is only good for generating initial deformation from earthquake tsunami but recent work will be integrated for eruption tsunami and reusing online database.


## Demo for generating a tsunami initial wave
 ```julia      
        using Tsunami.jl
        using NetCDF
        
        # Read bathy grid for a nearfield simulations
        x=NetCDF.ncread("Tohoku_bathy.nc","lon")
        y=NetCDF.ncread("Tohoku_bathy.nc","lat")
        zb=NetCDF.ncread("Tohoku_bathy.nc","z")
        
        # Convert topography to actual water depth assuming 0.0 mean sea level
        zs=0.0;
        H=max.(zs,zs.-zb)

        # Define fault parameters here relative to the bottom left reference point
        # fault=faultparam(lon,lat,length,width,depth strike,dip,rake,slip,tinit,trise)
        # Note: bottom left corner of fault with a strike of 192 means north east corner
        # Also we keep slip to 0.0 but we will calculate it after
        fault=faultparam(144.33,39.6,300.0,150.0,0.0,192.0,12.0,90.0,0.0,0.0,0.0);

        # Convert width, length, and depth from km to m
        # (It is easier to think in km but safer to keep all in a standard unit [m])
        faultkm2m!(fault)

        # Convert from bottom left reference point to centroid (This is what our Okada function expects)
        # some seismologist like bottom-left reference and some like centroid
        mvBLref2centroid!(fault)

        # Calculate the slip for a given earthquake magnitude and store it in our fault parameter
        Calcslip!(fault,9.1)

        # Calculate the initial tsunami wave
        dz=InitTsunamiGeo(x,y,H,fault)

        write2nc(x,y,dz,"Tohoku_fulldz.nc")
  ```
        
 Tohoku_fulldz.nc file can now be used to initialise a tsunami wave in a hydrodynamics model 

        
