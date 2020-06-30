
"""
 Tsunami module
 Collection of function useful for tsunami simulations.

     Cyprien Bosserelle 2020

"""
module Tsunami

    import Okada,NetCDF
    export InitTsunamiGeo,faultparam,faultkm2m!,rotatexy,unrotatexy,unrotatexyCompass,rotatexyCompass,sphericDist,sphericOffset,mvBLref2centroid!,emptygrid,cartsphdist2eq,cartdistance2eq,CalcMw,Mw2slip,Calcslip!,Okadavert

"""
Fault parameter structure to simplify tsunami generation from earthquake
"""
    mutable struct faultparam
        lon::Float64
        lat::Float64
        length::Float64
        width::Float64
        depth::Float64
        strike::Float64
        dip::Float64
        rake::Float64
        slip::Float64
        tinit::Float64
        trise::Float64
    end
"""
Generate tsunami initial wave for a Geographical domain (i.e. lat and lon coordinates)
    usage: dz=InitTsunamiGeo(xx,yy,H,fault)
    where xx and yy can be vector or a range
    H is an array of water depth [m]
    see tsunamidemo()

"""
    function InitTsunamiGeo(xx,yy,H,fault::faultparam)

        # Create the array of easting and northing where deformation will be calculated
        E,N=emptygrid(xx,yy)
        ef,nf=cartsphdist2eq(E,N,fault.lon,fault.lat);


        # Calculate dHdx and dHdy which we need to account for horizontal displacement comntribution
        dHdx=zeros(size(ef));
        dHdy=zeros(size(nf));

        dHdx[2:end,:]=diff(H,dims=1)./diff(ef,dims=1);
        dHdy[:,2:end]=diff(H,dims=2)./diff(nf,dims=2);

        # Calculate horizontal and vertical deformation
        uX,uY,uZ=Okada.Okada85Dis(ef,nf,fault.depth,fault.strike,fault.dip,fault.length,fault.width,fault.rake,fault.slip,0);

        dz = uZ .+ uX .* dHdx .+ uY .* dHdy

        return dz
    end
    function tsunamidemo(ncfile)

        #Read bathy grid for nearfield simulations
        x=NetCDF.ncread("C:\\Users\\bosserellec\\Documents\\Work\\Tohoku_cut.nc","lon")
        y=NetCDF.ncread("C:\\Users\\bosserellec\\Documents\\Work\\Tohoku_cut.nc","lat")
        zb=NetCDF.ncread("C:\\Users\\bosserellec\\Documents\\Work\\Tohoku_cut.nc","z")

        xx=x
        yy=y
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
    end


    function CalcMw(fault::faultparam)
        # Calculate Mw based on fault dimention and slip
        Mo=4.0e11*fault.slip*100*fault.length*100*fault.width*100;
        Mw=2/3*log10(Mo)-10.7;
        return Mw
    end

    function Mw2slip(Mw,length,width)
        #Calculate slip based on Mw and fault dimension in m
        Mo=10^((Mw+10.7).*(3.0/2.0))
        slip=Mo/(4.0e11*100*length*100*width*100);


        return slip
    end

    function Calcslip!(fault::faultparam,Mw)
        slip=Mw2slip(Mw,fault.length,fault.width)
        fault.slip=slip;
        return fault
    end


    function faultkm2m!(fault::faultparam)
        #modify fault length and width from km to m
        fault.length*=1000.0;
        fault.width*=1000.0;
        fault.depth*=1000.0;
        return fault
    end

    function rotatexy(x,y,xo,yo,theta)
        # Rotate x y coordinate by theta.
        # theta is in radian following the maths convention

        newx=cos(theta)*(x-xo)-sin(theta)*(y-yo);
        newy=sin(theta)*(x-xo)+cos(theta)*(y-yo);
        return newx,newy
    end

    function unrotatexy(x,y,xo,yo,theta)
        ## unRotate x y coordinate by theta.
        # theta is in radian following the maths convention
        newx=cos(-1*theta)*(x)-sin(-1*theta)*(y)+xo;
        newy=sin(-1*theta)*(x)+cos(-1*theta)*(y)+yo;

        return newx,newy
    end

    function unrotatexyCompass(x,y,xo,yo,alpha)
        ## Rotate x y coordinate by theta.
        # alpha is in degrees following the compass convention

        theta=-1.0*deg2rad(alpha);
        return unrotatexy(x,y,xo,yo,theta);
    end

    function rotatexyCompass(x,y,xo,yo,alpha)
        # Rotate x y coordinate by alpha.
        # alpha is in degrees following the compass convention
        theta=-1.0*deg2rad(alpha);
        return rotatexy(x,y,xo,yo,theta);

    end





    function sphericDist(lon1,lat1,lon2,lat2)
        # Calculate spherical distance between 2 pts
        # distance is in m

        earthradius=6356750.52;
        #  Determine proper longitudinal shift.
        delta=lon2-lon1;
        l=abs(delta);
        l=l>=180 ? 360-l : l;
        # Convert Decimal degrees to radians.
        beta1 = deg2rad(lat1);
        beta2 = deg2rad(lat2);
        l = deg2rad(l);
        # Calculate S/Bo subformulas.
        st = sqrt(((sin(l).*cos(beta2)).^2)+(((sin(beta2).*cos(beta1))-(sin(beta1).*cos(beta2).*cos(l))).^2));

        return asin(st)*earthradius
    end

    function sphericOffset(lon,lat,de,dn)
        # Offset lat and lon by given eats and north offsets meters
        # de: delta east;  dn: delta north  offsets are in meters



        earthradius=6356750.52;

        bearing=atan(dn,de);

        dist=hypot(dn,de);

        # this is a bit rubbish because it assumes a flat earth!
        # # Coordinate offsets in radians
        # dLat = dn/R;
        # dLon = de/(R*cos(deg2rad(lat)));
        #
        # # OffsetPosition, decimal degrees
        # latO = lat + rad2deg(dLat);
        # lonO = lon + rad2deg(dLon);
        beta1=deg2rad(lat);
        alpha1=deg2rad(lon);

        latO=asin(sin(beta1)*cos(dist/earthradius)+cos(beta1)*sin(dist/earthradius)*cos(bearing));
        lonO=lon+rad2deg(atan(sin(bearing)*sin(dist/earthradius)*cos(beta1),cos(dist/earthradius)-sin(beta1)*sin(latO)))

        return lonO,rad2deg(latO)
    end

    function moveRef2Centroid!(fault::faultparam)
        fault.depth=fault.depth+sind(fault.dip)*fault.width.*0.5;
        newX=0.5*fault.width
        newY=0.5*fault.length

        rotX,rotY=rotatexyCompass(newX,newY,0.0,0.0,fault.strike);
        fault.lon=fault.lon+rotX;
        fault.lat=fault.lat+rotY;
        return fault
    end


    function mvBLref2centroid!(fault::faultparam)
        # move reference coordinates and depth
        # from bottom left corner (relative to the strike: for stroke of zero it is the south west corner, for stike of 180 it would be the north east corner)
        # to centroid of the fault plane
        fault.depth=fault.depth+sind(fault.dip)*fault.width.*0.5;

        # Moving to the centroid is a simple rotation problem

        newX=0.5*fault.width
        newY=0.5*fault.length

        rotX,rotY=rotatexyCompass(newX,newY,0.0,0.0,fault.strike);

        # bearing=rad2deg(atan(rotY,rotX));

        # Warning this need to be in meters!!!
        newlon,newlat=sphericOffset(fault.lon,fault.lat,rotX,rotY);



        fault.lon=newlon;
        fault.lat=newlat;
        return fault


    end

    function emptygrid(xx::Vector{Float64},yy::Vector{Float64})
        # porduce arrays of lat,lon/easting,northing representing an empty grid



        E=xx*ones(size(yy))'
        N=collect(transpose(yy*ones(size(xx))'))

        return E,N
    end

    function emptygrid(x::StepRange,y::StepRange)
        # porduce arrays of lat,lon/easting,northing representing an empty grid



        xx=collect(xx);
        yy=collect(yy);



        return emptygrid(xx,yy)
    end

    function emptygrid(xo,xmax,yo,ymax,inc)
        # porduce arrays of lat,lon/easting,northing representing an empty grid

        x=xo:inc:xmax
        y=yo:inc:ymax



        return emptygrid(x,y)
    end



    function cartdistance2eq(e,n,eqx,eqy)
        de=e-eqx;
        dn=n-eqy;
    end

    function cartsphdist2eq(e,n,eqx,eqy)
        facE=ones(size(e));
        facE[e .< eqx] .= -1.0;
        facN=ones(size(n));
        facN[n .< eqy] .= -1.0;

        ef=sphericDist.(e,eqy,eqx,eqy).*facE;
        nf=sphericDist.(eqx,n,eqx,eqy).*facN;

        return ef,nf
    end
end

#write2nc writes a matrix to netcdf file
function write2nc(x,y,z,ncfile,varnames)

    xatts = Dict("longname" => "longitude",
      "units"    => "m")
    yatts = Dict("longname" => "latitude",
              "units"    => "m")
    varatts = Dict("longname" => "z",
              "units"    => "m")
    NetCDF.nccreate(ncfile,varnames[3],varnames[1],x,xatts,varnames[2],y,yatts,atts=varatts)
    NetCDF.ncwrite(x,ncfile,varnames[1]);
    NetCDF.ncwrite(y,ncfile,varnames[2]);
    NetCDF.ncwrite(z,ncfile,varnames[3]);
    NetCDF.ncclose();
end
function write2nc(x,y,z,ncfile)
    varnames=["x","y","z"];
    write2nc(x,y,z,ncfile,varnames)
end
