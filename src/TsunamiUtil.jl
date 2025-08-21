

"""
Fault parameter structure to simplify tsunami generation from earthquake
    faultparam(lon,lat,length,width,depth strike,dip,rake,slip,tinit,trise)
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
        CalcMw(fault::faultparam;rigidity=4.0e11)

    Calculate earthquake moment magnitude based on fault dimension and slip amount.

    """
    function CalcMw(fault::faultparam;rigidity=4.0e11)
        # Calculate Mw based on fault dimention and slip
        Mo=rigidity*fault.slip*100*fault.length*100*fault.width*100;
        Mw=2/3*log10(Mo)-10.7;
        return Mw
    end

    """
        Mw2slip(Mw,length,width;rigidity=4.0e11)

    Calculate earthquake slip amount on a fault based on fault dimension and moment magnitude
    """
    function Mw2slip(Mw,length,width;rigidity=4.0e11)
        #Calculate slip based on Mw and fault dimension in m
        Mo=10^((Mw+10.7).*(3.0/2.0))
        slip=Mo/(rigidity*100*length*100*width*100);


        return slip
    end

    """
        Calcslip!(fault::faultparam,Mw)

    Inplace calculation of fault slip amount given a fault and earthquake magnitude
    """
    function Calcslip!(fault::faultparam,Mw)
        slip=Mw2slip(Mw,fault.length,fault.width)
        fault.slip=slip;
        return fault
    end

    """
        faultkm2m!(fault::faultparam)
    Inplace conversion of fault dimension from kilometers to meters. Only applies to fault length, width and depth (does not applies to slip) 
    """
    function faultkm2m!(fault::faultparam)
        #modify fault length and width from km to m
        fault.length*=1000.0;
        fault.width*=1000.0;
        fault.depth*=1000.0;
        return fault
    end

    """
        rotatexy(x,y,xo,yo,theta)
    coordinate rotation function used to change between fault origin to centroid 
    """
    function rotatexy(x,y,xo,yo,theta)
        # Rotate x y coordinate by theta.
        # theta is in radian following the maths convention

        newx=cos(theta)*(x-xo)-sin(theta)*(y-yo);
        newy=sin(theta)*(x-xo)+cos(theta)*(y-yo);
        return newx,newy
    end

    """
        unrotatexy(x,y,xo,yo,theta)
    coordinate unrotation function used to change between fault origin to centroid 
    """
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



    """
        sphericDist(lon1,lat1,lon2,lat2)
    Harversine equation to calculate sperical distances
    """
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

    function sphericOffset(lon,lat,de,dn;earthradius=6356750.52)
        # Offset lat and lon by given eats and north offsets meters
        # de: delta east;  dn: delta north  offsets are in meters



        

        bearing=atan(de,dn);

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

    function mvBotref2centroid!(fault::faultparam)
        # move reference coordinates and depth
        # from  middle of bottom edge
        # to centroid of the fault plane
        fault.depth=fault.depth-sind(fault.dip)*fault.width.*0.5;

        # Moving to the centroid is a simple rotation problem

        newX=0.5*fault.width
        newY=0.0*fault.length

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

    function emptygrid(x::AbstractRange,y::AbstractRange)
        # porduce arrays of lat,lon/easting,northing representing an empty grid



        xx=collect(x);
        yy=collect(y);



        return emptygrid(xx,yy)
    end

    function emptygrid(xo,xmax,yo,ymax,inc)
        # porduce arrays of lat,lon/easting,northing representing an empty grid

        x=xo:inc:xmax
        y=yo:inc:ymax



        return emptygrid(x,y)
    end



    function cartdistance2eq(e,n,eqx,eqy)
        ef=e.-eqx;
        nf=n.-eqy;
        return ef,nf
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
        uX,uY,uZ = okada85(ef,nf,fault.depth,fault.strike,fault.dip,fault.length,fault.width,fault.rake,fault.slip,0.0; nargout=3)

        dz = uZ .+ uX .* dHdx .+ uY .* dHdy

        return dz
    end

    function zerobnd(dz;bndval=0.0)
        dz[1,:] .= bndval;
        dz[end,:] .= bndval;
        dz[:,1] .= bndval;
        dz[:,end] .= bndval;

        return dz
    end
 
    """
    Generate tsunami initial wave for a projected domain (i.e. Easting and Northing coordinates)
        usage: dz=InitTsunami(xx,yy,H,fault)
        where xx and yy can be vector or a range
        H is an array of water depth [m]
        see tsunamidemo()

    """
        function InitTsunami(xx,yy,H,fault::faultparam)

            # Create the array of easting and northing where deformation will be calculated
            E,N=emptygrid(xx,yy)

            ef,nf=cartdistance2eq(E,N,fault.lon,fault.lat);


            # Calculate dHdx and dHdy which we need to account for horizontal displacement comntribution
            dHdx=zeros(size(ef));
            dHdy=zeros(size(nf));

            dHdx[2:end,:]=diff(H,dims=1)./diff(ef,dims=1);
            dHdy[:,2:end]=diff(H,dims=2)./diff(nf,dims=2);

            # Calculate horizontal and vertical deformation
            uX,uY,uZ = okada85(ef,nf,fault.depth,fault.strike,fault.dip,fault.length,fault.width,fault.rake,fault.slip,0; nargout=3)

            dz = uZ .+ uX .* dHdx .+ uY .* dHdy

            return dz
        end
    function tsunamidemo()
        # First read a bathymetry
        # Read bathy grid from netcdf
        x=collect(135.0:0.1:150.0)
        y=collect(30.0:0.1:40.0)
        # create a dummy bathymetry
        zb=fill(-4000.0,(length(x),length(y)))
        # Put an island just off the center
        zb[Int(ceil(length(x)/3)):Int(ceil(length(x)/2)),Int(ceil(length(y)/3)):Int(ceil(length(y)/2))].=10.0

        # Convert topography to actual water dpth assuming 0.0 mean sea level
        zs=0.0;
        H=max.(0.0,zs.-zb)

        # Define fault parameters here relative to the bottom left reference point
        # here we put the ref point at the center of the grid
        # Normally the sismology determines the location of the reference poit of the fault
        flt_lon=minimum(x)+0.5*(maximum(x)-minimum(x));
        flt_lat=minimum(y)+0.5*(maximum(y)-minimum(y));
        # Note: bottom left corner of fault with a strike of 192 means north east corner
        fault=faultparam(flt_lon,flt_lat,300.0,150.0,0.0,192.0,12.0,90.0,0.0,0.0,0.0);

        # Convert km width length and depth to m
        # (It is easier to think in km but safer to keep all in a standard unit [m])
        faultkm2m!(fault)

        # Convert from bottom left reference point to centroid
        mvBLref2centroid!(fault)

        # Calculate the slip for a given earthquake magnitude
        Calcslip!(fault,9.1)

        # Calculate the initial tsunami wave
        dz=InitTsunamiGeo(x,y,H,fault)

        # Make sure we have zeros on the boundary (our model will extrapolate the bnd values to the rest of the grid)
        dz=zerobnd(dz);

        write2nc(x,y,dz,"Vertical_displacement.nc")

        # For comparison Let's redo the analysis but ignoring the
        # contribution the effect of coseismic horizontal displacements of
        # ocean bottom on the sea surface
        uz=InitTsunamiGeo(x,y,zeros(size(H)),fault)
        write2nc(x,y,uz,"Vertical_displacement_Flat_Bathymetry.nc")
    end

    function gaussian2dv(x,y,a,xo,yo,c)
    	return a * exp(-1.0 * ((x - xo) * (x - xo) + (y - yo) * (y - yo)) / (2.0 * c * c));
    end



    function gaussian2d(x,y,a,xo,yo,c)
    	zs=zeros(length(x),length(y))
    	for (i, v) in enumerate(x)
    		zs[i,:]=gaussian2dv.(v,y,a,xo,yo,c)
    	end
    	return zs
    end

    """
    sftstr="3.904*ki5a+8.345*ki6z+1.1*ki7b+11.433*ki7a+4.184*ki7z+1.847*ki8b+6.248*ki8a+5.957*ki8z+5.709*ki9z"
    """
    function SIFTTsunami(siftstring::String;dx=0.1)
        ################
        Fltseg,bbox=readSIFTstring(siftstring);

        x=collect(bbox[1]:dx:bbox[2])
        y=collect(bbox[3]:dx:bbox[4])
        # create a dummy bathymetry
        zb=fill(-4000.0,(length(x),length(y)))
        zs=0.0;
        H=max.(0.0,zs.-zb)

        dz=zeros(size(H));

        for i=1:length(Fltseg)

            dzthisflt=InitTsunamiGeo(x,y,H,Fltseg[i]);

            dz=dz.+dzthisflt;
        end

        dz=zerobnd(dz);

        write2nc(x,y,dz,"Vertical_displacement.nc")





    end

    function move_characters_to_end(s::String, idx1::Int)
    # Convert string to a mutable array of characters
        chars = collect(s)

        # Perform the swap if indices are valid
        if 1 <= idx1 <= length(chars) 
            inchar=chars[idx1];
            deleteat!(chars, idx1);
            push!(chars,inchar)
        else
            error("Invalid indices provided.")
        end

        # Convert the character array back to a string
        return String(chars)
    end

    function remove_chars_at_index(s::String, index::Int, count::Int=1)
        chars = collect(s) # Convert string to a character array
        if 1 <= index <= length(chars) && index + count - 1 <= length(chars)
            deleteat!(chars, index:index+count-1) # Remove characters at the specified range
        else
            println("Invalid index or count for removal.")
        end
        return String(chars) # Convert back to a string
    end

    function readSIFTstring(sftstr;siftfile="C:\\Users\\bosserellec\\Documents\\GitHub\\Tsunami.jl\\data\\SIFT_source_metadata.txt")
        ## Read the SIFT database
        
        ouflt=Vector{faultparam}();
        sdata=readdlm(siftfile,',',skipstart=2);

        ;

        segname=move_characters_to_end.(remove_chars_at_index.(String.(sdata[:,1]),3,2),3);


         nl,lc=size( sdata)

         siftindex=Dict{String,Int}()
        for i=1:nl
            siftindex[segname[i]]=i
        end

        ################
        ## Read the source String

        xmin=360.0;
        xmax=0.0;
        ymin=90.0;
        ymax=-90.0;


        #sftstr="3.904*ki5a+8.345*ki6z+1.1*ki7b+11.433*ki7a+4.184*ki7z+1.847*ki8b+6.248*ki8a+5.957*ki8z+5.709*ki9z"
        sftss=split(sftstr,"+");
        for i=1:length(sftss)
            sc=split(sftss[i],"*");
            slip=parse(Float64, sc[1])
            unit=sc[2];
            Long, Lat, Strike, Dip, Depth, Length, Width,Rake =  Float64.(sdata[siftindex[unit],[3,4,6,7,8,9,10,11]]);

            fltseg=faultparam(Long, Lat,Length, Width, Depth,Strike, Dip,Rake,slip,0.0,0.0);
            # Convert length width depth to km
            faultkm2m!(fltseg)
            # move fault lat lon and depth to segment centroid since that is what our Okada expects\
            # SIFT database seem to refer to the bottom centre of the fault.
            #mvBotref2centroid!(fltseg)
            push!(ouflt,fltseg);

            # Calculate the ideal bounding box for this string

            lonmax,latmax=sphericOffset(Long, Lat,5*Length*1000.0, 5*Length*1000.0);
            lonmin,latmin=sphericOffset(Long, Lat,-5*Length*1000.0, -5*Length*1000.0);

            xmin=min(xmin,lonmin);
            xmax=max(xmax,lonmax);

            ymin=min(ymin,latmin);
            ymax=max(ymax,latmax);

           

        end


        return ouflt,(xmin,xmax,ymin,ymax)
    end




