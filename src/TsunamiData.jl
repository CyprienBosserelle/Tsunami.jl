
using Wavelets

function interp1(xA, yA, newx)
        # Interpolate to time nx even if xA is not monotonic (although it has to be increasing)
        # This is a naive implementation and could do with improvements.
        # This function deals with extrapolation by padding the first/last known value
        #

        #selcasttimedep[indnonan], error[indnonan], ttt

        #xA=selcasttimedep[indnonan]
        #yA=error[indnonan]
        #newx=copy(ttt)

		xA=vec(xA);
	    yA=vec(yA);


        y = zeros(length(newx))
        for n = 1:length(newx)
            index = findfirst(xA .> newx[n]);
            if index==nothing
                if newx[n] <= xA[1]
                    index = 1;
                elseif newx[n] >= xA[end]
                    index = length(xA)+1;
                end
            end
            prev = yA[max(index[1] - 1, 1)];
            next = yA[min(index[1],length(xA))];

            time = newx[n] - xA[max(index[1] - 1, 1)];
            timenext = xA[max(min(index[1],length(xA)), 1)] - xA[max(index[1] - 1, 1)];

            if max(min(index[1],length(xA)), 1) == max(index[1] - 1, 1)

                y[n] = yA[max(min(index[1],length(xA)), 1)];
            else
                y[n] = prev + (time) / (timenext) * (next - prev);

            end
        end
        return y;
end

    """
Low pass filtering using bootstrapping of residual of higher order
x: signal to low pass
wl: wavelet type (db4 is default)
n: highest order of the wavelet decomposition
"""
function wlopass(x;wl=WT.db4,n=5)
	lx=length(x)
	l=Wavelets.nextpow(2, lx)
	la=max(Int64(floor(((l-lx)*0.5))),1);
	xdetrend=x.-mean(x)
	xpad=zeros(l);
	xpad[la:(la+lx-1)]=xdetrend;

	wt = wavelet(wl)
	
	ywt=wpt(xpad, wt,n);

	ywt[first(detailrange(xpad,n)):last(detailrange(xpad,1))].=0

	xiwt=iwpt(ywt,wt,n);

	xlo=xiwt[la:(la+lx-1)].+mean(x)

	return xlo
end


function readIOC(file)
	data=readdlm(file,'\t',skipstart=3)
	ddd=DateTime.(data[:,2],"yyyy-mm-dd HH:MM:SS")
	sl=Float64.(data[:,1])
	return ddd,sl
end

function readDARTNOAA(file)
	data=readdlm(file,' ',skipstart=2)
	ddd=reverse(DateTime.(data[:,1],data[:,2],data[:,3],data[:,4],data[:,5],data[:,6]))
	sl=reverse(Float64.(data[:,8]))
	return ddd,sl
end

function readDARTNZ(file)
	data=readdlm(file,',',skipstart=1)
	ddd=DateTime.(data[:,7],"yyyy-mm-ddTHH:MM:SSZ")

	ps=sortperm(ddd)
	sl=Float64.(data[:,8])
	return ddd[ps],sl[ps]
end

function DetideNoaaDART(d,wl;n=7)
	xa=Float64.(Dates.value.(Dates.Second.(d.-d[1])))
	newx=collect(xa[1]:15:xa[end])
	newdate=d[1].+Dates.Second.(Int64.(round.(newx)))


	#Despike
	stdwl=std(wl)

	selectdata=abs.(wl.-mean(wl)) .< 5*stdwl;
	
	
	wlint=interp1(xa[selectdata], wl[selectdata], newx)

	wloint=wlopass(wlint,n=n)

	return newdate,wloint,wlint.-wloint
end
	
function readBGFTS(file;var="zs_N1")
	wl=dropdims(ncread(file,var), dims=(1, 2))
	time=Dates.Second.(Int64.(ncread(file,"time"))).+DateTime.("2025-07-29 23:24:50","yyyy-mm-dd HH:MM:SS")
	return time,wl
end

function readBGFTSO(file)
	data=readdlm(file,'\t',skipstart=2)
	ddd=DateTime.("2025-07-29 23:24:50","yyyy-mm-dd HH:MM:SS").+Dates.Millisecond.(Int64.(round.(Float64.(data[:,1]).*1000.0)))

	
	sl=Float64.(data[:,2])
	return ddd,sl
end


