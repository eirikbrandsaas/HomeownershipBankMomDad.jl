"""
    emptymedian(v::AbstractVector)

Find the median, but allows for empty collection. If `isempty(v)==true` returns NaN.
"""
function emptymedian(v::AbstractVector)
    if isempty(v)
        return NaN
    else
        median(v)
    end
end

"""
    emptycor(v1::AbstractVector)

Find Pearsons correlations, but allows for empty collection. If `isempty(v)==true` for either vector it returns an NaN
"""
function emptycor(v1::AbstractVector,v2::AbstractVector)
    if isempty(v1) || isempty(v2)
        return NaN
    else
        cor(v1,v2)
    end
end

function update_prices!(price,np::NumPar)
  if np.ns == 1
    np.p_grd[1] = price
    @assert np.p_grd[1] == price
    @assert length(np.p_grd) == 1
  elseif np.ns==3
    np.p_grd = np.p_grd*price/np.p_grd[2]
    @assert np.p_grd[2] == price
    @assert length(np.p_grd) == 3
    @assert np.p_grd[1] < np.p_grd[2]
    @assert np.p_grd[2] < np.p_grd[3]
  else
   throw("No valid price grid")
  end 

  return nothing
    
end

function update_hgrd!(ho::Float64,np::NumPar)
  @assert length(np.h_grd)==2
  @assert np.h_grd[1]==1.0
  np.h_grd[2] = ho

  return nothing
end


function bequest(dfin::DataFrame,np::NumPar)
  df = dfin[dfin.iak .== np.nyoung,:]
  df.beq = df.bp + np.h_grd[df.ihp].*np.p_grd[df.is].*df.ownp + df.tp

  return df[!,[:id, :beq, :gen, :iak]]
end

function smm_obj(dfsim,dfdata,moments)

    val = 0
    for mom in moments
        val += (dfdata[!,mom][1] - dfsim[!,mom][1])^2/(dfdata[!,mom][1])^2
    end

    return val

end

function find_moments(np::NumPar,mp::ModPar,pan,period::String)
    parampath = joinpath(@__DIR__,"../../data/moments/")

    dmom=DataFrame(CSV.File(parampath * "moments_race_age"  * ".csv")) # Import datamoms
    select!(dmom,[:race,:owneryoung,:medwealthyoung,:medwealthold,:transfrateyoung,:tp2wealthpyoung,:h2wyoung,:wealthatpurchaseyoung, :r2incyoung],All()) # Change order of columns for readability
    dmom=dmom[dmom.race.=="All",:]


    dmom.type = ["Data"]
    smom = DataFrame(type="Model")
    
    pan = simpanel_addvars(np,mp,pan) # Add variables to panel
    df = DataFrame(pan)
    try 
      df=leftjoin(df,combine(groupby(df[df.ownk.==1,:],:id), :iak => first => :first_own),on=:id) # Find age at first own
    catch # When nobody evers owns you reduce over an empty collection so manually set equal to missing
        df[!,:first_own] .= missing
    end
    dfbeq = bequest(df,np)
    df = leftjoin(df,dfbeq,on=[:id,:iak,:gen])
    df = df[df[:,:age] .<=53,:]
    
    tmp = leftjoin(df,combine(groupby(df,:id),:xp=>first),on=:id) # Find first parental wealth obs
    df.wealthyp = ( quantile(tmp.xp_first,2/3) .<= tmp.xp_first ) # Divide sample into three based on initial parental wealth
    df.middlep = (quantile(tmp.xp_first,1/3) .<= tmp.xp_first .< quantile(tmp.xp_first,2/3))
    df.poorp = (tmp.xp_first .< quantile(tmp.xp_first,1/3))

    # Manipulate panel
    df.agep = df.age .+ 30

    ylo = 25
    yhi = 44
    olo = 55
    ohi = 70

    dfy = df[(df.age .>= ylo) .& (df.age .< yhi),:]
    dfo = df[(df.agep .>= olo) .& (df.agep .< ohi),:]
    
    smom[!,:owneryoung] .= mean(dfy.ownk)
    smom[!,:ownerold] .= mean(dfy.ownp)
    smom[!,:r2incyoung] .= mean(dfy.rent2inc[(dfy.ownk .== 0)  ,:])
    try
        smom[!,:first_ownyoung] .= mean(skipmissing(df.first_own[(df.age .>= ylo + 1) .& (df.age .< yhi),:]))*(np.age_grd[2]-np.age_grd[1]) + 23
    catch # When nobody evers owns you reduce over an empty collection so manually set equal to missing
        smom[!,:first_ownyoung] .= missing
    end
    smom[!,:wealthyoung] .= mean(dfy.xk)
    smom[!,:medwealthyoung] .= median(dfy.xk)
    smom[!,:wealthold] .= mean(skipmissing(dfo.xp))
    smom[!,:medwealthold] .= median(skipmissing(dfo.xp))
    smom[!,:tp2wealthpyoung] .= mean(dfo[(dfo.tp .>0.5) .& (dfo.tp2xp.<1.0),:tp2xp])
    smom[!,:h2wyoung] .= mean(skipmissing(dfy[(dfy.ownk .== 1) .& (dfy.xk .> 0)  ,:h2w]))
    smom[!,:bng2w] .= sum(skipmissing((df.savingsp + df.tp)[(df.agep .== 83)  ,:]))/sum(skipmissing(df.xk + df.xp))
    smom[!,:transfrateyoung] .= mean(dfy.tp .> 0)
    smom[!,:transfsizeyoung] .= mean(dfy.tp)/smom.transfrateyoung[1]
    smom[!,:ownerall] .= mean([df.ownk ; df[np.age_grd[df.iak.+np.nyoung].<75,:].ownp])
    df_owners = df[ismissing.(df[!,:first_own]) .== 0,:]
    df_owners = df_owners[(df_owners.age .>= ylo) .& (df_owners.age .< yhi),:] # Same age limitations!
    smom[!,:wealthatpurchaseyoung] .= mean(df_owners[df_owners.iak .== df_owners.first_own,:].xk)
    smom[!,:mortgyoung] .= abs(mean(dfy.mortg)/mean(dfy.ownk))
    smom[!,:medwealthpgradyoung] .=  emptymedian(dfy[(dfy.ownk .== 1)  ,:xp])/ emptymedian(dfy[(dfy.ownk .== 0)  ,:xp]) # avoid that the code crashes if nobody owns or rents (and generate NaNs instead)
    smom[!,:owneryoung_prich] .= mean(dfy.ownk[(dfy.wealthyp .== 1) ,:])
    smom[!,:owneryoung_pmiddle] .= mean(dfy.ownk[(dfy.middlep .== 1) ,:])
    smom[!,:owneryoung_ppoor] .= mean(dfy.ownk[(dfy.poorp .== 1) ,:])
    smom[!,:transfrateyoung_prich] .= mean(dfy.tp[(dfy.wealthyp .== 1),:] .> 0)
    smom[!,:transfrateyoung_pmiddle] .= mean(dfy.tp[(dfy.middlep .== 1),:] .> 0)
    smom[!,:transfrateyoung_ppoor] .= mean(dfy.tp[(dfy.poorp .== 1),:] .> 0)

    df_owners = df_owners[!,[:iak, :id, :ownk, :first_own, :tp]]
    df_owners = df_lag(df_owners,:tp,:iak,:id,forward=true)
    df_owners = df_lag(df_owners,:tp,:iak,:id,forward=false)
    df_owners = df_owners[df_owners.first_own .== df_owners.iak,:] # Keep only age at which they become first owners
    df_owners.transferbuyersamnt = df_owners.tp  + df_owners.Ftp + df_owners.Ltp
    df_owners.transferbuyersLamnt = df_owners.Ltp
    df_owners.transferbuyersNamnt = df_owners.tp
    df_owners.transferbuyersFamnt = df_owners.Ftp
    try
        smom[!,:transferbuyersyoung] .=  mean(skipmissing(df_owners.transferbuyersamnt .> 0))
        smom[!,:transferbuyersLyoung] .=  mean(skipmissing(df_owners.transferbuyersLamnt .> 0))
        smom[!,:transferbuyersNyoung] .=  mean(skipmissing(df_owners.transferbuyersNamnt .> 0))
        smom[!,:transferbuyersFyoung] .=  mean(skipmissing(df_owners.transferbuyersFamnt .> 0))

        smom[!,:transferbuyersamntyoung] .=  mean(skipmissing(df_owners.transferbuyersamnt))/smom.transferbuyersyoung[1]
        smom[!,:transferbuyersLamntyoung] .=  mean(skipmissing(df_owners.transferbuyersLamnt))/smom.transferbuyersLyoung[1]
        smom[!,:transferbuyersNamntyoung] .=  mean(skipmissing(df_owners.transferbuyersNamnt))/smom.transferbuyersNyoung[1]
        smom[!,:transferbuyersFamntyoung] .=  mean(skipmissing(df_owners.transferbuyersFamnt))/smom.transferbuyersFyoung[1]
    catch
        smom[!,:transferbuyersyoung] .= NaN
        smom[!,:transferbuyersamntyoung] .=  NaN
    end

    
    dfy_own = dfy[dfy.ownk .== 1,:]
    df_first = dfy_own[(dfy_own.iak .== dfy_own.first_own),:] # dataset of households at purchase
    smom[!,:LTVatpurchaseyoung] .= mean( -df_first.mortg ./ (np.p_grd[df_first.is].*np.h_grd[df_first.ihk]) )

    moms = vcat(dmom,smom,cols=:union)

    return moms

end


## Functions to find CDF's
function findoutcomeCDF(pr::Float64,CDF)
    @assert pr >= 0.0 && pr <= 1.0
    @assert CDF[1] >= 0.0
    @assert CDF[end] == 1.0

    outcome = 0
    for (i,cdf) in enumerate(CDF)
        if pr <= cdf
            outcome = i
            break
        end
    end

    @assert outcome != 0
    return outcome
end

function quick_findoutcomePDF(prob,probs)

  indx  = 0
  _CDF = probs[firstindex(probs)] # The value of the CDF and PDF is the same at the first grid point
  for idx in eachindex(probs)
    indx = idx
    if prob <= _CDF || idx == length(probs) # Stop searhing once CDF is higher than the random number OR you're at the final point
        break
    end
      _CDF += probs[idx+1] # Moving the CDF one grid point forward (by a manual cumulative sum)
      @assert 0 <=  _CDF <= 1.0 + eps(Float64)*4 "$_CDF is outside [0,1]"
  
  end 

  return indx
end

function CDF_1d(PDFvec::Vector)
    # Check that PDFs are actually PDFs
    @assert abs(sum(PDFvec) - 1) < 0.001
    @assert minimum(PDFvec) >= 0.0
    @assert maximum(PDFvec) <= 1.0

    CDF = cumsum(PDFvec) # Running sum to get CDF

    CDF ./= CDF[end]
    return CDF
end

function CDF_2d(PDFarr::Matrix)
    CDF = similar(PDFarr)
    [CDF[i,:] = CDF_1d(PDFarr[i,:]) for i = 1:size(PDFarr)[1]]
    return CDF
end


function itp_gp_preallocate(np::NumPar,gp)

  itertypec = Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 5}
  itertyped = Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Int64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 5}
  
  @views gtp_itp    = [linear_interpolation((np.xk_grd,np.xp_grd),   gp.t[:,:,ivk,ihk,iop,is,np.nyoung,ihp],extrapolation_bc=Flat()) for ivk=1:np.nv,ihk=1:np.no,iop=1:np.no,is=1:np.ns,ihp=1:np.nh]::itertypec
  @views gbp_itp    = [linear_interpolation((np.xk_grd,np.xp_grd),   gp.b[:,:,ivk,ihk,iop,is,np.nyoung,ihp],extrapolation_bc=Flat()) for ivk=1:np.nv,ihk=1:np.no,iop=1:np.no,is=1:np.ns,ihp=1:np.nh]::itertypec
  @views gdiscp_itp = [linear_interpolation((np.xk_grd,np.xp_grd),gp.disc[:,:,ivk,ihk,iop,is,np.nyoung,ihp],extrapolation_bc=Flat()) for ivk=1:np.nv,ihk=1:np.no,iop=1:np.no,is=1:np.ns,ihp=1:np.nh]::itertyped

  return gtp_itp::itertypec,gbp_itp::itertypec,gdiscp_itp::itertyped
end



function itp_gp!(np::NumPar,iak::Int,gp,Vknxt,Vpnxt,
  _tp_itp::Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 5},
  _bp_itp::Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 5},
  _discp::Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Int64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 5},
  iter)
  if iak < np.nyoung
    _iak = iak + 1   # Parents policy move in the next period, so must shift forward by one age
  else 
    _iak::Int = np.nyoung # But np.na + 1 doesnt exit, but these are never used, so just shift them back one period
  end        
  for ivk=1:np.nv,ihk=1:np.no,iop=1:np.no,is=1:np.ns,ihp=1:np.nh
     @views _tp_itp[ivk,ihk,iop,is,ihp].itp.coefs .=  gp.t[:,:,ivk,ihk,iop,is,_iak,ihp]
     @views _bp_itp[ivk,ihk,iop,is,ihp].itp.coefs .= gp.b[:,:,ivk,ihk,iop,is,_iak,ihp]
     @views _discp[ivk,ihk,iop,is,ihp].itp.coefs .= gp.disc[:,:,ivk,ihk,iop,is,_iak,ihp]
  end   

  itertypev = Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 4}

  if iak < np.nyoung # Note that iter+1 refers to using the value function from the current iteration (1st index is initial guess, which is never used for age 1:np.nyoung-1 since that's a life-cycle problem!)
    @views Vnxt = [linear_interpolation((np.xk_grd,np.bp_grd),Vknxt[:,:,ihpn,ivk,iok,is,iak+1,iter+1],extrapolation_bc=Flat()) for ihpn=1:np.nh,ivk=1:np.nv,iok=1:np.no,is=1:np.ns]::itertypev
  else # In next period kid becomes parent, so use parents value function (at first age of parenthood). Note that this is the value function from previous iteration
    @views Vnxt = [linear_interpolation((np.xk_grd,np.xp_grd),Vpnxt[:,:,ivk,ihk,iop,is,1,iter],extrapolation_bc=Flat()) for ihk=1:np.no,ivk=1:np.nv,iop=1:np.no,is=1:np.ns]::itertypev
  end

  return Vnxt::Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 4}
end

function itp_gk_preallocate(np::NumPar,gk,Vpnxt)
  itertypec = Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 5}
  itertyped = Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Int64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 5}
  itertypev = Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 4}

  @views    gck_itp = [linear_interpolation((np.xk_grd,np.bp_grd),      gk.c[:,:,ihp,ivk,iok,is,np.nyoung,ihk],extrapolation_bc=Flat()) for ihp=1:np.nh,ivk=1:np.nv,iok=1:np.no,is=1:np.ns, ihk=1:np.nh]::itertypec
  @views gdisck_itp = [linear_interpolation((np.xk_grd,np.bp_grd),   gk.disc[:,:,ihp,ivk,iok,is,np.nyoung,ihk],extrapolation_bc=Flat()) for ihp=1:np.nh,ivk=1:np.nv,iok=1:np.no,is=1:np.ns,ihk=1:np.nh]::itertyped
  @views Vnxt = [linear_interpolation((np.xk_grd,np.xp_grd),Vpnxt[:,:,ivk,iok,iop,is,np.nyoung]      ,extrapolation_bc=Flat()) for ivk=1:np.nv,iok=1:np.no,iop=1:np.no,is=1:np.ns]::itertypev
  
  return gck_itp,gdisck_itp,Vnxt
end

function itp_gk!(np::NumPar,iak,gk,Vpnxt,
  _ck_itp::Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 5},
  _disck_itp::Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Int64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 5},
  _Vpnxt_itp::Array{Interpolations.Extrapolation{Float64, 2, Interpolations.GriddedInterpolation{Float64, 2, Matrix{Float64}, Gridded{Linear{Throw{OnGrid}}}, Tuple{Vector{Float64}, Vector{Float64}}}, Gridded{Linear{Throw{OnGrid}}}, Flat{Nothing}}, 4},
  iter)

  for ihp=1:np.nh,ivk=1:np.nv,iok=1:np.no,is=1:np.ns, ihk=1:np.nh
    @views _ck_itp[ihp,ivk,iok,is,ihk].itp.coefs .= gk.c[:,:,ihp,ivk,iok,is,iak,ihk]
    @views _disck_itp[ihp,ivk,iok,is,ihk].itp.coefs .= gk.disc[:,:,ihp,ivk,iok,is,iak,ihk]
  end

  if iak < np.nyoung
    iakn = iak + 1
    itern = iter + 1 # Note that iter+1 refers to using the value function from the current iteration (1st index is initial guess, which is never used for age 1:np.nyoung-1 since that's a life-cycle problem!)
  else
    iakn::Int = 1
    itern = iter # Use final value function from previous iteration (initial guess)
  end          
  for ivk=1:np.nv,iok=1:np.no,iop=1:np.no,is=1:np.ns
    @views _Vpnxt_itp[ivk,iok,iop,is].itp.coefs .= Vpnxt[:,:,ivk,iok,iop,is,iakn,itern] 
  end

  return nothing
end

##
function init_sol(np::NumPar)
    g = Pol(np)

    g.gk.b .= 1.0
    g.gk.c .= 1.0
    g.gp.b .= 1.0
    g.gp.c .= 1.0

    for ih = 2:np.nh
      g.gk.disc[:,:,:,:,:,:,:,ih] .= 0.0
      g.gp.disc[:,:,:,:,:,:,:,ih] .= 0.0
    end
      
    return g
end

##
function find_deterministicwage(age,βvec,retage,retinc)

    if age<retage # For ages below retirement age it's a fifth degre polynomial
        y = βvec[1]
        for i = 1:length(βvec)-1
            y += βvec[i+1]*(age.^i)/(10^i)
        end

    else  # For retired households its the income in retired period multiplied by the λ retirement income shfiter
        y = find_deterministicwage(retage-1,βvec,retage,retinc)*retinc
    end
    return y
end


function fillΠs(ns::Int64)
  if ns == 1
    Πs = fill(1.0,ns,ns)
  elseif ns == 2
    Πs = [0.96 0.04;
          0.5 0.5]
  elseif ns == 3
    Πs = [0.9 0.1 0;
          0.02 0.96 0.02;
          0.0 0.25 0.75]
  elseif ns == 4
    Πs = [0.90 0.10 0.00 0.00;
          0.02 0.94 0.02 0.02;
          0.02 0.00 0.96 0.02;
          0.00 0.25 0.00 0.75]
  else
    throw("An invalid ns was inputed")
  end

  return Πs
end

## Functions to work with the initial distribution
function find_inidist(Θnxk::Int64,xk_grd::Vector{Float64},nv::Int64)
    nxk = length(xk_grd)
    parampath = joinpath(@__DIR__,"../../data/model/calibration")
    Θcsv = CSV.File(joinpath(parampath,"theta_"*@sprintf("%i",nv) * ".csv"))
    Θlimscsv = collect(skipmissing(CSV.File(joinpath(parampath,"thetaxplims_"*@sprintf("%i",nv) * ".csv"),header=false).Column1))
    @assert Θnxk == maximum(Θcsv.xk)
    @assert Θnxk == 4 # See collect(Θcsv.wealth[1:4]), which only counts to four for some reason

    # Allocate matrices
    Θ = zeros(Θnxk,nv,Θnxk,nv)

    for ind = 1:Θnxk*nv*Θnxk*nv
      ixn,ivn,ixk,ivk=Tuple(CartesianIndices((Θnxk,nv,Θnxk,nv))[ind])
      Θ[ixn,ivn,ixk,ivk] = Θcsv.prob[ind]
    end

    Θlims_tmp = vcat(Θlimscsv, Inf64) # Wealth thresholds
    Θlims = zeros(Θnxk,nv)

    θx_grd_tmp = collect(Θcsv.wealth[1:4]) # Imports the vector
    θx_grd = fill(0.0,(length(θx_grd_tmp),nv))
    for iv = 1:nv
        θx_grd[:,iv] = θx_grd_tmp
        Θlims[:,iv] = Θlims_tmp
    end

    Θgrid = zeros(Θnxk,nv,nxk,nv)
    for ixk = 1:nxk, ivk = 1:nv
        Θgrid[:,:,ixk,ivk] = Θpdf(xk_grd[ixk],ivk,Θ,Θlims[:,ivk],Θnxk)
        Θgrid[:,:,ixk,ivk] = Θgrid[:,:,ixk,ivk]/sum(Θgrid[:,:,ixk,ivk] )
    end



    return Θgrid,θx_grd
end

function Θpdf(xk::Real,iv::Int,Θ,θlims,Nxp::Int)
    ind = 0
    for ix=1:Nxp
        if xk<= θlims[ix]
            ind = ix
            break
        end
    end
    out = Θ[:,:,ind,iv]
    return out
end

function gridpower(length::Int64,gridmin::Real,gridmax::Real,coeff::Real,power::Real)
  x = range(0,stop=1,length=length)
  xx  = x + coeff*x.^power
  xmax = maximum(xx)
  xmin = minimum(xx)
  grid = (gridmax - gridmin)/(xmax - xmin)*xx .+ gridmin
end

## Random stuff


function df_lag(dfin::DataFrame,lagvar::Symbol,tvar::Symbol,idvar::Symbol;forward = false)
  df_lag = copy(dfin[!,[tvar, idvar]])
  if forward == false
      df_lag[!,tvar] = df_lag[!,tvar] .+ 1
      df_lag[!,Symbol("L"*string(lagvar))] = copy(dfin[!,lagvar])
  else
      df_lag[!,tvar] = df_lag[!,tvar] .- 1
      df_lag[!,Symbol("F"*string(lagvar))] = copy(dfin[!,lagvar])
  end
  dfout = leftjoin(dfin,df_lag,on=[idvar,tvar])
  return dfout
end


# Due to interpolation, sometimes the borrowing choice can be lower than the borrowing constraint, so manually reset it by lowering consumption and setting borrowing to the constraint
function enforce_feasibleborrowing(_c,_b,ih,is,ia,mp,np)
  maxborr = BorrCons(np.h_grd[ih],isowner(ih,np),np.p_grd[is],ia,mp)
  if _b < maxborr
      _b = maxborr
      _c -= (_b - maxborr)
  end
  @assert _c > 0 && _b ≥ maxborr

  return _c,_b
end


"""
    function benchpar()

Returns the benchmark parameters `mp` and `np` that should be used in all experiments.
"""
function benchpar(nstate,nchoice;robust_ν = false,robust_s=false)

    if robust_ν == false && robust_s == false ## Use benchmark numpar except in robustness exercises
        np = NumPar(nxk=nstate,nxp=nstate,nbp=nstate,
            ncp_choice=nchoice,nck_choice=nchoice,ntp_choice=nchoice,
            nxphat = nstate, ntp_choice2 = nchoice, ncp_choice2 = nchoice)
    elseif robust_ν == true && robust_s == false
        np = NumPar(nxk=nstate,nxp=nstate,nbp=nstate,
            ncp_choice=nchoice,nck_choice=nchoice,ntp_choice=nchoice,
            nν = 3)
    elseif robust_ν == false && robust_s == true
        np = NumPar(nxk=nstate,nxp=nstate,nbp=nstate,
            ncp_choice=nchoice,nck_choice=nchoice,ntp_choice=nchoice,
            ns = 3)
    else
        throw("No valid option inserted")
    end

    par = BestPars(nstate,nchoice)
    mp = ModPar(np,
        η = par.η,
        χ =  par.χ,
        ho = par.ho,
        )

    @assert np.nh == 2

    return mp, np
end

function simprobs(rng,nyoung,Ndyn,Ngen)
  simprob = ( # Probabilities
        ν = rand(rng,nyoung+1,Ndyn,Ngen),
        v = rand(rng,nyoung+1,Ndyn,Ngen),
        discbth = rand(rng,nyoung+1,Ndyn,Ngen),
        Θ = rand(rng,Ndyn,Ngen),
        δp = rand(rng,nyoung+1,Ndyn,Ngen),
        δk = rand(rng,nyoung+1,Ndyn,Ngen),
        s = rand(rng,nyoung+1,Ndyn,Ngen)
    ) 
end
