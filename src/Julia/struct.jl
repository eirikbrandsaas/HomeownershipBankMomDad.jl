mutable struct NumPar
  # Convergence criterias etc:
  vtol :: Float64
  viters :: Int64
  Ndyn :: Int
  Ngen :: Int

  # Grid sizes
  nxk :: Int64
  nxp :: Int64
  nbp :: Int64
  nv :: Int64
  nν :: Int64
  na :: Int64
  nyoung :: Int64
  nho :: Int64
  nhr :: Int64
  no :: Int64
  nh :: Int64
  ns :: Int64
  neh :: Int64
  xmaxk :: Float64
  xmaxp :: Float64
  xmin :: Float64
  Θnxk :: Int64
  nck_choice :: Int64
  ncp_choice :: Int64
  ntp_choice :: Int64
  nδ :: Int64

  # Exogenous grids
  h_grd :: Array{Float64,1}
  p_grd :: Array{Float64,1}
  y_grd :: Array{Float64,1}
  v_grd :: Array{Float64,2}
  ν_grd :: Array{Float64,2}
  θx_grd ::Array{Float64,2}
  age_grd :: Array{Float64,1}
  δ_grd :: Array{Float64,1}

  # Endogenous grids (states)
  xk_grd :: Array{Float64,1}
  xp_grd :: Array{Float64,1}
  bp_grd :: Array{Float64,1}

  # Choice grids
  ck_grd :: Array{Float64,1}
  cp_grd :: Array{Float64,1}
  tp_grd :: Array{Float64,1}

  #second stage grids for parents problem
  nxphat :: Int
  ntp_choice2 :: Int
  ncp_choice2 ::Int
  xphat_grd :: Array{Float64,1}
  tp_grd2 :: Array{Float64,1}
  cp_grd2 :: Array{Float64,1}

  # Transition matrices
  Πv :: Array{Float64,3}
  Πν :: Array{Float64,2}
  Πeh :: Array{Float64,1}
  Πs :: Array{Float64,2}
  ΠΘ ::Array{Float64,4}
  Πδ :: Array{Float64,1}

  # Transition CDFs
  ΠvCDF :: Array{Float64,3}
  ΠνCDF :: Array{Float64,2}
  ΠsCDF :: Array{Float64,2}
  ΠΘCDF ::Array{Float64,4}
  ΠδCDF :: Array{Float64,1}

  # Switches
  switch_str :: Dict{String,Any} # Not performance sensistive
  FinalTransfer :: Bool
  ParentsHousingChoice :: Bool
 
  # Other/random
  parvec :: DataFrame

  function NumPar(;vtol = 0.0000001, viters=35,Ndyn = 40000, Ngen = 5,nxk=10,nxp=10,nbp=10,nv=3,nν=1,na=30,nho=1,nhr=1,ns=1,neh=3,
                    xmaxk=1150.0,xmaxp=1150.0,xmin=-5.0,Θnxk=4,nck_choice=10,ncp_choice=10,ntp_choice=10,nδ=1,
                    nxphat = 20, ntp_choice2 = 60, ncp_choice2 = 60,
                    grdpwr = 8, grdcff = 2.0,period="benchmark",
                    parhousingchoice=true)

    switch_str = Dict("period" => period,
                  "IniDistr" => "constant",
                  )
          
    FinalTransfer = false # Turn of transfers in final period
    ParentsHousingChoice = parhousingchoice # Turn of parents housing choice

    parampath = joinpath(@__DIR__,"../../data/model/calibration")
    pvec = DataFrame(CSV.File(joinpath(parampath,"parameters.csv")))

    nh = nhr + nho
    no = nh
    nyoung = na÷2

    if nh == 2
      h_grd = [1.0, NaN64]
    end

    if ns == 3
      p_grd = [0.7, 1.0, 1.3]
    elseif ns == 1
      p_grd = [1.0]
    else 
      throw("No valid price grid size")
    end

    ## Income
    ageinc_betas = collect(CSV.File(joinpath(parampath,"ageincome_betas.csv"),header=false).Column1) #Import CSV, collect to get into a nice vector
    λ = pvec.lambda[1]
    retage = pvec.retage[1]

    agestep = round(Int64,(85-25)/na)
    age_grd = collect(25:agestep:83)
    y_grd = [find_deterministicwage(age_grd[ia],ageinc_betas,retage,λ) for ia = 1:na]

    nvstring = @sprintf("%i",nv)
    nνstring = @sprintf("%i",nν)

    Πv = fill(NaN64,(nyoung,nv,nv))
    v_grd = fill(NaN64,(nyoung,nv))

    ΠvDF = DataFrame(CSV.File(joinpath(parampath,"Piv_k_"*nvstring*".csv")))
    vDF = DataFrame(CSV.File(joinpath(parampath,"vgrid_k_"*nvstring*".csv")))
    for iak = 1:nyoung
      for iv = 1:nv
        Πv[iak,iv,:] = collect(ΠvDF[ΠvDF.agecl .== age_grd[iak],:][1,2+nv*(iv-1):1+nv+nv*(iv-1)])
        v_grd[iak,:]  = collect(vDF[vDF.agecl .== age_grd[iak],:][1,:2:nv+1])
      end
    end

    ν_grd = fill(NaN64,(na,nν))
    if nν > 1
      νDF = DataFrame(CSV.File(joinpath(parampath, "vpgrid_"*nνstring*".csv")))
      ν_grd[nyoung+1:na,:] = Matrix(νDF[:,2:end])

      iret = round(Int,(retage - 25)/2)
      Πν = fill(NaN64,na,nν)
      for ia = nyoung+1:iret
        Πν[ia,:]          = fill(1.0/nν,(nν))
      end
      for ia = iret+1:na
        Πν[ia,:]          = fill(1.0/nν,(nν))
      end

    else
      νDF = DataFrame(CSV.File(joinpath(parampath,"vpgrid_3.csv")))
      νDF = νDF[:,[1,2]]
      νDF[:,2] .= 1.0
      ν_grd[nyoung+1:na,:] = Matrix(νDF[:,2:end])
      Πν = fill(1.0,na,nν)
    end

    xk_grd = gridpower(nxk,xmin,xmaxk,grdcff,grdpwr)
    xp_grd = gridpower(nxp,xmin,xmaxp,grdcff,grdpwr)
    bp_grd = gridpower(nbp,-200,xmaxp,1,grdcff)*NaN # This grid is changed later, as the bottom limmit should depend on the maximum amount parents can borrow

    ck_grd = gridpower(nck_choice,1e-3,xmaxk,grdcff*2,grdpwr*1.5)
    cp_grd = gridpower(ncp_choice,1e-3,xmaxp,grdcff*2,grdpwr*1.5)
    tp_grd = gridpower(ntp_choice,0.0,xmaxp,grdcff,grdpwr) # Less dense at bottom grid for tp since whether people give 2 or 3 or 5 is not important, but you need density at the top for the richest households.

    # Grids only used in 2-stage parents problem
    xphat_grd = gridpower(nxphat,-35.0,xmaxp,grdcff*2,grdpwr*2)
    tp_grd2 = gridpower(ntp_choice2,0.0,xmaxp,grdcff,grdpwr)
    cp_grd2 = gridpower(ncp_choice2,0.01,xmaxp,grdcff*2,grdpwr*1.5)

    δ_grd = collect(range(0.0,stop=0.0,length=nδ))
    ## Initial distribution
    if Θnxk == 1
      Θnxk_bench = 4
      ΠΘtemp, θx_grdtemp = find_inidist(Θnxk_bench,xk_grd,nv)
      θx_grd = reshape(θx_grdtemp[2,:] .+ 27,(1,nv))
      ΠΘ = sum(ΠΘtemp;dims=1)
    else
      ΠΘ, θx_grd = find_inidist(Θnxk,xk_grd,nv)
    end

    Πeh = fill(0.0,neh)
    Πs = fill(0.0,ns,ns)
    # ΠΘ = fill(0.0,Θnxk,nv,nxk,nv)
    Πeh = fill(1.0,neh)*(1.0/neh)
    Πs = fillΠs(ns)
    # [ΠΘ[1,1,ixk,ivk] = 1.0 for ixk=1:nxk,ivk=1:nv]


    ΠvCDF = fill(NaN64,size(Πv))
    for ia = 1:nyoung
      ΠvCDF[ia,:,:] = CDF_2d(Πv[ia,:,:])
    end
    ΠsCDF = CDF_2d(Πs)
    if ns > 1
      @warn "In simulations, prices are constant at the medium level, though they expect prices can change. \n This is done by keeping the real PDF as-is but modifying the CDF, since the PDF is only used to solve the decision problem while the CDF only in simulation"
      ΠsCDF .= 0 
      ΠsCDF[:,div(ns,2,RoundUp):end] .= 1.0 # implies that price always goes to middle price level
    end

    ΠνCDF = fill(NaN64,size(Πν))
    for ia = nyoung+1:na
      ΠνCDF[ia,:] =  CDF_1d(Πν[ia,:])
    end
    ΠΘCDF = fill(0.0,size(ΠΘ))
    for ixk=1:nxk,ivk=1:nv
      ΠΘCDF[:,:,ixk,ivk] = reshape(CDF_1d(vec(ΠΘ[:,:,ixk,ivk])),(Θnxk,nv)) # Reshape and vectorize to make the 2-dimension into 1-dimension
    end

    Πδ = fill(1/nδ,nδ)
    ΠδCDF = CDF_1d(Πδ)

    new(vtol, viters,Ndyn, Ngen, nxk, nxp, nbp, nv, nν, na, nyoung,nho, nhr, no, nh, ns, neh, xmaxk, xmaxp, xmin, Θnxk,nck_choice, ncp_choice, ntp_choice, nδ,
          h_grd,p_grd,y_grd,v_grd,ν_grd,θx_grd,age_grd,δ_grd,xk_grd,xp_grd,bp_grd,ck_grd,cp_grd,tp_grd,
        nxphat, ntp_choice2, ncp_choice2, xphat_grd, tp_grd2, cp_grd2, 
        Πv,Πν,Πeh,Πs,ΠΘ,Πδ,ΠvCDF,ΠνCDF,ΠsCDF,ΠΘCDF,ΠδCDF,switch_str,FinalTransfer,ParentsHousingChoice,pvec)
  end
end


mutable struct ModPar
  β :: Float64
  γ :: Float64
  η :: Float64
  κ :: Float64
  χ :: Float64
  rf :: Float64
  rm :: Array{Float64,1}
  τb :: Float64
  τt :: Float64
  ξ :: Float64
  δ :: Float64
  LTV :: Array{Float64,1}
  PMIlim :: Float64
  PMIprem :: Float64
  ho :: Float64 # Size of owner-occupied home
  price :: Float64 # "benchmark" price

  mb :: Array{Float64,1}
  ms :: Array{Float64,1}
  rmp :: Float64
  rmk :: Float64

  # Various "loose ends"
  ihini :: Int64

  function ModPar(np::NumPar;β=0.92, γ=2.0, η=0.3, κ = 0.10, χ = 1.0 ,rf=0.04, τb=0.0, τt=0.0, ξ=0.2, δ=0.05, LTV=0.9, 
      PMIlim=0.78, PMIprem=0.03,ho = 2.5, mb = 0.02, ms = 0.075,rmk=0.03, rmp=0.0)
    ihini = 1
    mb = fill(mb,np.na)
    ms = fill(ms,np.na)
    if np.ParentsHousingChoice == false
      mb[np.nyoung+1:end-1] .= 1e30
      ms[np.nyoung+1:end-1] .= 1e30
    end

    LTV = fill(LTV,np.na)
    rm = [fill(rmk,np.nyoung);fill(rmk + rmp,np.nyoung)]
    
    update_hgrd!(ho,np)
    price = np.parvec.valhouse[1]/ho # Rescale so that price*ho = calibrated value
    update_prices!(price,np)
    
    # Set the lower limit of the parent's bond choice to be at the maximum LTV limit
    MaxBorr = maximum(LTV)*maximum(np.h_grd)*maximum(np.p_grd)
    np.bp_grd[:] = gridpower(np.nbp,-MaxBorr,np.xmaxp,1,3) # Somewhat unevenly spaced as others, since there is less curvature at the bottom, less mass of parents with very low bonds, and the curvature is more around 0 (when parents switch from being net borrower to net saver)
  
    new(β, γ, η, κ, χ, rf, rm, τb, τt, ξ, δ, LTV, PMIlim,PMIprem, ho, price, mb, ms, rmp, rmk,ihini)
  end

end

mutable struct Polp
  c :: Array{Float64,8}
  t :: Array{Float64,8}
  h :: Array{Float64,8}
  b :: Array{Float64,8}
  disc ::  Array{Int,8}


  function Polp(p)
    c = ones(p.nxk,p.nxp,p.nv,p.no,p.no,p.ns,p.nyoung,p.nh)
    t = zeros(p.nxk,p.nxp,p.nv,p.no,p.no,p.ns,p.nyoung,p.nh)
    h = ones(Int64,p.nxk,p.nxp,p.nv,p.no,p.no,p.ns,p.nyoung,p.nh)
    b = zeros(p.nxk,p.nxp,p.nv,p.no,p.no,p.ns,p.nyoung,p.nh)

    disc = ones(p.nxk,p.nxp,p.nv,p.no,p.no,p.ns,p.nyoung,1:p.nh) 
    new(c,t,h,b,disc)
  end
end

mutable struct Polk
  c :: Array{Float64,8}
  h :: Array{Float64,8}
  b :: Array{Float64,8}
  disc ::  Array{Int,8}

  function Polk(p)
    c = ones(p.nxk,p.nbp,p.nh,p.nv,p.no,p.ns,p.nyoung,p.nh)
    h = ones(Int64,p.nxk,p.nbp,p.nh,p.nv,p.no,p.ns,p.nyoung,p.nh) 
    b = zeros(p.nxk,p.nbp,p.nh,p.nv,p.no,p.ns,p.nyoung,p.nh)
    disc = ones(p.nxk,p.nbp,p.nh,p.nv,p.no,p.ns,p.nyoung,p.nh)

    new(c,h,b,disc)
  end
end

##
##
mutable struct Pol
  gk :: Polk
  gp :: Polp

  function Pol(p)
    gk = Polk(p)
    gp = Polp(p)
    new(gk,gp)
  end
end

mutable struct Val
  Vk :: Array{Float64,7}
  Vp :: Array{Float64,7}
  Vk_disc :: Array{Float64,8}
  Vp_disc :: Array{Float64,8}

  function Val(p)
    Vk = zeros(p.nxk,p.nbp,p.nh,p.nv,p.no,p.ns,p.nyoung)
    Vp = zeros(p.nxk,p.nxp,p.nv,p.no,p.no,p.ns,p.nyoung)

    Vk_disc = zeros(p.nxk,p.nbp,p.nh,p.nv,p.no,p.ns,p.nyoung,p.nh)
    Vp_disc = zeros(p.nxk,p.nxp,p.nv,p.no,p.no,p.ns,p.nyoung,p.nh)
    new(Vk, Vp, Vk_disc, Vp_disc)
  end
end

mutable struct NumDet
  vdists :: Array{Float64,1}
  simdists :: Array{Float64,1}
  dist  :: Float64
  cnvrgd :: Bool
  time  :: Float64
  smmobj :: Float64
  simdist :: Float64
  previter :: Int64
  cyclical :: Bool
  cycle_length :: Int64

  function NumDet(p::NumPar)
    vdists = fill(NaN64,100)
    simdists = fill(NaN64,100)
    simdists .= NaN
    dist = 0.0
    cnvrgd = false
    time = 0
    smmobj = -Inf
    simdist = NaN
    previter = 0
    cyclical = false
    cycle_length = 0

    new(vdists,simdists,dist,cnvrgd,time,smmobj,simdist,previter,cyclical,cycle_length)
  end
end


mutable struct Model
  V :: Val
  g :: Pol
  mp :: ModPar
  np :: NumPar
  numdet  :: NumDet
  simpan :: DataFrame
  moms :: DataFrame
  
  function Model(mp::ModPar,np::NumPar)
    g = Pol(np)
    V = Val(np)
    g =     init_sol(np)
    numdet = NumDet(np)

    simpan = DataFrame()
    moms = DataFrame()
    
    new(V,g,mp,np,numdet,simpan,moms)

  end
end

##
