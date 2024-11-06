function itp_gpage(np::NumPar,gtp,gcp,gdiscp)
    itptp = Array{Interpolations.Extrapolation{Float64,2,Interpolations.GriddedInterpolation{Float64,2,Matrix{Float64},Gridded{Linear{Throw{OnGrid}}},Tuple{Array{Float64,1},Array{Float64,1}}},Gridded{Linear{Throw{OnGrid}}},Flat{Nothing}},6}
    @views gtp_itp = [linear_interpolation((np.xk_grd,np.xp_grd),            gtp[:,:,ivk,ihk,iop,is,iak,ihp],extrapolation_bc=Flat()) for ivk=1:np.nv,ihk=1:np.no,iop=1:np.no,is=1:np.ns,iak=1:np.nyoung,ihp=1:np.nh]
    @views gbp_itp  = [linear_interpolation((np.xk_grd,np.xp_grd),           gcp[:,:,ivk,ihk,iop,is,iak,ihp],extrapolation_bc=Flat()) for ivk=1:np.nv,ihk=1:np.no,iop=1:np.no,is=1:np.ns,iak=1:np.nyoung,ihp=1:np.nh]
    @views gdiscp_itp = [linear_interpolation((np.xk_grd,np.xp_grd),      gdiscp[:,:,ivk,ihk,iop,is,iak,ihp],extrapolation_bc=Flat()) for ivk=1:np.nv,ihk=1:np.no,iop=1:np.no,is=1:np.ns,iak=1:np.nyoung,ihp=1:np.nh]
  
    return gtp_itp,gbp_itp,gdiscp_itp
end
  
  
function itp_gkage(np::NumPar,gck,gdisck)
    t1 = Array{Array{Interpolations.Extrapolation{Float64,2,Interpolations.GriddedInterpolation{Float64,2,Matrix{Float64},Gridded{Linear{Throw{OnGrid}}},Tuple{Array{Float64,1},Array{Float64,1}}},Gridded{Linear{Throw{OnGrid}}},Flat{Nothing}},5},1}
    t2 = Array{Array{Interpolations.Extrapolation{Float64,2,Interpolations.GriddedInterpolation{Float64,2,Matrix{Int64},Gridded{Linear{Throw{OnGrid}}},Tuple{Array{Float64,1},Array{Float64,1}}},Gridded{Linear{Throw{OnGrid}}},Flat{Nothing}},5},1}
    @views gbk_itp = [linear_interpolation((np.xk_grd,np.bp_grd),     gck[:,:,ihp,ivk,iok,is,iak,ihk],extrapolation_bc=Flat()) for ihp=1:np.nh,ivk=1:np.nv,iok=1:np.no,is=1:np.ns,iak=1:np.nyoung, ihk=1:np.nh]
    @views gdisck_itp = [linear_interpolation((np.xk_grd,np.bp_grd),gdisck[:,:,ihp,ivk,iok,is,iak,ihk],extrapolation_bc=Flat()) for ihp=1:np.nh,ivk=1:np.nv,iok=1:np.no,is=1:np.ns,iak=1:np.nyoung, ihk=1:np.nh]
  
    return gbk_itp,gdisck_itp
end

function draw_discretechoice!(np,mp,ivk,iok,iop,is,xk,xp,iak,gdiscp,gtp,gcp,gdisck,probs,prob)
    # This function is used to draw the discrete choices of the kids and parents, based on the "probs" matrix, where:
        # sum of first row = prob of parent choosing to rent,
        # sum of first column = prob of kid choosing to rent
        # The first row, first column, is the prob of both choosing to rent
        # The first row, second column, is the prob of parent choosing to rent, and kid choosing to own etc
    for ihp = 1:np.nh
        probp = gdiscp[ivk,iok,iop,is,iak,ihp](xk,xp)
        hp = np.h_grd[ihp]
        ownp = isowner(ihp,np)
        if probp > 0 # No reason to calculate numbers if probp = 0.0!
            tp = gtp[ivk,iok,iop,is,iak,ihp](xk,xp)*(1.0 - mp.τt )
            cp = gcp[ivk,iok,iop,is,iak,ihp](xk,xp)
            bp = BCp_netshock(xp,cp,tp,hp,ownp,iak+np.nyoung,iop,np.p_grd[is],np,mp) 
            xkcash = xk + tp

            probk = 0.0 # Prob that kid choses housing index ihk
            # Loop over kids housing choices
            for ihk = 1:np.nh
                probk = gdisck[ihp,ivk,iok,is,iak,ihk](xkcash,bp)

                if probk > 0
                    probs[ihp,ihk] = probp*probk
                    
                end
            end
        end
       
    end

    indx = quick_findoutcomePDF(prob,probs)
    ihp,ihk = Tuple(CartesianIndices((np.nh,np.nh))[indx])    

    return ihp,ihk
end


function simulate_oneperiod(np,mp,pr,xp,xk,ivk,iok,iop,is,iak,ind,gen,gdiscp,gtp,gcp,gdisck,gck,_probs)
    iap = iak + np.nyoung
    _probs.=0.0
    ihp,ihk = draw_discretechoice!(np,mp,ivk,iok,iop,is,xk,xp,iak,gdiscp,gtp,gcp,gdisck,_probs,pr.discbth[iak,ind,gen])

    # Note that in this function you already find tp,cp,bp,xkcash (for each discrete combo), so you probably don't have to refind them here below...
    hp = np.h_grd[ihp]
    ownp = isowner(ihp,np)
    hk = np.h_grd[ihk]
    ownk = isowner(ihk,np)

    tp = gtp[ivk,iok,iop,is,iak,ihp](xk,xp)*(1.0 - mp.τt)
    cp = gcp[ivk,iok,iop,is,iak,ihp](xk,xp)
    bp = BCp_netshock(xp,cp,tp,hp,ownp,iap,iop,np.p_grd[is],np,mp) 
    cp, bp =  enforce_feasibleborrowing(cp,bp,ihp,is,iak,mp,np)

    xkcash = xk + tp
    ck = gck[ihp,ivk,iok,is,iak,ihk](xkcash,bp)
    bk = BCk(xkcash,ck,0.0,hk,ownk,iak,iok,np.p_grd[is],np,mp,ivk)
    ck,bk = enforce_feasibleborrowing(ck,bk,ihk,is,iak,mp,np)

    # Next period wealth (note that they depend on the stochasitc (ivk,is,iν)

    # Need to draw these stochastic guys

    isn = findoutcomeCDF(pr.s[iak+1,ind,gen],@view np.ΠsCDF[is,:])
    iνn = findoutcomeCDF(pr.ν[iak+1,ind,gen],@view np.ΠνCDF[iap,:])
    iδp = findoutcomeCDF(pr.δp[iak+1,ind,gen],np.ΠδCDF)
    iδk = findoutcomeCDF(pr.δk[iak+1,ind,gen],np.ΠδCDF)
    ivkn = findoutcomeCDF(pr.v[iak+1,ind,gen],@view np.ΠvCDF[iak,ivk,:])
    if iak < np.nyoung
        xpn = xp′_scalar(bp,hp,ownp,iap,np,mp,isn,iδp,iνn)
        xkn = xk′_scalar(bk,hk,ownk,iak,np,mp,isn,iδk)    
    else
        xpn = xpdead′_scalar(bp,hp,ownp,iap,np,mp,isn,iδp)*(1.0 - mp.τb)
        xkn = xkdead′_scalar(bk,hk,ownk,iak,np,mp,isn,iδk,iνn)
    end

    return xpn,xkn,ivkn,ihk,ihp,iνn,isn,iδp, iδk,tp,cp,bp,ck,bk
end






function initial_states(np,mp)
        ivk = div(np.nv,2,RoundUp)
        iok = mp.ihini
        iak = 1
        iop = 1
        iν = 1
        is = 1
        xk = mean(np.xk_grd)
        xp = mean(np.xp_grd)

        iap = np.nyoung+1
        xp = xp + wp_netshock(iap,np) - wp(iap,iν,np)    

    return ivk, iok, iak, iop, iν, is, xk, xp    
end


function simulate_dynasty!(np,mp,simprob,ind,gdiscp,gtp,gcp,gdisck,gck,ΠΘCDF,
    xpmat, xkmat, ivkmat, iokmat, iopmat, ihkmat, ihpmat, iνmat, ismat, iδpmat, iδkmat, tpmat, cpmat, bpmat, ckmat, bkmat, iakmat, genmat, idmat,
    _probs,_CDF)

    ivk, iok, iop, is, xk, xp = initial_states(np,mp)

    idmat .= ind

    for gen = 1:np.Ngen
        genmat[:,gen] .= gen
        
        ivkmat[1,gen] = ivk
        iokmat[1,gen] = iok
        iopmat[1,gen] = iop
        ismat[1,gen] = is
        xkmat[1,gen] = xk
        xpmat[1,gen] = xp

        if gen == 1 
            ivkmat[1,gen] = ivk = 1
            xkmat[1,gen] = xk = 100.0
            xpmat[1,gen] = xp = 1000.0
        end

        for iak = 1:np.nyoung
            
            # display([gen,gen,iak,xk])
            iakmat[iak,gen] = iak
            xp, xk, ivk, ihk, ihp, iν, is, iδp, iδk, tp, cp, bp, ck, bk =  simulate_oneperiod(np,mp,simprob,xp,xk,ivk,iok,iop,is,iak,ind,gen,gdiscp,gtp,gcp,gdisck,gck,_probs)

            # choices 
            tpmat[iak,gen] = tp
            cpmat[iak,gen] = cp
            bpmat[iak,gen] = bp
            ckmat[iak,gen] = ck
            bkmat[iak,gen] = bk
            ihkmat[iak,gen] = ihk
            ihpmat[iak,gen] = ihp

            # Shocks to next-period wealth 
            iνmat[iak,gen] = iν
            ismat[iak,gen] = is
            iδpmat[iak,gen] = iδp
            iδkmat[iak,gen] = iδk
            # next-period states
            if iak < np.nyoung
                    # update housing choice to state
                iokmat[iak+1,gen] = iok = ihk
                iopmat[iak+1,gen] = iop = ihp
                xpmat[iak+1,gen] = xp
                xkmat[iak+1,gen] = xk
                ivkmat[iak+1,gen] = ivk
            else
                xp,iop,ixn,iok,ivk = nextgen_state(mp,np,ΠΘCDF,xk,xp,ihk,xkmat[iak,gen],ivk,simprob.Θ[ind,gen],_CDF)
                xk = np.θx_grd[ixn,ivk]                
            end
        end       
    end

end



function simpar(M)
    np = M.np
    @views ΠΘCDF = [linear_interpolation((np.xk_grd,),M.np.ΠΘCDF[ixn,ivn,:,ivk],extrapolation_bc=Flat()) for ixn = 1:np.Θnxk, ivn = 1:np.nv,ivk=1:np.nv]
    return ΠΘCDF
end


function simulate(M;gs)
    np = M.np
    mp = M.mp
    g = M.g
    numdet= M.numdet


    ## If you have more than one cycle, you want to ensure that different cycles have the same shocks, so the only difference is the policy functions!
    ## So we will simulate each dynasty cycle_length times.
    _Ndyn_old = np.Ndyn
    rng = StableRNG(4006381) # Barcode on a marker
    _simprob =  simprobs(rng,np.nyoung,np.Ndyn,np.Ngen)
    np.Ndyn = np.Ndyn*numdet.cycle_length
    simprob =  simprobs(rng,np.nyoung,np.Ndyn,np.Ngen)
    for cycle = 1:numdet.cycle_length # Loop over cycles
        offset = _Ndyn_old * (cycle - 1) # Pushes the individual one "cycle" out. (e.g., if Ndyn = 2, and you have three cycles, you want idx  = 1, 4, and 7 for the first individual, 2,5, 8, for the second, 3,6,9 for the third and last)
        for ind = 1:_Ndyn_old
            idx = ind + offset  # Compute index once and reuse it
            @views simprob.ν[:, idx, :] = _simprob.ν[:, ind, :]
            @views simprob.v[:, idx, :] = _simprob.v[:, ind, :]
            @views simprob.discbth[:, idx, :] = _simprob.discbth[:, ind, :]
            @views simprob.Θ[idx, :] = _simprob.Θ[ind, :]
            @views simprob.δp[:, idx, :] = _simprob.δp[:, ind, :]
            @views simprob.δk[:, idx, :] = _simprob.δk[:, ind, :]
            @views simprob.s[:, idx, :] = _simprob.s[:, ind, :]
        end
    end

    # Preallocate
    ΠΘCDF = simpar(M)
    xpmat, xkmat, ivkmat, iokmat, iopmat, ihkmat, ihpmat, iνmat, ismat, iδpmat, iδkmat, tpmat, cpmat, bpmat, ckmat, bkmat, iakmat, genmat, idmat,cycmat = preallocate_simmats(np)
    _probs = fill(0.0,(np.nh,np.nh,np.Ndyn)) # pre-allocate PDF for housing choices
    _CDF = fill(0.0,(np.Θnxk,np.nv,np.Ndyn)) # Pre-allocate CDF for initial wealth and productivity

    # Draw which cycle a dynasty belongs
    dyn_cycle = floor.(Int, range(start = 1, stop = numdet.cycle_length+1, length = np.Ndyn)) # Note: Final element is too high
    dyn_cycle[end] = numdet.cycle_length # Correct final element
    cycleinds = [findall(x->x==cycle,dyn_cycle) for cycle = 1:numdet.cycle_length]
 
    for cycle = 1:numdet.cycle_length

        g = gs[end-cycle+1]    
        gtp,gcp,gdiscp = itp_gpage(np,g.gp.t,g.gp.c,g.gp.disc);
        gck,gdisck = itp_gkage(np,g.gk.c,g.gk.disc);

        Threads.@threads for ind in cycleinds[cycle]
            @views simulate_dynasty!(np,mp,simprob,ind,gdiscp,gtp,gcp,gdisck,gck,ΠΘCDF,
                xpmat[:,:,ind], xkmat[:,:,ind], ivkmat[:,:,ind], iokmat[:,:,ind], iopmat[:,:,ind], ihkmat[:,:,ind], ihpmat[:,:,ind], iνmat[:,:,ind], ismat[:,:,ind], iδpmat[:,:,ind], iδkmat[:,:,ind], tpmat[:,:,ind], cpmat[:,:,ind], bpmat[:,:,ind], ckmat[:,:,ind], bkmat[:,:,ind], iakmat[:,:,ind], genmat[:,:,ind], idmat[:,:,ind],
                _probs[:,:,ind],_CDF[:,:,ind])
            @views cycmat[:,:,ind] .= cycle
        end

        
        for gen = 2:np.Ngen
            _do = mean(xkmat[:,gen,:];dims=2)
            _dn = mean(xkmat[:,gen-1,:];dims=2)
            numdet.simdists[gen] = mean((abs.(_dn - _do)./_do))/2
        end
    end
    
    np.Ndyn = _Ndyn_old # Reset counter
    println("Average percentage difference between generations ",round.(numdet.simdists[2:np.Ngen];digits=5))

    # Create DataFrame of last generation
    pan = calc_ind_pan_lastgen(xpmat, xkmat, ivkmat, iokmat, iopmat, ihkmat, ihpmat, iνmat, ismat, iδpmat, iδkmat, tpmat, cpmat, bpmat, ckmat, bkmat, iakmat, idmat, genmat, cycmat)

    return pan
end


function nextgen_state(mp::ModPar,np::NumPar,ΠΘCDF,xk,xp,ihk,xkstart,ivk,prΘ,_ΠΘCDF)
    xpn = xk + xp
    iop = ihk
    iok = mp.ihini

    for ixn=1:np.Θnxk, ivn = 1:np.nv
        _ΠΘCDF[ixn,ivn] = ΠΘCDF[ixn,ivn,ivk](xkstart)
    end
    index = findoutcomeCDF(prΘ,_ΠΘCDF)
    ixn, ivn = Tuple(CartesianIndices((np.Θnxk,np.nv))[index]) # Find the correct starting wealth and productvity of new child

    return xpn,iop,ixn,iok,ivn
end

function preallocate_simmats(np::NumPar)
    xpmat = fill(-Inf64,np.nyoung,np.Ngen,np.Ndyn)
    xkmat = fill(-Inf64,np.nyoung,np.Ngen,np.Ndyn)
    ivkmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    iokmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    iopmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    ihkmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    ihpmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    iνmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    ismat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    iδpmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    iδkmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    tpmat = fill(-Inf64,np.nyoung,np.Ngen,np.Ndyn)
    cpmat = fill(-Inf64,np.nyoung,np.Ngen,np.Ndyn)
    bpmat = fill(-Inf64,np.nyoung,np.Ngen,np.Ndyn)
    ckmat = fill(-Inf64,np.nyoung,np.Ngen,np.Ndyn)
    bkmat = fill(-Inf64,np.nyoung,np.Ngen,np.Ndyn)
    iakmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    genmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    idmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)
    cycmat = fill(0,np.nyoung,np.Ngen,np.Ndyn)

    return xpmat, xkmat, ivkmat, iokmat, iopmat, ihkmat, ihpmat, iνmat, ismat, iδpmat, iδkmat, tpmat, cpmat, bpmat, ckmat, bkmat, iakmat, genmat, idmat, cycmat
end

function calc_ind_pan_lastgen(xpmat, xkmat, ivkmat, iokmat, iopmat, ihkmat, ihpmat, iνmat, ismat, iδpmat, iδkmat, tpmat, cpmat, bpmat, ckmat, bkmat, iakmat,idmat,genmat, cycmat)
    ind = size(idmat,2)
    pan = (
        id = vec(view(idmat,:,ind,:)),
        gen = vec(view(genmat,:,ind,:)),
        xp = vec(view(xpmat,:,ind,:)),
        xk = vec(view(xkmat,:,ind,:)),
        ivk = vec(view(ivkmat,:,ind,:)),
        iok = vec(view(iokmat,:,ind,:)),
        iop = vec(view(iopmat,:,ind,:)),
        ihk = vec(view(ihkmat,:,ind,:)),
        ihp = vec(view(ihpmat,:,ind,:)),
        iν = vec(view(iνmat,:,ind,:)),
        is = vec(view(ismat,:,ind,:)),
        iδp = vec(view(iδpmat,:,ind,:)),
        iδk = vec(view(iδkmat,:,ind,:)),
        tp = vec(view(tpmat,:,ind,:)),
        cp = vec(view(cpmat,:,ind,:)),
        bp = vec(view(bpmat,:,ind,:)),
        ck = vec(view(ckmat,:,ind,:)),
        bk = vec(view(bkmat,:,ind,:)),
        iak = vec(view(iakmat,:,ind,:)),
        cyc = vec(view(cycmat,:,ind,:))
       )
end

# Simulates and individual when you feed in a state
function simpanel_addvars(np,mp,pan)
    
    nrows=length(pan.id)# np.Ndyn*np.nyoung# nrow(df)::Int

    wkv = zeros(nrows)
    wpv = zeros(nrows)
    hexpk = zeros(nrows)
    utilk = zeros(nrows)
    utilp = zeros(nrows)
    ownk = falses(nrows)
    ownp = falses(nrows)
    hkv = zeros(nrows)
    hpv = zeros(nrows)
    price = zeros(nrows)
    rentexpense = zeros(nrows)
    rent2inc = zeros(nrows)
    h2w = zeros(nrows)
    savingsp = zeros(nrows)
    BorrConsYoung = zeros(nrows)
    cashk = zeros(nrows)
    liq = zeros(nrows)
    mortg = zeros(nrows)
    tp2xp = zeros(nrows)
    tp2xk = zeros(nrows)
    tp2inck = zeros(nrows)
    for i = 1:nrows
        iak = pan.iak[i]::Int
        iok = pan.iok[i]::Int
        ihk = pan.ihk[i]::Int
        ihp = pan.ihp[i]::Int
        is = pan.is[i]::Int
        ivk = pan.ivk[i]::Int
        iν = pan.iν[i]::Int
        xk = pan.xk[i]
        
        
        hkv[i] = hk = np.h_grd[ihk]::Float64
        hpv[i] = hp = np.h_grd[ihp]::Float64
        ownk[i] = _ownk = isowner(ihk,np)::Bool
        ownp[i] = _ownp = isowner(ihp,np)::Bool
        price[i] = pp = np.p_grd[is] :: Float64
        
        ck = pan.ck[i]::Float64
        cp = pan.cp[i]::Float64

        wkv[i] = _wk = wk(iak,ivk,np)
        wpv[i] = wp(iak+np.nyoung,iν,np)
        hexpk[i] = hexp_func(hk,_ownk,iak,iok,pp,np,mp)
        utilk[i] = _utilk = mp.β^(iak-1)*util(ck,hk,_ownk,mp)
        utilp[i] = mp.β^(iak-1)*util(cp,hp,_ownp,mp) + mp.η*_utilk

        rentexpense[i] = _rentexpense = mp.κ*pp.*hk.*(1.0 .- _ownk)
        rent2inc[i] = _rentexpense./_wk
        h2w[i] = pp.*hk.*_ownk ./xk
        savingsp[i] = pan.bp[i] + pp.*hk.*_ownp

        BorrConsYoung[i] = BorrCons(hk,_ownk,pp,iak,mp)
        cashk[i] = xk + pan.tp[i]
        liq[i] = pan.bk[i].*(pan.bk[i] .> 0)
        mortg[i] = pan.bk[i].*(pan.bk[i] .< 0)
        tp2xp[i] = pan.tp[i] ./ pan.xp[i]
        tp2xk[i] = pan.tp[i] ./ (xk)
        tp2inck[i] = pan.tp[i] ./ (_wk)
    end
 
    rentk = (1.0 .- ownk)*mp.κ.*np.h_grd[pan.ihk].*np.p_grd[pan.is]
    rcver = pan.tp .> 0
    pan = merge(pan,(
        wk = wkv,
        wp = wpv,
        hexpk = hexpk,
        utilk = utilk,
        utilp = utilp,
        ownk = ownk,
        ownp = ownp,
        hk = hkv,
        hp = hpv,
        price = price,
        rentk = rentk,
        rcver = rcver,
        rentexpense = rentexpense,
        rent2inc = rent2inc,
        h2w = h2w,
        savingsp = savingsp,
        BorrConsYoung = BorrConsYoung,
        cashk = cashk,
        liq = liq,
        mortg = mortg,
        tp2xp = tp2xp,
        tp2xk = tp2xk,
        tp2inck = tp2inck,
        age = np.age_grd[pan.iak]))

    return pan
end

