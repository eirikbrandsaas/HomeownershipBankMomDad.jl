## Find MPCs
function find_choices(ivk,iok,iop,is,iak,xk,xp,np,mp,gdiscp,gtp,gcp,gdisck,gck)
    ihp = Int(gdiscp[ivk,iok,iop,is,iak,1](xk,xp) <= gdiscp[ivk,iok,iop,is,iak,2](xk,xp) ) + 1
    hp = np.h_grd[ihp]
    ownp = isowner(ihp,np)
    tp = gtp[ivk,iok,iop,is,iak,ihp](xk,xp)*(1.0 - mp.τt )
    cp = gcp[ivk,iok,iop,is,iak,ihp](xk,xp)
    bp = BCp_netshock(xp,cp,tp,hp,ownp,iak+np.nyoung,iop,np.p_grd[is],np,mp) 
    cp, bp =  enforce_feasibleborrowing(cp,bp,ihp,is,iak,mp,np)
    
    xkcash = xk + tp
    ihk = Int(gdisck[ihp,ivk,iok,is,iak,1](xkcash,bp) <gdisck[ihp,ivk,iok,is,iak,2](xkcash,bp)) + 1
    hk = np.h_grd[ihk]
    ownk = isowner(ihk,np)

    ##
    
    
    xkcash = xk + tp
    ck = gck[ihp,ivk,iok,is,iak,ihk](xkcash,bp)
    bk = BCk(xkcash,ck,0.0,hk,ownk,iak,iok,np.p_grd[is],np,mp,ivk)
    ck,bk = enforce_feasibleborrowing(ck,bk,ihk,is,iak,mp,np)
    
    return ihp,hp,cp,bp,tp,ihk,hk,ck,bk
end

function find_mpcs(M,Δ)
    np = M.np
    mp = M.mp
    df = deepcopy(M.simpan)

    gtp,gcp,gdiscp = itp_gpage(np,M.g.gp.t,M.g.gp.c,M.g.gp.disc);
    gck,gdisck = itp_gkage(np,M.g.gk.c,M.g.gk.disc);


    mpc = fill(NaN64,nrow(df))
    mpt = fill(NaN64,nrow(df))
    
    for i = 1:nrow(df)
        dfr = df[i,:]
        ivk = dfr.ivk
        iok = dfr.iok
        iop = dfr.iop
        is = dfr.is
        iak = dfr.iak

        xk = dfr.xk
        xp = dfr.xp

        # Note: The find_choices function is based on the simulate_oneperiod() and draw_discretechoice!()
        ihp1,hp1,cp1,bp1,tp1,ihk1,hk1,ck1,bk1, = find_choices(ivk,iok,iop,is,iak,xk  ,xp,np,mp,gdiscp,gtp,gcp,gdisck,gck)

        ihp2,hp2,cp2,bp2,tp2,ihk2,hk2,ck2,bk2, = find_choices(ivk,iok,iop,is,iak,xk+Δ,xp,np,mp,gdiscp,gtp,gcp,gdisck,gck)
        
        mpc[i] = (ck2-ck1)/Δ
        mpt[i] = (tp2-tp1)/Δ

    end

    h2mnw, h2mliq = identify_h2m(df,np)

    df.mpc = mpc
    df.mpt = mpt
    df.h2mnw = h2mnw
    df.h2mliq = h2mliq
    return df
end

function identify_h2m(df,np)
    threshold_nw = 2/12*0.5
    threshold_liq = 0.165*0.5


    _wk = fill(NaN64,nrow(df))
    for i = 1:nrow(df)
        iak = df.iak[i]::Int
        ivk = df.ivk[i]::Int
        _wk[i] = wk(iak,ivk,np)
    end

    h2mnw = df.xk.<threshold_nw*_wk
    h2mliq = (abs.(df.bk).*(df.bk.<0).<threshold_liq*_wk)

    return h2mnw, h2mliq
end


function create_MPC_table(Mb,Mn,store=false)

    df = DataFrame() # Table to hold results
    agelo = 25
    agehi = 44


    iaklo = argmin(abs.(NumPar().age_grd .- agelo))
    iakhi = argmin(abs.(NumPar().age_grd .- agehi))
    Δ = +5.0
    
    ##
    # MPC with altruism
    M2 = deepcopy(Mn);
    M2.simpan = deepcopy(Mb.simpan)
    df2 = find_mpcs(M2,Δ)
    df2 = df2[iaklo.<= df2.iak .<= iakhi,:]
    append!(df,DataFrame(
        type = "No Altruism",
        median = median(df2.mpc),
        mean = mean(df2.mpc),
        h2mliqavgy = mean(df2[df2.h2mliq.==1,:mpc]),
        h2mliqavgn = mean(df2[df2.h2mliq.==0,:mpc]),
        h2mnwavgy = mean(df2[df2.h2mnw.==1,:mpc]),
        h2mnwavgn = mean(df2[df2.h2mnw.==0,:mpc]),
        wh2mnwavg = mean(df2[df2.h2mnw.==0 .&& df2.h2mliq.==1,:mpc]),
        )    
    )

    # Households (unexpectdly) doesn't reciver a transfer this period (set transfer tax = 1.0). You get almost identical results if you instead use the no-altruism policy functions of the parent
    M3 = deepcopy(Mb);
    M3.mp.τt = 1
    df3 = find_mpcs(M3,Δ)
    df3 = df3[iaklo.<= df3.iak .<= iakhi,:]
    append!(df,DataFrame(
        type = "\\;Future transfers only",
        median = median(df3.mpc),
        mean = mean(df3.mpc),
        h2mliqavgy = mean(df3[df3.h2mliq.==1,:mpc]),
        h2mliqavgn = mean(df3[df3.h2mliq.==0,:mpc]),
        h2mnwavgy = mean(df3[df3.h2mnw.==1,:mpc]),
        h2mnwavgn = mean(df3[df3.h2mnw.==0,:mpc]),
        wh2mnwavg = mean(df3[df3.h2mnw.==0 .&& df3.h2mliq.==1,:mpc]),
        )    
    )
    
    # Household thinks this is the last transfers (i.e., use kid functions without altruism and parent function with altruism)
    M4 = deepcopy(Mb);
    M4.g.gk = deepcopy(Mn.g.gk)
    df4 = find_mpcs(M4,Δ)
    df4 = df4[iaklo.<= df4.iak .<= iakhi,:]
    append!(df,DataFrame(
        type = "\\;Current transfer only",
        median = median(df4.mpc),
        mean = mean(df4.mpc),
        h2mliqavgy = mean(df4[df4.h2mliq.==1,:mpc]),
        h2mliqavgn = mean(df4[df4.h2mliq.==0,:mpc]),
        h2mnwavgy = mean(df4[df4.h2mnw.==1,:mpc]),
        h2mnwavgn = mean(df4[df4.h2mnw.==0,:mpc]),
        wh2mnwavg = mean(df4[df4.h2mnw.==0 .&& df4.h2mliq.==1,:mpc]),
        )    
    )



    # MPC with altruism
    dfb = find_mpcs(Mb,Δ)
    dfb = dfb[iaklo.<= dfb.iak .<= iakhi,:]
    append!(df,DataFrame(
        type = "Altruism",
        median = median(dfb.mpc),
        mean = mean(dfb.mpc),
        h2mliqavgy = mean(dfb[dfb.h2mliq.==1,:mpc]),
        h2mliqavgn = mean(dfb[dfb.h2mliq.==0,:mpc]),
        h2mnwavgy = mean(dfb[dfb.h2mnw.==1,:mpc]),
        h2mnwavgn = mean(dfb[dfb.h2mnw.==0,:mpc]),
        wh2mnwavg = mean(dfb[dfb.h2mnw.==0 .&& dfb.h2mliq.==1,:mpc]),
        )    
    )


    ## Create a nice latex table
    par_latex = latexify(df[!,[:type,:mean,:h2mliqavgy,:h2mliqavgn]]; env=:table, latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true)
    par_latex = replace(par_latex,"cccc" => "lrrr") # Change Alignement
    par_latex = replace(par_latex,"type" => " ") # Change Alignement
    par_latex = replace(par_latex," mean " => " \\multicolumn{1}{c}{Mean} ") # Change Alignement
    par_latex = replace(par_latex," h2mliqavgy " => " \\multicolumn{1}{c}{Liq. Constrained} ") # Change Alignement
    par_latex = replace(par_latex," h2mliqavgn" => " \\multicolumn{1}{c}{Liq. Unconstrained} ") # Change Alignement

    println(par_latex)

    if store == true
            filename = "tabfig/postest/MPCs" 
            write(filename*".tex", par_latex)
    end

end