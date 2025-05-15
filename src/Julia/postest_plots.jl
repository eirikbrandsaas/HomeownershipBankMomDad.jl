## Perfrom the "Chetty-Experiment" in the simulated panel
function ChettyEventStudy(M,wealthypar_lim::Float64;store=false)
    @assert wealthypar_lim >= 0.0
    @assert wealthypar_lim <= 1.0
    np = M.np
    mp = M.mp
    df = copy(M.simpan)
    df.ownk = df.ihk.>np.nhr

    sort!(df,[:id,:iak])
    df.hcons = np.h_grd[df.ihk].*np.p_grd[df.is]*mp.κ.*((1 .-df.ownk) .+  df.ownk)
    df.hcons = log.(df.hcons) # Log variables to get log growth rates
    select!(df,[:iak,:id,:ivk,:hcons,:xp,:ihk,:is,:ownk])
    lag!(df,:id,:iak,:ivk,name="Livk")
    lead!(df,:id,:iak,:ivk,name="Fivk")
    diff!(df,:id,:iak,:hcons,name="Dhousecons")


    # Find income status
    # Cut first and last obs (need to observe lags and forwards)
    df = df[(df.iak .>=2) .& (np.age_grd[df.iak] .<= 44),:]
    df[!,:unemp] .= false
    df[(df.ivk .== 1) .& (df.Livk .>= 2) .& (df.Fivk .!= 1),:unemp] .= true # Unemployment if you lose income
    df = leftjoin(df,combine(groupby(df,:id),:unemp => sum),on=:id) # Keeps those with only one spell
    df = df[df.unemp_sum .== 1,:]

    # Find first "unemployment" spell
    df=leftjoin(df,combine(groupby(df[df.unemp.==true,:],:id),:iak => first => :first_unemp ),on=:id)
    df.relyear = (df.iak - df.first_unemp)*2

    df.wpq = (df.xp .>= quantile(df.xp,wealthypar_lim)) # Wealthy if belong to top quantile
    df=df[abs.(df.relyear) .<= 4,:] # Only keep +/- 4 years

    # Agregate
    dfagg = combine(groupby(df,[:relyear,:wpq]), :Dhousecons => mean, :Dhousecons => sem, :Dhousecons => length)
    sort!(dfagg,[:wpq, :relyear])
    dfagg.σ = dfagg.Dhousecons_sem *1.96

    # Plot settings
    upscale = 80
    ylim = (-0.2,0.2)
    dcolor = :orange
    scolor = :black
    dfill = 0.2
    dfd = DataFrame(CSV.File("data/PSID/eventstudy.csv"))
    rename!(dfd,[:wealthy_prnt, :yearrelunemp, :mean_Dlhousing] .=> [:wpq, :relyear, :Dhousecons_mean])
    dfd.σ = dfd.sem_Dlhousing
    hline([0],style=:solid,color=:pink,label="",alpha=1)
    plot!(dfd[dfd.wpq.==false,:relyear],dfd[dfd.wpq.==false,:Dhousecons_mean],label="Data",
        ribbon = dfd[dfd.wpq.==false,:σ],seriescolor=dcolor,fillalpha=dfill)
    plot!(ylim=ylim,title=@sprintf("Parents in bottom %i",wealthypar_lim*100)*"% of Wealth",legend=:bottomleft,legendfontsize=10,xlabel="Years Relative to Unemployment",ylabel="Housing Growth Rate")
    plot!(size=(6*upscale,3.5*upscale))
    pdata1 = deepcopy(plot!())
    plot!(dfagg[dfagg.wpq.==false,:relyear],dfagg[dfagg.wpq.==false,:Dhousecons_mean],label="Model",
        ribbon = dfagg[dfagg.wpq.==false,:σ],seriescolor=scolor,fillalpha=0.25)
    p1 = plot!()

    hline([0],style=:solid,color=:pink,label="",alpha=1)
    plot!(dfd[dfd.wpq.==true,:relyear],dfd[dfd.wpq.==true,:Dhousecons_mean],label="Data",
        ribbon = dfd[dfd.wpq.==true,:σ],seriescolor=dcolor,fillalpha=dfill)
    plot!(ylim=ylim,title=@sprintf("Parents in top %i",(1-wealthypar_lim)*100)*"% of Wealth",legend=false,xlabel="Years Relative to Unemployment",ylabel="Housing Growth Rate")
    plot!(size=(6*upscale,3.5*upscale))
    pdata2 = deepcopy(plot!())
    plot!(dfagg[dfagg.wpq.==true,:relyear],dfagg[dfagg.wpq.==true,:Dhousecons_mean],label="Model",
        ribbon = dfagg[dfagg.wpq.==true,:σ],seriescolor=scolor,fillalpha=0.25)
    p2 = plot!()

    pdata1 = plot(pdata1,legend=false,ylim=ylim)
    pdata2 = plot(pdata2,legend=false,ylim=ylim)
    if store == true
        Plots.pdf(p1,"tabfig/model_housinggrowthpoor_both")
        Plots.pdf(p2,"tabfig/model_housinggrowthrich_both")
        Plots.pdf(pdata1,"tabfig/model_housinggrowthpoor_bothdata")
        Plots.pdf(pdata2,"tabfig/model_housinggrowthrich_bothdata")
    end


    return plot(p1,p2)
end


## Maintain ownership
function find_firsttime_owner!(df)
    # Find first-time owners
    transform!(groupby(df, :id),
    :own => (v -> begin
        idx = findfirst(identity, v)     # first TRUE in that household
        out  = falses(length(v))
        if !isnothing(idx)               # mark only that one row
            out[idx] = true
        end
        out
    end) => :first_own)
end

function maintaining_ownership_plots(M;store=false)
    
    np = M.np
    df = deepcopy(M.simpan)
    @warn "Not obvious what is right here... Should it be: you count as owner the period you own (like I do now) or should it be the period you wake up as owner?"
    df.own = df.ihk.==2

    
    df = df[:,[:id,:iak,:own, :xk, :xp, :ivk, :iop, :iok,:iν]]
    find_firsttime_owner!(df)
    df.wk .= NaN64
    df.wp .= NaN64
    for i = 1:nrow(df)
        iak = df.iak[i]::Int
        ivk = df.ivk[i]::Int
        iν = df.iν[i]::Int
        df.wk[i] = wk(iak,ivk,np)
        df.wp[i] = wp(iak+np.nyoung,iν,np)
    end

    #
    maxlead = 15
    dff = DataFrame(Lead=Int[],first=Float64[],all=Float64[])
    df.age=Int.(np.age_grd[df.iak])
    for F in 0:maxlead
        # 1.  Create a lead‑F version of :own within each household
        leadcol = Symbol("Fown", F)
        transform!(groupby(df, :id), :own => (v -> ShiftedArrays.lag(v, -F)) => leadcol)
    end

    agelo = 25
    agehi = 44
    for F = 0:maxlead
        _mfirst = mean(skipmissing(df[(df.first_own.==true) .&& (agelo .<= df.age .<= agehi),Symbol(:Fown,F)]))
        _mall = mean(skipmissing(df[(df.own.==true)  .&& (agelo .<= df.age .<= agehi),Symbol(:Fown,F)]))
        push!(dff,[Int(F),_mfirst,_mall])
    end



    #
    df2 = DataFrame()
    print("running regressions with leads")
    for F = 1:1:maximum(dff.Lead)-2
        print("$F, ")
        leadcol = Symbol("Fown", F)
        df1 = dropmissing(df, leadcol)
        df1 = df1[df1.xk.>0 .&& df1.xp.>0,:]
        df1.lxk = log.(df1.xk)
        df1.lxp = log.(df1.xp)
        df1.lwk = log.(df1.wk)
        df1.lwp = log.(df1.wp)
        df1.ivk = categorical(df1.ivk)
        df1.iak2 = df1.iak.^2

        mod = lm(@eval(@formula($leadcol ~ lxp + lwp + lxk + lwk + iak + iak2)),df1[df1.first_own.==true  .&& (agelo .<= df1.age .<= agehi),:])
        _idk = findfirst(==("lxk"), coefnames(mod))
        _idp = findfirst(==("lxp"), coefnames(mod))
        _df2 = DataFrame(Lead=F,first_k = coef(mod)[_idk], firstse_k = stderror(mod)[_idk], first_p = coef(mod)[_idp], firstse_p = stderror(mod)[_idp])
        mod = lm(@eval(@formula($leadcol ~ lxp + lwp + lxk + ivk + iak + iak2)),df1[df1.own.==true  .&& (agelo .<= df1.age .<= agehi),:])
        _idk = findfirst(==("lxk"), coefnames(mod))
        _idp = findfirst(==("lxp"), coefnames(mod))
        leftjoin!(_df2,DataFrame(Lead=F,all_k = coef(mod)[_idk], allse_k = stderror(mod)[_idk], all_p = coef(mod)[_idp], allse_p = stderror(mod)[_idp]),on=:Lead)
        df2 = vcat(df2,_df2)
    end

    df2
    tstatadj = 1.68
    dfe = DataFrame(CSV.File("data/PSID/maintain.csv"))
    dff.fv = dff.Lead*2 # In model, period is two years so each lag is two years
    df2.fv = df2.Lead*2 # In model, period is two years so each lag is two years
    sfill = 0.25
    dfill = 0.2
    dcol = :orange
    scol = :black
    upscale = 75
    # ## First-time buyers
    pfo = plot(title="Own Wealth",xlabel="Years Since Purchase",ylabel="Coefficient")
    hline!([0],style=:solid,color=:pink,label="",alpha=1)
    plot!(dfe.fv, dfe.first_k,label="Data",
        ribbon = tstatadj*replace(dfe.firstse_k,missing=>NaN),seriescolor=dcol,fillalpha=dfill,ylims=(-0.05,0.08))
    plot!(df2.fv, [df2.first_k ],label="Model",
        ribbon = (tstatadj*df2.firstse_k),seriescolor=scol,fillalpha=sfill,style=:solid)

    pfp = plot(title="Parental Wealth",xlabel="Years Since Purchase",ylabel="Coefficient")
    hline!([0],style=:solid,color=:pink,label="",alpha=1)
    plot!(dfe.fv, [dfe.first_p ],label="Data",
        ribbon = tstatadj*replace(dfe.firstse_p,missing=>NaN),seriescolor=dcol,fillalpha=dfill,ylims=(-0.05,0.08))
    plot!(df2.fv, [df2.first_p ],label="Model",style=:solid,
        ribbon = (tstatadj*df2.firstse_p),seriescolor=scol,fillalpha=sfill )

    pfm = plot(title="Maintaing ownership",ylims=(0.70,1.0),xlabel="Years (leads)")
    plot!(dfe.fv,dfe.m_first,color=dcol,style=:solid,label="Data")
    plot!(dff.fv,dff.first,color=scol,style=:solid,label="Model")
        
    plot(pfm,pfo,pfp,xlims=(0,10),layout=(1,3),xlabel="Lead (years)", ylabel="Coefficient")

    # ## all owners

    pao = plot(title="Own Wealth",xlabel="Years Since Purchase",ylabel="Coefficient")
    hline!([0],style=:solid,color=:pink,label="",alpha=1)
    plot!(dfe.fv, dfe.all_k,label="Data",
        ribbon = tstatadj*replace(dfe.allse_k,missing=>NaN),seriescolor=dcol,fillalpha=dfill,ylims=(-0.05,0.08))
    plot!(df2.fv, [df2.all_k ],label="Model",
        ribbon = (tstatadj*df2.allse_k),seriescolor=scol,fillalpha=sfill,style=:solid)

    pap = plot(title="Parental Wealth",xlabel="Years Since Purchase",ylabel="Coefficient")
    hline!([0],style=:solid,color=:pink,label="",alpha=1)
    plot!(dfe.fv, [dfe.all_p ],label="Data",
        ribbon = tstatadj*replace(dfe.allse_p,missing=>NaN),seriescolor=dcol,fillalpha=dfill,ylims=(-0.05,0.08))
    plot!(df2.fv, [df2.all_p ],label="Model",style=:solid,
        ribbon = (tstatadj*df2.allse_p),seriescolor=scol,fillalpha=sfill )

    pam = plot(title="Maintaing ownership",ylims=(0.70,1.0),xlabel="Years (leads)")
    plot!(dfe.fv,dfe.m_owner,color=dcol,style=:solid,label="Data")
    plot!(dff.fv,dff.all,color=scol,style=:solid,label="Model")

    plot(pam,pao,pap,layout=(1,3),xlims=(0,12),xlabel="Lead (years)", ylabel="Coefficient")

    

    ##
    p1 = plot(pfo,size=(6*upscale,3.5*upscale),xlims=(0,12))
    p2 = plot(pfp,size=(6*upscale,3.5*upscale),xlims=(0,12))

    if store == true
        Plots.pdf(p1,"tabfig/model_remain_own")
        Plots.pdf(p2,"tabfig/model_remain_parent")
    end


    plot(pfm,pfo,pfp,pam,pao,pap,xlims=(0,12))

end