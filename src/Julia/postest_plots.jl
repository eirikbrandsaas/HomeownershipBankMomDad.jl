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