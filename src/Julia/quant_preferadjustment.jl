function prefer_adjustment(dfout::DataFrame,Vbk_itp,Vnk_itp,Vbp_itp,Vnp_itp,gnp_t_itp,gnp_b_itp,gnp_disc_itp,np::NumPar)
    # To find whether kids prefer adjustment costs you must also update the parent's gift giving choices, which is why this is not straightforward.
    # This also requires some interpolation and so on.
    @assert maximum(dfout.ihk) == 2
    preferk = fill(false,size(dfout)[1])
    preferp = fill(false,size(dfout)[1])
    wtpk = fill( -1.0,size(dfout)[1])
    wtpp = fill( -1.0,size(dfout)[1])
    tpncol = fill(0.0,size(dfout)[1])
    for ind = 1:size(dfout)[1]
        df = dfout[ind,:]
        xk = df.xk
        xp = df.xp
        bp = df.bp
        ivk = df.ivk
        iok = df.iok
        iop = df.iop
        is = df.is
        ihpn = df.ihp
        iak = df.iak
        tp = df.tp

        ## Choices if illiquid
        # Here the choices in the simulated panel is already optimal 
        Vbk = Vbk_itp[ihpn,ivk,iok,is,iak](xk + tp,bp)
        Vbp = Vbp_itp[ivk,iok,iop,is,iak](xk,xp)

        ## Chocies if liquid
        # Here you have to find the new choices, since the simulated choices are suboptimal with liquid housing
        probp = [gnp_disc_itp[ihpnn,ivk,iok,iop,is,iak](xk,xp) for ihpnn = 1:np.nh]
        ihpnn = argmax(probp)
        tpn = gnp_t_itp[ihpnn,ivk,iok,iop,is,iak](xk,xp)
        bpn = gnp_b_itp[ihpnn,ivk,iok,iop,is,iak](xk,xp)
        Vnk = Vnk_itp[ihpnn,ivk,iok,is,iak](xk + tpn,bpn)
        Vnp = Vnp_itp[ivk,iok,iop,is,iak](xk,xp)
        tpncol[ind] = tpn

        if Vbk > Vnk
            preferk[ind] = true
            wtpk[ind] = (Vbk - Vnk)
        end
        if Vbp > Vnp
            preferp[ind] = true
            wtpp[ind] = (Vbp - Vnp)
        end
    end
    dfout.wtpk = wtpk
    dfout.preferk = preferk
    dfout.preferp = preferp
    dfout.tpn = tpncol
    dfout.rcvern =  dfout.tpn.> 0

    return dfout
end

function solve_adjustmentcost(Mb)
    mp = Mb.mp
    np = Mb.np
    msbench = copy(mp.ms)
    mp.ms[1:np.nyoung] .=  0.0
    Mn = Model(mp,np)
    solvemodel2nd!(Mn,np,mp)
    mp.ms = msbench

    return Mn
end
function quant_adjustmentcost(Mb::Model,Mn::Model)
    @assert Mb.np.nh == 2 # COde now only works with two housing sizes
    np = Mb.np


    ## Interpolation objects
    Vbk = [linear_interpolation((np.xk_grd,np.bp_grd),Mb.V.Vk[:,:,ihpn,ivk,iok,is,iak],extrapolation_bc=Flat()) for ihpn=1:np.nh,ivk=1:np.nv,iok=1:np.no,is=1:np.ns,iak=1:np.nyoung]
    Vnk = [linear_interpolation((np.xk_grd,np.bp_grd),Mn.V.Vk[:,:,ihpn,ivk,iok,is,iak],extrapolation_bc=Flat()) for ihpn=1:np.nh,ivk=1:np.nv,iok=1:np.no,is=1:np.ns,iak=1:np.nyoung]
    Vbp = [linear_interpolation((np.xk_grd,np.xp_grd),Mb.V.Vp[:,:,ivk,iok,iop,is,iak],extrapolation_bc=Flat()) for ivk=1:np.nv,iok=1:np.no,iop=1:np.no,is=1:np.ns,iak=1:np.nyoung]
    Vnp = [linear_interpolation((np.xk_grd,np.xp_grd),Mn.V.Vp[:,:,ivk,iok,iop,is,iak],extrapolation_bc=Flat()) for ivk=1:np.nv,iok=1:np.no,iop=1:np.no,is=1:np.ns,iak=1:np.nyoung]
    gnp_disc =[linear_interpolation((np.xk_grd,np.xp_grd),Mn.g.gp.disc[:,:,ivk,iok,iop,is,iak,ihpn],extrapolation_bc=Flat()) for ihpn=1:np.nh,ivk=1:np.nv,iok=1:np.no,iop=1:np.no,is=1:np.ns,iak=1:np.nyoung]
    gnp_b = [linear_interpolation((np.xk_grd,np.xp_grd),Mn.g.gp.b[:,:,ivk,iok,iop,is,iak,ihpn],extrapolation_bc=Flat()) for ihpn=1:np.nh,ivk=1:np.nv,iok=1:np.no,iop=1:np.no,is=1:np.ns,iak=1:np.nyoung]
    gnp_t = [linear_interpolation((np.xk_grd,np.xp_grd),Mn.g.gp.t[:,:,ivk,iok,iop,is,iak,ihpn],extrapolation_bc=Flat()) for ihpn=1:np.nh,ivk=1:np.nv,iok=1:np.no,iop=1:np.no,is=1:np.ns,iak=1:np.nyoung]

    ## Find if they prefer adjustment costs or not
    df = copy(prefer_adjustment(Mb.simpan,Vbk,Vnk,Vbp,Vnp,gnp_t,gnp_b,gnp_disc,np))

    ## Manipulate/calculate the numers you care about

    df.preferboth = (df.preferp + df.preferk) .== 2
    display([mean(df.preferk), mean(df.preferp), mean(df.preferboth)])

    dfini = df[df.iak.==1,:]
    wp_lim = percentile(dfini.xp,75)
    dfini = dfini[dfini.xp .>=wp_lim,[:id]]
    df = leftjoin(df,dfini,on=:id)
    df.age = np.age_grd[df.iak]
    df = df[df.age .<=45,:]
    df.ownk = (df.ihk.>np.nhr)
    df.rcver = (df.tp.>0)
    dfall = copy(df)


    df1,df_all = create_summary_stats(df,:preferk)
    df2,_ = create_summary_stats(df,:preferp)
    df3,_ = create_summary_stats(df,:preferboth)

    if size(df2)[1] == 1
        push!(df2,df2[1,:])
        df2[2,:preferp] = abs(df2[1,:preferp] - 1.0)
    end


    dispvars = [:age_mean, :xk_median,:xp_median,:ownk_mean, :rcver_mean,:tp_mean,:percent,:preferp_mean]
    df_stats= vcat(df1[!,dispvars],df2[!,dispvars],df_all[!,dispvars])
    display(df_stats)   
    display(df1)
    return df, df_stats, dfall
end

function create_summary_stats(df::DataFrame,byvar::Symbol)
    meanvars = [:age,:xk,:xp,:ownk,:rcver,:tp,:rcvern, :tpn,:preferk,:preferp,:preferboth]
    medvars = [:xk, :xp,:preferk,:preferp,:preferboth]
    dfout = combine(groupby(df,[byvar]), meanvars .=> mean)
    dfout = leftjoin(dfout,combine(groupby(df,[byvar]), medvars .=> median),on=byvar)
    dfout = leftjoin(dfout,combine(groupby(df[:,[:age, byvar]],[byvar]), [:age, byvar].=> length),on=byvar)

    rename!(dfout,:age_length => :percent)
    dfout[!,:percent] ./= nrow(df)
    dfout[!,:smpl] .= [string(byvar)]

    agg = combine(df, meanvars .=> mean)
    agg = hcat(agg,combine(df, medvars .=> median))
    agg.percent = [1.00]
    agg[!,:smpl]  .= ["All"]

    return dfout, agg
end

function adjustmentcost_plot(dfall;store=false)
    szfac = (900,400).*0.43            # plt size
    tallsize = (700,800).*0.3

    colors = [:black, :orange]
    pattern = [:solid :dash]
    # Preference for liquidity
    temp = [combine(groupby(dfall,:iak,sort=true),:preferk => func) for func in [mean, sem]] # Creates a dataframe with data by age
    ages = np.age_grd[temp[1].iak]

    plot(ages,temp[1].preferk_mean.*(120),legend=false,ylim=(0,100),color=colors[1],style=pattern[1])

    p1 = plot!(xlabel="Child Age",size=szfac,ylabel="Prefer Illiquid (%)")
    pout = plot!(p1)
    pb = plot!(deepcopy(pout),size=tallsize)
    if store == true
        Plots.pdf(pout,"tabfig/preferliq/prefer")
        Plots.pdf(pb,"tabfig/preferliq/prefer_pres")
        pb = plot!(deepcopy(pb),title = "No",ylabel="%")
        Plots.pdf(pb,"tabfig/preferliq/prefer_pres2")
    end

    
    temp = combine(groupby(dfall,[:preferk,:iak],sort=true),[:xk,:xp, :ownk] .=> mean)
    plot(ages,temp[temp.preferk.==1,:ownk_mean].*100,label="Prefer Illiquid",legend=:bottomright,color=colors[1],style=pattern[1])
    plot!(ages,temp[temp.preferk.==0,:ownk_mean].*100,label="Prefer Liquid",color=colors[2],style=pattern[2],legendfontsize=10)
    p2 = plot!(ylim=(0,100),xlabel="Child Age",size=szfac,ylabel="Homeowners %")
    pout = plot!(p2)
    pb = plot!(deepcopy(pout),size=tallsize)
    if store == true
        Plots.pdf(pout,"tabfig/preferliq/prefer_owner")
        Plots.pdf(pb,"tabfig/preferliq/prefer_owner_pres")
    end

    temp = combine(groupby(dfall,[:preferk,:iak],sort=true),[:xk,:xp, :ownk] .=> mean)
    plot(ages,temp[temp.preferk.==1,:ownk_mean].*100,label="No",legend=:bottomright,color=colors[1],style=pattern[1])
    plot!(ages,temp[temp.preferk.==0,:ownk_mean].*100,label="Yes",color=colors[2],style=pattern[2],legendfontsize=10)
    pb = plot!(ylim=(0,100),xlabel="Child Age",size=szfac)
    pb = plot!(deepcopy(pb),title = "Homeownership",ylabel="%",size=tallsize)
    if store == true
        Plots.pdf(pb,"tabfig/preferliq/prefer_owner_pres2")
    end


    plot(ages,temp[temp.preferk.==1,:xk_mean]./temp[temp.preferk.==0,:xk_mean],label="Child Wealth Ratio",color=colors[1],style=pattern[1])
    plot!(ages,temp[temp.preferk.==1,:xp_mean]./temp[temp.preferk.==0,:xp_mean],label="Parent Wealth Ratio",color=colors[2],style=pattern[2],legendfontsize=10)
    p3 = plot!(ylim=(0,5.0),xlabel="Child Age",size=szfac,legend=:bottomleft,ylabel="Wealth Ratio")
    pout = plot!(p3)
    pb = plot!(deepcopy(pout),size=tallsize)
    if store == true
        Plots.pdf(pout,"tabfig/preferliq/prefer_wealth")
        Plots.pdf(pb,"tabfig/preferliq/prefer_wealth_pres")
    end

    temp = combine(groupby(dfall,[:preferk,:iak],sort=true),[:xk,:xp, :ownk] .=> median)
    plot(ages,temp[temp.preferk.==1,:xp_median],label="No",color=colors[1],style=pattern[1])
    plot!(ages,temp[temp.preferk.==0,:xp_median],label="Yes",color=colors[2],style=pattern[2])
    p4 = plot!(xlabel="Child Age",size=szfac,legend=:topright,title="Parent Wealth",ylabel="1000s USD")
    pout = plot!(p4)
    pb = plot!(deepcopy(pout),size=tallsize)
    if store == true
        Plots.pdf(pb,"tabfig/preferliq/prefer_wealthp_pres")
        plot(ages,temp[temp.preferk.==1,:xp_median],label="Prefer Illiquid",color=colors[1],style=pattern[1])
        plot!(ages,temp[temp.preferk.==0,:xp_median],label="Prefer Liquid",color=colors[2],style=pattern[2],legendfontsize=10)
        p4 = plot!(xlabel="Child Age",size=szfac,legend=:topright,title="Parent Wealth",ylabel="1000s USD")
        Plots.pdf(plot!(p4,title=""),"tabfig/preferliq/prefer_wealthp")
    end

    temp = combine(groupby(dfall,[:preferk,:iak],sort=true),[:xk,:xp, :ownk] .=> mean)
    plot(ages,temp[temp.preferk.==1,:xk_mean],label="No",color=colors[1],style=pattern[1])
    plot!(ages,temp[temp.preferk.==0,:xk_mean],label="Yes",color=colors[2],style=pattern[2],legendfontsize=10)
    p5 = plot!(xlabel="Child Age",size=szfac,legend=:topright,title="Child Wealth",ylabel="1000s USD")
    pout = plot!(p5)
    pb = plot!(deepcopy(pout),size=tallsize)
    if store == true
        Plots.pdf(pb,"tabfig/preferliq/prefer_wealthk_pres")
        plot(ages,temp[temp.preferk.==1,:xk_mean],label="Prefer Illiquid",color=colors[1],style=pattern[1])
        plot!(ages,temp[temp.preferk.==0,:xk_mean],label="Prefer Liquid",color=colors[2],style=pattern[2],legendfontsize=10)
        p5 = plot!(xlabel="Child Age",size=szfac,legend=:topleft,title="",ylabel="1000s USD")

        Plots.pdf(p5,"tabfig/preferliq/prefer_wealthk")
    end

    p4 = plot(p1,p2,p3,p4,p5,size=szfac.*(3,1))
    return p4
end

function adjustmentcost_table(df_in;store=false)
    df2 = copy(df_in)
    select!(df2,[:percent,:age_mean,:xk_median, :xp_median, :ownk_mean, :rcver_mean, :tp_mean, :preferp_mean])
    df = DataFrame(tmp= ["Fraction of Children", "Age","Child Wealth", "Parent Wealth", "Child Ownership Rate", "Transfer Rate", "Transfer Size","Parents Prefering Costs"],
        DislikeCost = collect(df2[1,:]),
        LikeCost = collect(df2[2,:]),
        All = collect(df2[3,:]))

    df_latex = latexify(df; env=:table,latex=false,fmt="%1.2f",escape_underscore=true,booktabs = true,convert_unicode=true,
        head = [" ", " \\multicolumn{1}{c}{Dislike Costs }", "\\multicolumn{1}{c}{Prefer Costs }", "\\multicolumn{1}{c}{All Children}"])
    df_latex = replace(df_latex,"cccc" => "l rrr")
    append_info = @sprintf("m%i_d%i",Dates.month(Dates.today()),Dates.day(Dates.today()))
    filename = "tabfig/preferliq/tab_"*append_info*".tex"
    filename2 = "tabfig/preferliq/tab_mostrecent.tex"
    write(filename, df_latex)
    if store == true
        write(filename2, df_latex)
    end
end
